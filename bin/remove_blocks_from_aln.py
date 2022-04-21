#!/usr/bin/env python
import re, os, getopt, sys

def Usage():
	print 'remove_blocks_from_aln.py Usage:'
	print 'Removes regions from a fasta alignment based on a tab file input. You can choose to remove or keep the regions in the tab file'
	print 'remove_blocks_from_aln.py [options]'
	print 'Options:'
	print '-a|align <file name>\talignment file name'
	print '-o|out <file name>\toutput file name'
	print '-t|tab <file name>\ttab file name (containing regions to keep/remove)'
	print '-c|cut\t\t\toutput alignment with regions cut out (default is to mask them)'
	print '-k|keep\t\t\toutput alignment containing only regions defined in tab file'
	print '-r|reference <name>\treference name (optional, but required if there are gaps in the reference sequence relative to the tab file)'
	print '-R|refrem\t\tDo not remove blocks from reference sequence (default is to remove from all sequences)'
	print '-s|symbol <char>\tSymbol to use for removed regions (default = N)'
	print '-h|help\t\t\tshow this help'

# Define functions
def parse_tab( tabfile ):
    regions = []
    for line in open(tabfile, 'r'):
	    if ".." in line:
		    regex = re.search("(\d+)\.\.(\d+)", line)
		    region = [int(regex.group(1))-1, int(regex.group(2))]
			
		    region.sort()

		    if "complement" in line:
			    region.append('r') # reverse
		    else:
			    region.append('f') # forward
			
		    regions.append(region)
			
    print "Found", len(regions), "regions"
    regions.sort()

    return regions
			
def die( message ):
    sys.stderr.write("!!%s!!\n" % message)
    sys.exit(2)

def getOptions(arg):
	try:
		opts, args = getopt.getopt(argv, "ho:a:t:kr:Rs:c", ["align=", "out=", "=tab", "keep", "reference=", "refrem", "symbol=", "cut"])
	except getopt.GetoptError:
		print "Option Error!", argv
		Usage()
		sys.exit(2)

	outfile=''
	alnfile=''
	tabfile=''
	keepremove='r'
	reference=''
	refrem=True
	symbol="N"

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			Usage()
			sys.exit()
		elif opt in ("-o", "--out"):
			outfile=os.path.abspath(arg)
		elif opt in ("-a", "--align"):
			alnfile=os.path.abspath(arg)
		elif opt in ("-t", "--tab"):
			tabfile=os.path.abspath(arg)
		elif opt in ("-k", "--keep"):
			keepremove='k'
		elif opt in ("-c", "--cut"):
			keepremove='c'
		elif opt in ("-r", "--reference"):
			reference=arg
		elif opt in ("-R", "--refrem"):
			refrem=False
		elif opt in ("-s", "--symbol"):
			symbol=arg

	if alnfile=='' or not os.path.isfile(alnfile):
		die("Error: Alignment file '%s' not found" % alnfile)
	elif tabfile=='' or not os.path.isfile(tabfile):
		die("Error: Tab file '%s' not found" % tabfile)
	elif outfile=='':
		die('Error: No output file specified')
	symbol=symbol.upper()
	if symbol not in ["N", "X", "?", "-"]:
		die('Error: Symbol must be N, X, ? or - ')

	return alnfile, outfile, tabfile, keepremove, reference, refrem, symbol
		
def get_ref_sequence( refname, alnfile ):
    refseq = ''
    capture = False
    for line in open(alnfile):
        if capture:
            if ">" in line:
                break
            else:
                refseq += line.rstrip()
        elif refname in line:
            capture = True
    if refseq != '':
        return refseq
    else:
        die("Reference could not be found in alignment")

def adjust_regions( regions, ref ):
    reftoaln = {}
    refnum = 0
    for alnnum, base in enumerate(ref):
        if base != "-":
            reftoaln[refnum]=alnnum
            refnum+=1
    
    adj_regions = []
    for region in regions:
        try:
            start = reftoaln[region[0]]
            stop = reftoaln[region[1]]
            adj_regions.append([start, stop, region[2]])
        except KeyError:
            die("Region has locations outside of the reference sequence length: %d - %d" % (region[0], region[1]))
    
    return adj_regions

def mask_regions( regions, seq, s_char ):
    for r in regions:
        start = r[0]
        end = r[1]
        diff = end - start
        try:
            seq = "%s%s%s" % (seq[:start], s_char*diff, seq[end:])
        except IndexError:
            die("Region has locations outside of the reference sequence length: %d - %d" % (start, end))
    return seq

def cut_regions( regions, seq ):
    m = mask_regions(regions, seq, "!")
    return re.sub("!", "", m)
    
def keep_regions( regions, seq ):
    keep = ""
    for r in regions:
        start = r[0]
        end = r[1]
        direction = r[2]
        try:
            if direction == 'f':
                keep += seq[start:end]
            elif direction == 'r':
                keep += rev_comp(seq[start:end])
        except IndexError:
            die("Region has locations outside of the reference sequence length: %d - %d" % (start, end))
    return keep

def rev_comp( seq ):
    rc = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    rev_comp = ""
    for base in seq[::-1]:
        try:
            rev_comp += rc[base]
        except KeyError:
            rev_comp += base
    return rev_comp

def validate_aln( alnfile ):
    lens = []
    seq = ""
    for line in open(alnfile):
        if ">" in line:
            if seq != "":
                lens.append(len(seq))
            seq = ""
        else:
            seq += line.strip()
            
    guide = lens[0]
    for l in lens[1:]:
        if l != guide:
            die("Input sequences are not of equal lengths.")
    return guide

def run(alnfile, outfile, tabfile, keepremove, reference, refrem, symbol):
    aln_len = validate_aln(alnfile)
    
    # get regions from the tab file
    regions = parse_tab(tabfile)

    # adjust the regions to the reference if one is given
    if reference != '':
        refseq = get_ref_sequence(reference, alnfile)
        if "-" in refseq:
            regions = adjust_regions(regions, refseq)

    # write new aln to file
    out = open(outfile, 'w')
    current = ""
    process = True
    no_seqs = 0
    new_aln_len = 0
    for line in open(alnfile):
        if ">" in line:
            no_seqs += 1
            if current != "" and process:                
                if keepremove == 'r':
                    new_seq = mask_regions(regions, current, symbol)
                elif keepremove == 'c':
                    new_seq = cut_regions(regions, current)
                elif keepremove == 'k':
                    new_seq = keep_regions(regions, current)
                new_aln_len = len(new_seq)
                out.write(new_seq + "\n")
            elif current != "" and not process:
                out.write(current + "\n")
            
            # check for reference. skip masking if refrem is false
            if ">"+reference == line.strip() and not refrem:
                process = False
            else:
                process = True
                
            out.write(line)
            current = ""
        else:
            current += line.strip()
    # do last sequence
    if keepremove == 'r':
        new_seq = mask_regions(regions, current, symbol)
    elif keepremove == 'c':
        new_seq = cut_regions(regions, current)
    elif keepremove == 'k':
        new_seq = keep_regions(regions, current)
    out.write(new_seq + "\n")
    out.close()
    
    print "Adjusted %d sequences" % no_seqs
    print "Original alignment length: %d\tNew alignment length:%d" % (aln_len, new_aln_len)
    print "Done."

if __name__ == "__main__":
    argv=sys.argv[1:]
    alnfile, outfile, tabfile, keepremove, reference, refrem, symbol = getOptions(argv)
    run(alnfile, outfile, tabfile, keepremove, reference, refrem, symbol)