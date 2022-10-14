import csv, pickle, operator, os, sys         
import collections
from functools import reduce

#   parsing the arguments
args=sys.argv

if len(args)<2:
    pathTBRuns=(os.path.dirname(os.getcwd()))
    instats="*.csv"
    pathPatterns="/Users/andries/Documents/HPC/Pipelines/btb-seq"
    refName="Mycbovis-2122-97_LT708304.fas"
    TBRun ="test"
    qth=8
    thMinGoodCov=2
    thCovProp=0.2
    thqualsnp=150
    thqualnonsnp=0
    strainVCF="*.vcf"
else:    
    pathTBRuns=(os.path.dirname(os.getcwd()))
    instats=sys.argv[1]
    pathPatterns=sys.argv[2]
    refName=sys.argv[3]
    TBRun=sys.argv[4]
    qth=int(sys.argv[5])
    thMinGoodCov=int(sys.argv[6])
    thCovProp=float(sys.argv[7])
    thqualsnp=int(sys.argv[8])
    thqualnonsnp=int(sys.argv[9])
    strainVCF=sys.argv[10]

patternsDetailsFile="CSSnewclusters_LT708304_230119.csv" #"CSSnewclusters_181115.csv" #"patterns20131220.csv" "CSSnewclusters_LT708304_181217.csv"
patternsBritishBTBFile="patternsBritishBTB_LT708304.csv"
patternsPinnipediiFile="patternsPinnipedii_LT708304.csv"
patternsMic_PinFile="patternsMic_Pin_LT708304.csv"
patternsMicrotiFile="patternsMicroti_LT708304.csv"
patternsBTBFile="patternsBTB_LT708304.csv"


# reads a csv file
# return a list where each element is a list containing the element of a row.
def readTable(fname,ch):
    infile=open(fname,"r")
    data = csv.reader(infile, delimiter=ch)
    dataOut = [row for row in data]
    infile.close()
    return dataOut


# writes a list into a csv file
# each element in the list is written as a row in the csv file
def writeCSV(fname,matrix):
    with open(fname, "w") as fileOut:
        writer = csv.writer(fileOut)
        writer.writerows(matrix)
        print("file "+fname+" saved.")

# transposes a python list 
def listT(matrix):
    zipped = zip(*matrix)
    output = [list(item) for item in zipped]
    return output

# compares a strain gSS base calls (strpat) to the a genotype group pattern give by groPat
# strPat is the reference patern.
# it returns two values, a list with the percentage of matches (M), mismatches (MM), notCovered (N) and anomalous (A).
#     The second list is a vector with the same length as strPat with Ms, MMs,Ns and As as corresponding.       

#This part does the matching - finds the closest reference SNP combination on a per sample basis     
def comparePatterns(refPat,strPat,groPat):
    lenPat=len(refPat)
    if lenPat!=len(strPat) or lenPat!=len(groPat):
        print("Different values refPat,strPat,groPat: "+str(lenPat)+"  "+str(len(strPat))+"  "+str(len(groPat)))
    res=[]
    for i in range(lenPat):
        if strPat[i].upper()=="N":
                res.append("N")
        else:
            if strPat[i].upper() not in ['A','C','G','T']:
                res.append("A")
            else:
                if strPat[i].upper()==groPat[i]: 
                    res.append("M")
                else:
                    if strPat[i].upper()==refPat[i] or groPat[i]==refPat[i]:
                        res.append("MM")
                    else:
                        res.append("A")
    counts=collections.Counter(res)
    counts=[round(100*float(counts['M'])/lenPat,2),round(100*float(counts['MM'])/lenPat,2),round(100*float(counts['N'])/lenPat,2),round(100*float(counts['A'])/lenPat,2)]
    return [counts,res]


#this is the key part that runs the per sample stage1 genotyping

def findGenotypeOneSample(strainsDetailsTittle,strainDetails,pathTBRuns,patternsDetails,patternsBritishBTBDetails,patternsBTBDetails,patternsMic_PinDetails,patternsMicrotiDetails,patternsPinnipediiDetails,refName,qth,pathAux,thMinGoodCov,thCovProp,thqualsnp,thqualnonsnp):    
    pmeanCov=strainsDetailsTittle.index('MeanDepth')
    pfileName=strainsDetailsTittle.index('Sample')
    name=[strainDetails[pfileName]]
    meanCov=float(strainDetails[pmeanCov])
 
    strainStatsFileName=strainDetails[pfileName]+".pileup.vcf"
    
    posToExtract=map(int,patternsDetails[0][1:])
    posToExtractBritishBTB=map(int,patternsBritishBTBDetails[0][1:])
    posToExtractBTB=map(int,patternsBTBDetails[0][1:])
    posToExtractMic_Pin=map(int,patternsMic_PinDetails[0][1:])
    posToExtractMicroti=map(int,patternsMicrotiDetails[0][1:])
    posToExtractPinnipedii=map(int,patternsPinnipediiDetails[0][1:])


# change flags in this section
# If the match to the first reference dataset is above a given threshold the process moves to the next
# This is iterated depending on the outcomes
    [strainGSSInfo,strainGSSBritishBTBInfo,strainGSSBTBInfo,strainGSSMic_PinInfo,strainGSSMicrotiInfo,strainGSSPinnipediiInfo]=getSnpsStatsStrain(strainStatsFileName,[posToExtract,posToExtractBritishBTB,posToExtractBTB,posToExtractMic_Pin,posToExtractMicroti,posToExtractPinnipedii],pathAux,thMinGoodCov,thCovProp,thqualsnp,thqualnonsnp)  
    if meanCov >=qth:
        BTB=getBestMatchPattern(patternsBTBDetails,strainGSSBTBInfo)[0]
        if BTB[0]=="bTB" and BTB[2]>=70:
            flag="bTB"
            BritishBTB=getBestMatchPattern(patternsBritishBTBDetails,strainGSSBritishBTBInfo)[0]
            if BritishBTB[0]=="BritishbTB" and BritishBTB[2]>=70:
                flag="BritishbTB"
            else:
                flag="nonBritishbTB"
            [maxPat,strainQ]=getBestMatchPattern(patternsDetails,strainGSSInfo)
            maxPat=strainDetails+[flag]+maxPat
            strainQ=[name]+strainQ
            return [maxPat,strainQ]
        else:
            Mic_Pin=getBestMatchPattern(patternsMic_PinDetails,strainGSSMic_PinInfo)[0]
            if Mic_Pin[0]=="MicPin" and Mic_Pin[2]>=70:
                Microti=getBestMatchPattern(patternsMicrotiDetails,strainGSSMicrotiInfo)[0]
                if Microti[0]=="Microti" and Microti[2]>=70:
                    flag="Microti"
                    maxPat=strainDetails+[flag]+Microti
                    return [maxPat,"NA"]
                else:
                    Pinnipedii=getBestMatchPattern(patternsPinnipediiDetails,strainGSSPinnipediiInfo)[0]
                    if Pinnipedii[0]=="Pinnipedii" and Pinnipedii[2]>=70:
                        flag="Pinnipedii"
                        maxPat=strainDetails+[flag]+Pinnipedii
                        print(maxPat)
                        return [maxPat,"NA"]
                    else:
                        flag="MicPin"
                        maxPat=strainDetails+[flag]+Mic_Pin
                        return [maxPat,"NA"] 
            else:
                flag="nonbTB"
                maxPat=strainDetails+[flag]+BTB
                return [maxPat,"NA"]
    else:
        flag="LowCoverage"
        maxPat=strainDetails+[flag]+6*["NA"]
        print(maxPat)
        return [maxPat,"NA"] 


def getSnpsStatsStrain(strainStatsFileName,listas,pathAux,thMinGoodCov,thCovProp,thqualsnp,thqualnonsnp):
    print("loading "+ strainStatsFileName)
    
    fileIn = open(strainVCF, 'r')
    csv=fileIn.readlines()
    fileIn.close()
    csv=dict([(t.split()[1],t.split()) for t in csv if t[0]!="#" and "INDEL" not in t])
    listasOut=[]
    for lista in listas:
        out=[]
        for i in lista:
            nu=[i,"n","NA","NA","NA","NA","NA","NA"]
            try:
                line=csv[str(i)]
            except:
                out=out+[nu]
                continue
            ref=line[3]
            if line[4]==".":
                alt=line[3]
            else:
                alt=line[4]
            qual=float(line[5])
            det=line[7].split(";")
            if "DP4=" in line[7]:
                call="n"
                gcovRF,gcovRR,gcovAF,gcovAR=map(int,[s for s in det if "DP4=" in s][0].split("DP4=")[1].split(","))
                if ref.upper() in ["A","C","G","T"] and alt.upper() in ["A","C","G","T"]:
                    if gcovRF+gcovAF>0 and gcovRR+gcovAR>0:
                        if ref.upper()==alt.upper() and qual>thqualnonsnp and gcovRF+gcovRR>thMinGoodCov and float(gcovAF)/(gcovRF+gcovAF)<thCovProp and float(gcovAR)/(gcovRR+gcovAR)<thCovProp:
                            call=ref.upper()
                        if ref.upper()!=alt.upper() and qual>thqualsnp and gcovAF+gcovAR>thMinGoodCov and float(gcovRF)/(gcovRF+gcovAF)<thCovProp and float(gcovRR)/(gcovRR+gcovAR)<thCovProp:
                            call=alt.upper()
                nu=[i,call,qual,gcovRF+gcovRR+gcovAF+gcovAR,gcovRF,gcovRR,gcovAF,gcovAR]
            out=out+[nu]
        listasOut=listasOut+[out]
    print("snps Extracted in lists of length:")
    itemLen = [len(item) for item in listasOut]
    print(itemLen)
    os.system("rm "+strainVCF)
    return listasOut


def isACGT(c):
    return c.upper() in ["A","C","G","T"]
    
def readTVSFile(fname):
    fileIn = open(fname, 'rb')
    data = csv.reader(fileIn)
    dataOut = [row for row in data]
    fileIn.close()
    return dataOut

def getBestMatchPattern(patternsDetails,strainGSSInfo):
    print("matching positions:")
    hal=[int(x) for x in patternsDetails[0][1:]]==[x[0] for x in strainGSSInfo]
    if hal==False:sys.exit("matching positions is false")
    strainGSS=[x[1] for x in strainGSSInfo]
    refPat=patternsDetails[2][1:]
    maxVal=0
    maxPat=[]
    maxgroPat=[]
    maxgroRes=[]
    for pattern in patternsDetails[3:]:
        groPat=pattern[1:]
        groRes=comparePatterns(refPat,strainGSS,groPat)
        comp=[pattern[0],len(groPat)]+groRes[0]
        if comp[2]>=maxVal:
            maxVal=comp[2]
            maxPat=comp
            maxgroPat=groPat
            maxgroRes=groRes
        if maxPat==[]:
            maxPat=comp
            maxgroPat=groPat
            maxgroRes=groRes
    strainQ=[[maxPat[0]]]+[maxPat[1:]]+[reduce(operator.add,[[y],[z],[x[1]],x[2:]]) for y,x,z in zip(maxgroPat,strainGSSInfo,maxgroRes[1])]
    return [maxPat,strainQ]


print(TBRun)
refName=refName.split(".")[0]
strainDetailsFile=instats
pathResutls=os.path.join(TBRun,"Stage1")
if not os.path.exists(pathResutls): os.makedirs(pathResutls)
pathAux=os.path.join(pathResutls,"Aux")
if not os.path.exists(pathAux): os.system("mkdir "+ pathAux)


strainsInfo=readTable(strainDetailsFile,',')
pfileName=strainsInfo[0].index('Sample')
pmeanCov=strainsInfo[0].index('MeanDepth')
ppermap=strainsInfo[0].index('%Mapped')
totalReads=strainsInfo[0].index('NumRawReads')
genomeCov=strainsInfo[0].index('GenomeCov')
outcome=strainsInfo[0].index('Outcome')
strainsInfo=listT(strainsInfo)
strainsDetails=listT([strainsInfo[pfileName][1:],strainsInfo[genomeCov][1:],strainsInfo[pmeanCov][1:],strainsInfo[totalReads][1:],strainsInfo[ppermap][1:],strainsInfo[outcome][1:]])

strainsDetails=[['Sample','GenomeCov','MeanDepth','NumRawReads','pcMapped','Outcome',]]+strainsDetails
print("Processing "+ str(len(strainsDetails))+" samples")

patternsDetails=listT(readTable(os.path.join(pathPatterns,patternsDetailsFile),","))
patternsBritishBTBDetails=listT(readTable(os.path.join(pathPatterns,patternsBritishBTBFile),","))
patternsPinnipediiDetails=listT(readTable(os.path.join(pathPatterns,patternsPinnipediiFile),","))
patternsMic_PinDetails=listT(readTable(os.path.join(pathPatterns,patternsMic_PinFile),","))
patternsMicrotiDetails=listT(readTable(os.path.join(pathPatterns,patternsMicrotiFile),","))
patternsBTBDetails=listT(readTable(os.path.join(pathPatterns,patternsBTBFile),","))


maxPats=[strainsDetails[0]+["flag","group","CSSTested","matches","mismatches","noCoverage","anomalous"]]
maxPatsQ=[[[patternsDetails[0][0]]]+[["PredGenotype"],["M-MM-N-A"]]+[[x] for x in patternsDetails[0][1:]],[[patternsDetails[1][0]]]+[[""],[""]]+[[x] for x in patternsDetails[1][1:]],[[patternsDetails[2][0]]]+[[""],[""]]+[[x] for x in patternsDetails[2][1:]]]

outFileName="_stage1.csv"

for strainDetails in strainsDetails[1:]:
    print(strainDetails)
    [maxPat,strainQ]=findGenotypeOneSample(strainsDetails[0],strainDetails,pathTBRuns,patternsDetails,patternsBritishBTBDetails,patternsBTBDetails,patternsMic_PinDetails,patternsMicrotiDetails,patternsPinnipediiDetails,refName,qth,pathAux,thMinGoodCov,thCovProp,thqualsnp,thqualnonsnp)
    maxPats=maxPats+[maxPat]
    if strainQ!="NA":
        maxPatsQ=maxPatsQ+[strainQ]

os.system("rm -R "+pathAux)

writeCSV(outFileName,maxPats)
maxPats=[maxPats[0]]+sorted(maxPats[1:],key=lambda x: x[0])
writeCSV(outFileName,maxPats)