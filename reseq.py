from optparse import OptionParser
from multiprocessing import Process, Manager
import subprocess
import pysam
import numpy
import os
import psutil
from getSamData import getSam
from statProcess import statSummary
import os


#insertPath (pls only use abosulte paths)

#bowtie path to index files +name
bowtieIndex="/PATH/TO/bowtieIndex/virusIndex" 

# reference multi fasta file
source="/PATH/TO/fullVirusReference.fna"

# ncbi taxonomy file
taxSource="/PATH/TO/virusTaxonomy.dmp"

#option parser 
parser = OptionParser()
parser.add_option("-i", dest="inputR1",default="fastQR1.fastq",help="fastQ input file R1")
parser.add_option("-j", dest="inputR2",default="fastQR2.fastq",help="fastQ input file R2")
parser.add_option("-n",dest="loopCount",default="5", help="number of repetitions")
parser.add_option("-o",dest="output",default="out", help="output name")
parser.add_option("-s",dest="sequenceCount",default="5000",help="number of sequences extracted")
parser.add_option("-r",dest="replace",default="y", help="\"n\": to draw with out replacement Default=draw with replacement")
parser.add_option("-d",dest="deletetmp",default="y",help="\"n\": to keep all temporary files. Default value deletes all tmp files")
parser.add_option("-S",dest="deleteSam",default="y",help="\"n\": to keep all sam / bam files. Default value deletes all files")
parser.add_option("-R",dest="useR",default="n",help="\"y\": use R to apply mixture model for False/True positive likelihood predictions" )
parser.add_option("-t",dest="threads",default="30",help="number of threads")
(options, args) = parser.parse_args()

print("processing input")
fastQDict={}
R2Dict={}
reference={}

#normal Run. Standart mapping of fastQ files for later references
normal="y" # set to n to skip original mapping
print("mapping Original fastQ")
if normal=="y":
	bowtieCMD="bowtie2 -p "+options.threads+" -x "+bowtieIndex+" -1 "+options.inputR1+" -2 "+options.inputR2+" -S normalRun.sam"
	bowtie=subprocess.Popen([bowtieCMD],stdout=subprocess.PIPE,shell=True)
	stdout,error=bowtie.communicate()
	getSam()



if options.deletetmp == "y":
	os.remove("normalRun.sam")
	os.remove("normalRun.bam")
	os.remove("normalRunsorted.bam")
	os.remove("normalRunsorted.bam.bai")


#reference genome.Processing of fullVirusReference.fna file 
print("building reference dictionary")
with open(source,"r") as infile:
	for line in infile:
		if ">" in line:
			sequence=infile.next()
			reference[line.split(" ")[0].strip(">")]=[line.split(" ")[0].strip(">"),line[len(line.split(" ")[0])+1:].strip("\n"),len(sequence),sequence.strip("\n")]
infile.close()

# building ncbi taxonomy dictionary from provided file
taxonomy={}
with open(taxSource,"r") as infile:
	for line in infile:
		if len(line.split(";"))>2:
			tax=line.split("|")[2].strip("\t")
			taxonomy[line.split("|")[1].strip("\t")]=[tax.split(";")[0],tax.split(";")[1]]
infile.close()

# Save both fast Qfiles as dictionary. Later the bootstrap funktion can draw reads from these dictionarys 
c=0 
with open (options.inputR1,"r") as infile:
	lines=[]
	for line in infile:
		lines.append(line.replace("\r","")) # remove newlines on windows machine
		if len(lines)==4: # one read is represented by 4 lines in a fastQ file
			lines[0]=lines[0].split(" ")[0]+"\n" # split at header info. CHANGE HERE if depending on header delimiter of fastQ file.
			x="".join(lines[0:4]) # join 4 lines together (1 read)
			fastQDict[str(c)]= x # save read in dictionary counter as dict. Key
			lines=[]
			c+=1
infile.close()
#R2 fastQ
with open (options.inputR2,"r") as infile:
	lines=[]
	for line in infile:
		lines.append(line.replace("\r","")) 
		if len(lines)==4:
			x="".join(lines[0:4]) # here use header ID as dict. Key
			R2Dict[lines[0].strip("\n").split(" ")[0]]= x #key for dict has to match dictionary above. look into header id and split at information of forward/backward
			lines=[]
infile.close()

#bootstrap function
def BT(dictA,dictB,loop,count,out,core,delTmp,replace):  # bootstrap
	numpy.random.seed(loop) # set seed --> for small dataSets --> README
	if replace != "y":  # generate a list of random numbers between one and total read count
		sample=numpy.random.choice(len(fastQDict),int(count),replace=False)
	else: 
		sample=numpy.random.choice(len(fastQDict),int(count))

	outputR1= open(out+str(loop)+"R1.fastq","w") 
	outputR2= open(out+str(loop)+"R2.fastq","w")
	for i in sample:  # with the list of random numbers (matches dict. key) write the reads into output file. Use R1 dictionary to find paired end read from R2 dictionary
		if dictA[str(i)].split("\n")[0].split("/")[0] in dictB: # matches the keys due to header norm this part is special for example data (" ") split here
			outputR1.write(dictA[str(i)])
			outputR2.write(dictB[dictA[str(i)].split("\n")[0].split("/")[0]])
		else:
			outputR1.write(dictA[str(i)])

	outputR1.close()
	outputR2.close()
	

# process mapping out to indexed and sorted bam file by using samtools. Calculate stats for each mapped fastQ pair. Complement information with taxonomy file
def stats(out,loop,core,taxonomyDict,refDict,delSam):

	bamCmd="samtools view -bS -@ "+core+" "+out+str(loop)+".sam > "+str(loop)+"tmp.bam"
	sortCmd="samtools sort -@ "+core+" "+str(loop)+"tmp.bam > "+str(loop)+"sorted.bam"
	indexCmd="samtools index -@ "+core+" "+str(loop)+"sorted.bam"
	
	bam=subprocess.Popen([bamCmd],shell=True,stdout=subprocess.PIPE)
	stdout,error=bam.communicate()
	sort=subprocess.Popen([sortCmd],shell=True,stdout=subprocess.PIPE)
	stdout,error=sort.communicate()
	index=subprocess.Popen([indexCmd],shell=True,stdout=subprocess.PIPE)
	stdout,error=index.communicate()

	statOut=open(out+str(loop)+"_stats","w")
	statOut.write("Species_ID\tSpecies\tPhylum\tSubphylum\tDNA.read_count\tgenome_length\tmean_coverage_per_base\tmedian_coverage_per_base\tcoverage.percent\n")
	
	depth=pysam.depth("-a",str(loop)+"sorted.bam").split("\n")
	depth.pop(-1)
	zeroCounter=0	#count positions with zero reads
	ID="null"
	seq=0
	cov=[]

	for line in depth:
		if line.split("\t")[0] != ID:
			if ID != "null":
			
				if refDict[ID][1].split(",")[0] in taxonomyDict:
					tax=taxonomyDict[refDict[ID][1].split(",")[0]]
				else:
					tax=["NA","NA"]
			
			
				statOut.write(ID+"\t"+reference[ID][1].split(",")[0]+"\t"+tax[0]+"\t"+tax[1]+"\t"+pysam.view("-c",str(loop)+"sorted.bam",ID).strip("\n")+"\t"+str(seq)+"\t"+str(numpy.mean(cov))+"\t"+str(numpy.median(cov))+"\t"+str((float(seq)-float(zeroCounter))/float(seq)*float(100))+"\n")
			ID=line.split("\t")[0]
			seq=0
			zeroCounter=0
			
		elif line.split("\t")[2] == "0":
			zeroCounter+=1
			
		cov.append(int(line.split("\t")[2]))
		seq+=1
	statOut.close()


	if delSam == "y":
		os.remove(str(loop)+"sorted.bam")
		os.remove(str(loop)+"tmp.bam")
		os.remove(str(loop)+"sorted.bam.bai")
		os.remove(out+str(loop)+".sam")


	
print("bootstrap")
# generate for each bootstrap step and indipendent proccess. --> parallelized 
process=[]
mem=(psutil.virtual_memory().available >> 30) / 20  #{ get ram and use bit shift operator in gb
for i in range(int(options.loopCount)):
	    
	while mem == len(process):								# this part stops ram from overflowing. one process here is assumed to need 20 gb of data and wiht our ram available we can only start "so many"
		process=[j for j in process if j.is_alive()==True]  #} this part can be comment out to remove the limiter (same for below)
	p=Process(target=BT,args=(fastQDict,R2Dict,i,options.sequenceCount,options.output,options.threads,options.deletetmp,options.replace)) 
	p.start()
	process.append(p)
for p in process: # continue only if all bootstrap processes are finished 
	p.join()

print("mapping")
# mapp all generated fastQ files with bowtie
for i in range(int(options.loopCount)):
	bowtieCMD="bowtie2 -p "+options.threads+" -x "+bowtieIndex+" -1 "+options.output+str(i)+"R1.fastq"+" -2 "+options.output+str(i)+"R2.fastq"+" -S "+options.output+str(i)+".sam"
	bowtie=subprocess.Popen([bowtieCMD],stdout=subprocess.PIPE,shell=True)
	stdout,error=bowtie.communicate()
	if options.deletetmp == "y":
		os.remove(options.output+str(i)+"R1.fastq")
		os.remove(options.output+str(i)+"R2.fastq")	


	
print("calculating stats")

mem=(psutil.virtual_memory().available >> 30) / 28     #{ get ram and use bit shift operator in gb

process=[]
for i in range(int(options.loopCount)):
	while mem == len(process):								# this part stops ram from overflowing. one process here is assumed to need 20 gb of data and wiht our ram available we can only start "so many"
		process=[j for j in process if j.is_alive()==True]  #} this part can be comment out to remove the limiter (same for below)

	p=Process(target=stats,args=(options.output,i,options.threads,taxonomy,reference,options.deleteSam))
	p.start()
	process.append(p)
for p in process:
	p.join()
	
#stat summary over all Runs	
statSummary()

#mixtools
if options.useR =="y":
	R= subprocess.call("Rscript mixtools.R --input bootsTrapSummary --output bootsTrapSummaryMixtureModel",shell=True)

print("done")
