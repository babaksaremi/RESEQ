def getSam():
	import subprocess
	import numpy
	import pysam
	import subprocess

	
	bamCmd="samtools view -bS -@ 20 normalRun.sam > normalRun.bam"
	sortCmd="samtools sort -@ 20 normalRun.bam > normalRunsorted.bam"
	indexCmd="samtools index -@ 20 normalRunsorted.bam"
			
	bam=subprocess.Popen([bamCmd],shell=True,stdout=subprocess.PIPE)
	stdout,error=bam.communicate()
	sort=subprocess.Popen([sortCmd],shell=True,stdout=subprocess.PIPE)
	stdout,error=sort.communicate()
	index=subprocess.Popen([indexCmd],shell=True,stdout=subprocess.PIPE)
	stdout,error=index.communicate()



	depth=pysam.depth("-a","normalRunsorted.bam").split("\n")
	depth.pop()
	zeroCounter=0	#count positions with zero reads
	ID="null"
	seq=0

	out=open("normalRun","w")
	out.write("reference\tDNA.read_count\tsequence_length\tmeanCov\tmedianCov\tcovPer\n")
	cov=[]

	for line in depth:
		if line.split("\t")[0] != ID:
			if ID != "null":
				out.write(ID+"\t"+pysam.view("-c","normalRunsorted.bam",ID).strip("\n")+"\t"+str(seq)+"\t"+str(numpy.mean(cov))+"\t"+str(numpy.median(cov))+"\t"+str((float(seq)-float(zeroCounter))/float(seq)*float(100))+"\n")
			ID=line.split("\t")[0]
			seq=0
			zeroCounter=0
			
		elif line.split("\t")[2] == "0":
			zeroCounter+=1
		cov.append(int(line.split("\t")[2]))
		seq+=1
	out.close()