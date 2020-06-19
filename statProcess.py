def statSummary():
	import subprocess
	import numpy
	import statistics
	from math import sqrt



	ref=open("normalRun","r") #bowtie standart run
	refDict={}
	# eventuell header beim einlesen skippen
	for line in ref.readlines():
		refDict[line.split("\t")[0]]=line.split("\t")[1] # redDict contains IDs of original bowtie run as keys and dna read count as value 
	ref.close()

	out=open("bootsTrapSummary","w")  # output file
	#out.write("ID\tMin-MaxCount;medianCount\tMinMean-MaxMean;medianMean\tminMedian-MaxMedian;medianMedian\tMinPer-MaxPer;medianPer\n")
	out.write("ID\tOriginalReadCount\tMinBTReadCount\tMaxBTReadCount\tMedianBTReadCount\tMeanBTReadCount\tMinBTPer\tMaxBTPer\tMedianBTPer\n")
	for i in refDict.keys():
		cmd="grep "+i+" *_stats"+" | cut -f5-"
		grep=subprocess.Popen([cmd],shell=True,stdout=subprocess.PIPE)
		stdout,error=grep.communicate()
		stdout=stdout.split("\n")
		stdout.pop(-1)
		stdout=[g.split("\t") for g in stdout]

		if len(stdout)==0:
			continue
		count=[]
		per=[]

		for t in range(len(stdout)):
			count.append(stdout[t][0])
			per.append(stdout[t][3])

		count=map(float,count)
		per=map(float,per)
		
		out.write(i+"\t"+refDict[i]+"\t"+str(min(count))+"\t"+str(max(count))+"\t"+str(numpy.median(count))+"\t"+str(numpy.mean(count))+"\t"+str(min(per))+"\t"+str(max(per))+"\t"+str(numpy.median(per))+"\n")