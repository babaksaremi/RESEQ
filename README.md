# RESEQ

RESEQ is a python script that combines mapping tools and bootstrap to access reproducibility in next generation sequencing data. It uses paired end FASTQ-files as input.

## Dependencies and Installation

### Dependencies:
- Unix based operating system  
- Python 2.7 (recommended)  
- Bowtie2  
- Samtools  
- R 3.3 (optional)  

## Taxonomy and reference Data.

### Reference Sequences:
For virus detection you need the complete reference sequences which are available at ncbi:  
https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/  
Simply download the sequences, concatenate them and if needed remove new lines from sequences.  
_Example Template:_  

\>genome1  
AATTGGCC  
\>genome2  
GGTTAAAC  

**EDIT YOUR PATH IN reseq.py LINE 19**

### Taxonomy Data:
You can download the taxonomy dumb here:
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/  
note that to speed up analysis you can alter the file to only include entries which you are mapping to

**EDIT YOUR PATH IN reseq.py IN LINE: 22**


### Bowtie2:
Follow installation here:  
http://bowtie-bio.sourceforge.net/bowtie2/index.shtml  

Bowtie needs to create index files for mapping:
```
$: bowtie2-build yourReferenceData.fna indexOut
```

**EDIT THIS PATH IN reseq.py IN LINE 16 WITH THE PREFIX**

###   Samtools:
Follow installation here:
http://www.htslib.org/


## Using RESEQ
```
python reseq.py -i _input1.fastq_ -j _input2.fastq 
```
reseq provides parameters to optimize your run. To view all parameters type:
````
python reseq.py -h
````

````
  -h, --help        show this help message and exit
  -i INPUTR1        fastQ input file R1
  -j INPUTR2        fastQ input file R2
  -n LOOPCOUNT      number of repetitions
  -o OUTPUT         output name
  -s SEQUENCECOUNT  number of sequences extracted
  -r REPLACE        "n": to draw with out replacement Default=draw with
                    replacement
  -d DELETETMP      "n": to keep all temporary files. Default value deletes
                    all tmp files
  -S DELETESAM      "n": to keep all sam / bam files. Default value deletes
                    all files
  -R USER           "y": use R to apply mixture model for False/True positive
                    likelihood predictions
  -t THREADS        number of threads
````

### known issues
- if data sets are too small reseq can sometimes draw the same reads from the fastQ files.
To prevent this reseq is implemented with a seed for each bootstrap repitition. To remove this seed you can comment out line 105
- FASTQ files sometimes have different header conventions to distinguish between paired end files. In this version reseq is hard coded for the following:  

@xxxxxx 1:N:0:5  

to change this you can edit how the headers should be distinguished at line 99. For example if your header strukture looks like this: 

@xxxxxx/1  
change line 86 and 99 from:  
````
R2Dict[lines[0].strip("\n").split(" ")[0]]= x
````
to  
````
R2Dict[lines[0].strip("\n").split("/")[0]]= x
````
