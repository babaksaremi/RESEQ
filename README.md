# RESEQ

RESEQ is a python script that combines mapping tools and bootstrap to access reproducibility in next generation sequencing data. It uses paired end FASTQ-files as input.

## Dependencies and Installation

### Dependencies:
Unix based operating system
Python 2.7 (recommended)
Bowtie2
Samtools
R (optional)

## Taxonomy and reference Data.

### Reference Sequences:
For virus detection you need the complete reference sequences which are available at ncbi:
https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/
Simply download the sequences, concatenate them and if needed remove new lines from sequences.
Example Template:
>genome1
AATTGGCC
>genome2
GGTTAAAC
…
**EDIT YOUR PATH IN reseq.py LINE 19

### Taxonomy Data:
You can download the taxonomy dumb here:
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/
note that to speed up analysis you can alter the file to only include entries which you are mapping to

**EDIT YOUR PATH IN reseq.py IN LINE: 22


### Bowtie2:
Follow installation here:
bowtie-bio.sourceforge.net/bowtie2/index.shtml

Bowtie needs to create index files for mapping:
```
>user$: bowtie2-build yourReferenceData.fna indexOut
```

**EDIT THIS PATH IN reseq.py IN LINE 16 WITH THE PREFIX

###   Samtools:
Follow installation here:
http://www.htslib.org/

