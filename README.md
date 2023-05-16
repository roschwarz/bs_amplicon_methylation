# bs_amplicon_methylation

Scripts to call methylation rates of CpGs in amplicons. 

**pool_separator.py** 

Separates amplicon reads stored in mixed .fastq-files based on 3 bp-tags added 
to the amplicon reads. The tool takes as input a .fastq-file and an info file 
(.csv-format) containing row-separated information about the samples by a comma-
separated list of pool, tag, forward-primer, reverse-primer, and sample name,
which identifies each sample.

```
python3 pool_separator.py <.fastq-file> <info-file.csv>
```

**methi**

Script that determines the methylation rate of each named CpG contained in the
respective amplicon, separating X and Y chromosomal reads. The script requires 
the forward and reverse primer, the X- and Y-specific amplicon length, the 
tag length, the position of all CpGs within the amplicon (comma-separated), 
a location of a non-CpG cytosine as well as 9 to 10 bp X (p1) and Y (p2) 
signature sequences (9 – 10 bp) in a .txt-file in the following pattern:

forward_primer:
reverse_primer:
expected_length_p1:
expected_length_p2:
tag_length:3
cpg_pos_p1:
cpg_pos_p2:
cx_pos_p1:
cx_pos_p2:
p1:
p2:

The methylation value of the non-CpG cytosine is used to correct the calculated
methylation values.

```
python3 methi.py -d fastq-files/ -i sample.txt -s
```
-s - separates the reads into sex-specific sets
