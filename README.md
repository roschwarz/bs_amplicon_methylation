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

Script that determines the methylation rate of each defined CpG contained in the
respective amplicon for the X and Y chromosome respectively. The script requires 
a .txt file with information for the forward and reverse primer, the X- and 
Y-specific amplicon length, the tag length, the position of all CpGs within 
the amplicon (comma-separated), a location of a non-CpG cytosine as well as 9
to 10 bp X (p1) and Y (p2) signature sequences (9 – 10 bp) in a in the 
following pattern:

forward_primer:<br/>
reverse_primer:<br/>
expected_length_p1:<br/>
expected_length_p2:<br/>
tag_length:<br/>
cpg_pos_p1:<br/>
cpg_pos_p2:<br/>
cx_pos_p1:<br/>
cx_pos_p2:<br/>
p1:<br/>
p2:<br/>

The methylation value of the non-CpG cytosine is provided in the 
column cx_methylation which can be used to correct calculated methylation 
values.

```
python3 methi.py -d fastq-files/ -i sample.txt -s
```

-d - directory with all fastq-files
-i - the information file with the patter like described above
-s - separates the reads into sex-specific sets
