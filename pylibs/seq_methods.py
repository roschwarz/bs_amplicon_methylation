########
# Library for sequence handling
#
#
# author: Robert Schwarz
#
#########

import random
import re
import pylibs.administration as admin
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq


def date():
    return admin.create_time_stamp()

def read_fastq_qual(fastq_file):

    print("read quality strings...")

    quality_strings = []
    
    for name, seq, qual in FastqGeneralIterator(open(fastq_file)):
            
            quality_strings.append(qual)

    return quality_strings

    print("read quality strings done...")

def load_seq_file(seq_file):

    print("load %s" %seq_file) 
    records = [] 


    if seq_file.split(".")[-1] in ["fa", "fasta"]:
        
        for rec in SeqIO.parse(seq_file, "fasta"):
            records.append(rec)

    elif seq_file.split(".")[-1] in ["fq", "fastq"]:

        for rec in SeqIO.parse(seq_file, "fastq"):
            records.append(rec)

    else:
        print("no seq file")

    return records

def read_fasta(fasta_file):

    for rec in SeqIO.parse(fasta_file, "fasta"):
        yield rec

def read_fastq(fastq_file):

    for rec in SeqIO.parse(fastq_file, "fastq"):
        yield rec


def read_seq_file(seq_file):

    if seq_file.split(".")[-1] in ["fa", "fasta"]:
        
        for rec in SeqIO.parse(seq_file, "fasta"):
            yield rec

    elif seq_file.split(".")[-1] in ["fq", "fastq"]:

        for rec in SeqIO.parse(seq_file, "fastq"):
            yield rec

    else:
        print("no seq file")

def p_error(phred_score):

    return 10 **(-(ord(phred_score)-33)/10)

def write_fastq(list_of_records, seed, args, time_stamp):

    fastq_file_name = "simulated_reads_c%d_l%d_s%d_%s" %(args.coverage,
            args.read_length, seed, time_stamp) 

    with open("%s.fq" %fastq_file_name, "w") as output_handle:
        for rec in list_of_records:
            output_handle.write("@%s\n" %rec.name)
            output_handle.write("%s\n" %rec.seq)
            output_handle.write("+\n")
            output_handle.write("%s\n" %rec.qual)

def mutate_read(read, random_qual_string):

    mutated_read = ''
    mutated_bases = 0

    for idx, base in enumerate(read):

        if base in "ACGT" and random.random() < p_error(random_qual_string[idx]):

            mutated_read = mutated_read + random.choice("ATGC".replace(base, ""))
            mutated_bases += 1
        
        else:
            mutated_read = mutated_read + base
    
    return mutated_read, mutated_bases

def random_seq_position(seq_length):

    return random.choice(list(range(0,seq_length)))

def calc_coverage_dependent_readcount(seq_length, read_length, coverage):
    
    return round(coverage*seq_length/read_length)

def draw_sub_sequences(seq, sub_seq_count, sub_seq_length):
   
    sub_seq_list = [] 
    if len(seq)-sub_seq_length >= 0:
        
        for read in range(1,sub_seq_count):

            start_position = random_seq_position(len(seq)-sub_seq_length)

            sub_seq = seq[start_position:start_position+sub_seq_length]
           #Is yield also possible? 
            sub_seq_list.append(sub_seq)
    else:
       my_error = ValueError("read length longer than the reference sequence")
       raise my_error

    return sub_seq_list

def change_strand_chance(rec, possibility):
    # To-Do: If fwd is not stored in the name
    if random.random() > possibility:
        rec.name = rec.name.replace("fwd", "rev")
        rec.seq =Seq(str(rec.seq)).reverse_complement()


def find_motif(seq, motif):
    
    hits = []
    
    for hit in re.finditer(r"%s" %str(motif), str(seq)):
        hits.append(hit)

    return hits 

def primer_check(seq, primer_1, primer_2):

    fwd = find_motif(seq, primer_1)
    fwd_comp = find_motif(seq, Seq(str(primer_1)).reverse_complement())
    rev = find_motif(seq, primer_2)
    rev_comp = find_motif(seq, Seq(str(primer_2)).reverse_complement())
    
    if fwd and rev_comp:
        return True, 0
    elif rev and fwd_comp:
        return True, 1
    else:
        return False, 0

def rev_comp(record):

    record.seq = record.seq.reverse_complement()
