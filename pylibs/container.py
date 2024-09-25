import pylibs.administration as admin
from Bio import SeqIO
from pylibs.seq_methods import load_seq_file, primer_check, find_motif
#from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

class fasta_records():

    def __init__(self):

        self.fasta_records = []
        self.number = 0
    
    def add_fasta_rec(self, rec):

        self.fasta_records.append(rec)
        self.number += 1

    def write_records():
        
        pass

class fasta_record():

    def __init__(self, seq):

        self.id = seq_id
        self.seq = seq
    
    def __repr__(self):

        pass

class amplicon():

    def __init__(self, name, fwd_primer, rev_primer):

        self.name = name
        self.fwd_primer = fwd_primer
        self.rev_primer = rev_primer
        self.fastq_records = []

    def add_fastq_record(self, record):

        self.fastq_records.append(record)
    
    def write_amplicon_rec(self):
        
        print("write records for amplicon %s..." %self.name)
        with open("%s.fq" %self.name, "w") as output_handle:
            SeqIO.write(self.fastq_records, output_handle, "fastq")

class log_file():

    def __init__(self, tool, running_command):
        self.tool = tool
        self.cmd = running_command
        self.beginning = admin.create_time_stamp()

    def add_storage(self, aim_path):

        aim_path = aim_path + "%s_%s" %(self.tool, self.beginning)
        self.data_storage = aim_path 

    def add_error():
        pass

    def write_log(self):
        
        with open("%s_%s.log" %(self.tool, self.beginning), "w+") as output_file:

            output_file.write(' '.join(self.cmd))

            output_file.write('\n%s' %self.data_storage)

class seq_file():

    def __init__(self, file_name):

        self.file = file_name.split("/")[-1]
        self.path = file_name
        self.entries = 0
        self.records = []
        
    def read(self):
        
        self.records = load_seq_file(self.path)
        self.entries = len(self.records)

    def add_entrie(self):
        
        self.entries += 1

    def filter_primer(self, primer_1, primer_2):
        
        temp = []

        for rec in self.records:

             if primer_check(rec.seq, primer_1, primer_2):
                 temp.append(rec)
                 
        self.records = temp
        self.entries = len(self.records)

    def split_by_motif(self, motif_1, motif_2):
        
        self.records_motif_1 = []
        self.records_motif_2 = []
        self.records_leftover = []

        for rec in self.records:
            if find_motif(rec.seq, motif_1):
                self.records_motif_1.append(rec)
            elif find_motif(rec.seq, motif_2):
                self.records_motif_2.append(rec)
            else: 
                self.records_leftover.append(rec)

    def write_subsets(self):
      
        aim_folder = self.path[:-(len(self.file))] +"subsets"
        admin.create_folder(aim_folder)
        file_1 = aim_folder +"/" + (self.file).split(".")[0] + "_p1.fq"
        file_2 = aim_folder + "/" + (self.file).split(".")[0] + "_p2.fq"
        
        SeqIO.write(self.records_motif_1, "%s" %file_1, "fastq")
        SeqIO.write(self.records_motif_2, "%s" %file_2, "fastq")
    
    def seq_log(self):
        pass

    def tag_cut(self, tag_length):

        #for rec in self.records:
        #    rec = rec[3:-3]
            #rec.qual = str(rec.seq)[tag_length:-tag_length]
        pass

