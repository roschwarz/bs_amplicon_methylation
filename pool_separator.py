from pylib.container import seq_file
from pylib.container import amplicon
from pylib.seq_methods import primer_check, rev_comp
import pylib.administration as admin
import sys

def read_amplicon_informations(info_file):
    
    amplicon_list = []
    print("read amplicon information...")
    with open(info_file, "r") as input_file:

        for line in input_file:
            
            ampl = amplicon(line.split(";")[-1].strip(), line.split(";")[2],line.split(";")[3])
            amplicon_list.append(ampl)
            
    return amplicon_list
            
def demultiplex(pool, amplicon_list):
    
    print("start to demultiplex the pool...")
    for record in pool.records:

        for ampl in amplicon_list:
            
            check_return = primer_check(record.seq, ampl.fwd_primer, ampl.rev_primer)
            
            # bring all in one direction 
            if check_return[0] and check_return[1]==0:   
                ampl.add_fastq_record(record)
            elif check_return[0] and check_return[1]==1:
                rev_comp(record)
                ampl.add_fastq_record(record)
                


def read_pool(fastq_file):

    pool = seq_file(fastq_file)
    pool.read()
    return pool
    
def main():

    pool = read_pool(sys.argv[1])
    
    amplicons = read_amplicon_informations(sys.argv[2])
    #admin.create_folder("separated_pool")
    demultiplex(pool, amplicons) 
    file_list = []

    print("write amplicon separated files...")

    for ampl in amplicons:

        ampl.write_amplicon_rec()
        file_list.append("%s.fq" %ampl.name) 
    
    admin.move_files(file_list, "separated_pool_%s" %admin.create_time_stamp())
    
if __name__ == "__main__":
    main()
