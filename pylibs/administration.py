import time
import re
import random
import subprocess
import shlex

def create_time_stamp():

    stamp = time.strftime("%Y%m%d_%H%M%S", time.gmtime())

    return stamp

def set_seed(seed):

    if seed == -1:
        return random.choice(list(range(1,1000)))
    else:
        return seed

class simulator_data_collector():
    
    def __init__(self):
        
        self.time_stamp = create_time_stamp()
        self.processed_references = 0
        self.simulated_reads = 0
        self.mutated_bases = 0

    def add_values(self, read_record):

        self.simulated_reads = self.simulated_reads +\
        read_record.simulated_reads

        self.mutated_bases = self.mutated_bases + read_record.mutated_bases 

    def add_ref(self):

        self.processed_references += 1

    def write_logfile(self):

        file_name = "run_%s.log" %self.time_stamp
        
        with open(file_name, "w") as output_handle:
            output_handle.write("Information about run %s\n" %self.time_stamp)
            output_handle.write("processed references: %d\n" \
                    %self.processed_references)
            output_handle.write("amount of simulated reads: %d\n" \
                    %self.simulated_reads)
            output_handle.write("mutated bases: %d\n" %self.mutated_bases)


def remove_multi_blanks(line):

    return re.sub('\s+', ' ', line).strip()

#To-Do check if the aim_path already exists
def move_files(file_list, aim_path):
    
    mv_command = "mv"
     
    create_folder(aim_path) 
    for f in file_list: 
        mv_command = mv_command +" %s" %f
    
    mv_command = mv_command + " %s" %aim_path
    
    print(mv_command)
    subprocess.run(shlex.split(mv_command))

def create_folder(aim_dict):
    
    mkdir_command = "mkdir -p "
    mkdir_command = mkdir_command + aim_dict
    
    print("create the storage directory\n%s" %mkdir_command)
    
    subprocess.run(shlex.split(mkdir_command))

def read_infofile(info_file):

    fobj=open(info_file)
    info={}

    for line in fobj:
        
        if re.search(r",", line):

            line=line.strip()
            array=line.split(":")[1].split(",")
            info[line.split(":")[0]]=array

        else:

            line=line.strip()
            info[line.split(":")[0]]=line.split(":")[1]
    
    return info
