import argparse
import pylibs.administration as admin
import os
import pylibs.container as con
import re
import csv

class compressed_file():

    def __init__(self, sample, cpg_list, cx_position):

        self.sample = sample.split("/")[-1].split(".")[0]
        self.cpg_positions = cpg_list
        self.cx_position = cx_position
        self.cpg_methylation = []
        self.read_count = 0
        self.file = sample 

    def add_cpg_meth(self, cpg, methylation):

        cpg_meth = (cpg, methylation)
        self.cpg_methylation.append(cpg_meth)

    def add_cx_meth(self, cx_position, cx_methylation):

        self.cx_methylation = (cx_position, cx_methylation)

def parse_options():

    parser = argparse.ArgumentParser(prog="cpg_count", description="soon")

    parser.add_argument("-d", type=dir_path, dest="data_path")

    parser.add_argument("-i", type=str, dest="info")

    parser.add_argument("-s",
                        action="store_true",
                        dest="sex",
                        help="divorces the reads into sex specific files \
                                before the cpg will be count")

    args = parser.parse_args()

    return parser, args

def dir_path(path):

    if os.path.isdir(path):

        return path
    else:

        raise NotADirectoryError(path)


def run_cpg_count(info, sample_file, result_path, args):
        
    script_path = os.path.dirname(os.path.realpath(__file__))
    #cpg_count = "/home/lakatos/nielsj/projects/steve/bucksad/cpg_count "
    cpg_count = script_path + "/cpg_count "
    if args.sex:
        
        pattern = sample_file.split(".")[0].split("_")[-1]
        print("pattern: %s" %pattern)
        cmd =cpg_count + \
             " -q %s" %(sample_file)  + \
             " -t %d" %int(info["tag_length"]) + \
             " -F %s" %info["forward_primer"] + \
             " -R %s" %info["reverse_primer"]+ \
             " -l %d" %int(info["expected_length_%s" %pattern]) + \
             " -o %s/%s.csv" %(result_path, sample_file.split("/")[-1].split(".")[0])
        print(cmd)
        os.system(cmd)

    else:
        cmd =cpg_count + \
             " -q %s" %(sample_file)  + \
             " -t %d" %int(info["tag_length"]) + \
             " -F %s" %info["forward_primer"] + \
             " -R %s" %info["reverse_primer"]+ \
             " -l %d" %int(info["expected_length"]) + \
             " -o %s/%s.csv" %(result_path, sample_file.split("/")[-1].split(".")[0])
        print(cmd)
        os.system(cmd)

def split_data_set(path, amplicon_information):

    for file in os.listdir(path):

        sample = con.seq_file(path +file)
        sample.read()
        sample.filter_primer(amplicon_information["forward_primer"],
            amplicon_information["reverse_primer"])
        sample.split_by_motif(amplicon_information["p1"],
            amplicon_information["p2"])
        sample.write_subsets()

def parse_cpg_update(output, cpgs, cx_position):

    compressed_cpg_out = compressed_file(output, cpgs, cx_position)
    fobj=open(output)

    for line in fobj:

        if(line.split("\t")[0].strip() in cpgs):

            cpg = line.split("\t")[0].strip() 
            cpg_meth = line.split("\t")[7].strip() 
            
            compressed_cpg_out.add_cpg_meth(cpg, cpg_meth)

        elif(line.split("\t")[0].strip() == cx_position):
            
            cpx = cx_position
            cpx_meth = int(line.split("\t")[2])/int(line.split("\t")[5])
            compressed_cpg_out.read_count = int(line.split("\t")[5])
            compressed_cpg_out.add_cx_meth(cpx, cpx_meth)

    return compressed_cpg_out


def compress_cpg_output(results, amplicon_information):
    
    parsed_information = []

    for cpg_out in os.listdir(results):
        
        pattern = cpg_out.split(".")[0].split("_")[-1]
        
        #parsed_cpg_out = parse_cpg("%s/%s" %(results,cpg_out),\
        #        amplicon_information["cpg_pos_%s" %pattern],\
        #        amplicon_information["cx_pos_%s" %pattern])
        
        parsed_cpg_out = parse_cpg_update("%s/%s" %(results,cpg_out),\
                amplicon_information["cpg_pos_%s" %pattern],\
                amplicon_information["cx_pos_%s" %pattern])

        #parsed_cpg_out.append(cpg_out.split(".")[0])
        parsed_information.append(parsed_cpg_out)
        
    return parsed_information


def write_compressed_cpg_update(compressed_info, file_name):
    
    with open("final_%s.csv" %file_name, "a") as f:

        writer = csv.writer(f) #, delimiter=";")
 
        for entry in compressed_info:
            
            for cpg in entry.cpg_methylation:

                writer.writerow(['%s' %entry.sample,
                    '%s' %cpg[0],
                    '%s' %cpg[1],
                    "%s" %entry.cx_methylation[1],
                    "%s" %entry.read_count])


def file_init(file_name):

    with open("final_%s.csv" %file_name, "a") as f:

        writer = csv.writer(f, delimiter=";")
        writer.writerow(["sample", "cpg","methylation", "cx_methylation",
            "read_count"])


def main():
    
    parser, args = parse_options()

    amplicon_information = admin.read_infofile(args.info)

    result_dic = "result_%s" %admin.create_time_stamp()
    admin.create_folder(result_dic) 
    
    sample_set = [ ]

    if args.sex: 
        
        print("divorce alleles")

        split_data_set(args.data_path, amplicon_information)
        
        for file in os.listdir(args.data_path +"subsets"):
            
            if re.search("\.fq", file):
            
                sample = con.seq_file(args.data_path +"subsets/" +file)
                run_cpg_count(amplicon_information, sample.path, result_dic,
                        args)
                sample_set.append(sample)
                
    else:
        
        print("run without divorcing")
        
        for file in os.listdir(args.data_path):
            
            if re.search("\.fastq", file):
            
                sample = con.seq_file(args.data_path +file)
                run_cpg_count(amplicon_information, sample.path, result_dic,
                        args)
                sample_set.append(sample)
    
    # write result 
    file_init(result_dic)
    write_compressed_cpg_update(compress_cpg_output(result_dic, amplicon_information),
            result_dic)

if __name__ == "__main__":
    main()
