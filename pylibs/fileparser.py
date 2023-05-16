class align_entry():

    def __init__(self):
        
        self.query_seq = "NA"
        self.start = "NA"
        self.end = "NA"
        self.family = "NA"
        self.order = "NA"
        self.superorder = "NA"
        self.entry_id = "NA"
        self.strand = "NA"
        self.kimura_distance = "NA"
        self.seq_id = "NA"
        self.swscore = "NA"
    
    def set_seq_id(self):

        self.seq_id = ("{}|{}|{}|{}|{}|{}|{}|{}|{}".format(self.query_seq,
                                                self.strand,
                                                self.start,
                                                self.end,
                                                self.order,
                                                self.superorder,
                                                self.family,
                                                self.entry_id,
                                                self.kimura_distance))

    def replace_kimura_seq_id(self, kimura):

        self.seq_id = self.seq_id[:-2] + kimura
    
    def get_strand(self):

        if self.strand == "C":
            return "-"
        else:
            return "+"

    def add_kimura(self, line):
        
        line = line.strip()
        kimura_distance = line.split("=")[1][1:]
        self.kimura_distance = kimura_distance
        self.replace_kimura_seq_id(kimura_distance)

def parse_align_file(file):
    
    align_entries_list = []
    i = -1

    with open(file, "r") as input_handle:

        for line in input_handle:
        
            if re.search(r'^[0-9]', line):

                i += 0

                te_instance = align_entry()
                
                line = line.strip()
                line_split = line.split(" ")
                te_instance.swscore = line_split[0]
                te_instance.query_seq = line_split[4]
                te_instance.start = line_split[5]
                te_instance.end = line_split[6]
                
                te_instance.entry_id = line_split[-1]
                
                # the positive strand is not named. If a C is stored, the length of the list is changing. 
                if line_split[8] == "C":
                    
                    te_instance.strand = "-"
                    te_instance.family = line_split[9].split("#")[0]

                    if re.search(r"/", line_split[9].split("#")[1]):
                        
                        te_instance.order = line_split[9].split("#")[1].split("/")[0]
                        te_instance.superorder = line_split[9].split("#")[1].split("/")[1]
                    
                    else:
                        te_instance.order = line_split[9].split("#")[1]
                        te_instance.superorder = line_split[9].split("#")[1]


                else:
                    
                    te_instance.strand = "+"
                    te_instance.family = line_split[8].split("#")[0]
                    if re.search(r"/", line_split[8].split("#")[1]):

                        te_instance.order = line_split[8].split("#")[1].split("/")[0]
                        te_instance.superorder = line_split[8].split("#")[1].split("/")[1]
                    
                    else:
                        te_instance.order = line_split[8].split("#")[1]
                        te_instance.superorder = line_split[8].split("#")[1]

                te_instance.set_seq_id()
                
                align_entries_list.append(te_instance)
            
            if re.search(r'Kimura', line):
                
                te_instance = align_entries_list[i]
                te_instance.add_kimura(line)

                        
    
    return align_entries_list


