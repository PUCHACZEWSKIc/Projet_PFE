class Fasta_extract : 
    
    def __init__(self, file) : 
        self.__file = file

    def sequences(self) : 
        sequence_list = []
        act_line = ''
        with open(self.__file) as read : 
            for seq in read : 
                seq = seq.strip("\n")
                if seq == "" : 
                    pass
                elif seq[0] == '>' :
                    sequence_list.append(act_line)
                    act_line = ''
                else :
                    act_line += seq
        sequence_list.append(act_line)
        return sequence_list[1:len(sequence_list)]
