#Ce fichier a pour but d'aligner des séquences sur une séquence consensus (fournit par le logiciel Jalview)

#Séquence consensus : 

seq_cons = ""

with open("Sequence_consensus.txt") as file :
    for line in file : 
        seq_cons += line

class align : 
    
    def __init__(self, seq_cons, dict_kmer) : 
        self.__seq_cons = seq_cons
        self.__dict_kmer = dict_kmer
    

 
