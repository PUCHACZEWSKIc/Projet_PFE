#Ce fichier a pour but d'aligner des séquences sur une séquence consensus (fournit par le logiciel Jalview)

#Séquence consensus : 

import parasail

seq_cons = ""

with open("Sequence_consensus.txt") as file :
    for line in file : 
        seq_cons += line

class align : 
    
    def __init__(self, seq_cons, seq_list) : 
        self.__seq_cons = seq_cons
        self.__seq_list = seq_list
    
    def align(self) :
        cigar_list = []
        for seq in self.__seq_list :
            test = parasail.nw_trace(seq, self.__seq_cons, 10, 1, parasail.dnafull)
            cigar_list.append(test.cigar.decode)
            traceback = test.get_traceback()
            print(traceback)
        return cigar_list

from Extract_fasta import Fasta_extract 
import cigar

testaille = Fasta_extract('sequences_SARS-CoV-2_ech.fasta')

data_test = testaille.sequences()

testouille = align(seq_cons, data_test)

test = testouille.align()
for i in test : 
    cig = cigar.Cigar(i)
    print(cig)
