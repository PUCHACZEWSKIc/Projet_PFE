#Ce fichier a pour but d'aligner des séquences sur une séquence consensus (fournit par le logiciel Jalview)

#On va utiliser le module parasil, qui nous permet d'aligner nos séquences sur la séquence, et nous fournis la séquence alignée (traceback)
import parasail

#Séquence consensus : 

seq_cons = ""

with open("Sequence_consensus.txt") as file : #Temporaire, ici pour des raisons de test
    for line in file : 
        seq_cons += line

class align : 
    #Notre classe prends en paramètre la séquence consensus et la liste des séquences que l'on veux aligner dessus.
    def __init__(self, seq_cons, seq_list) : 
        self.__seq_cons = seq_cons
        self.__seq_list = seq_list
    
    def align(self) :
        align_list = []
        for seq in self.__seq_list :
            test = parasail.nw_trace_scan(seq, self.__seq_cons, 10, 1, parasail.dnafull) #On utilise parasail pour aligner nos séquences. 
            traceback = test.get_traceback().query #On récupère la séquence alignée et on la stocke dans une nouvelle liste
            align_list.append(traceback)
        return align_list

