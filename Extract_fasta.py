#Ce fichier permet d'extraire les séquences d'un fichier fasta sous forme d'une liste contrenant des string, chaque string étant une séquence.
#Ce procédé efface les noms des séquences.

class Fasta_extract : 
    
    def __init__(self, file) :
        #La classe ne prends en entrée que le nom du fichier fasta dont on veux extraire les séquences. 
        self.__file = file

    def sequences(self) : 
        sequence_list = []
        act_line = ''
        with open(self.__file) as read : 
            for seq in read : 
                seq = seq.strip("\n") #On retire les sauts de ligne
                if seq == "" : 
                    pass
                elif seq[0] == '>' : #Ici, "if" et "elif" on pour but de passer les lignes vides et les intitulés de séquences. 
                    sequence_list.append(act_line) #Si l'on trouve un intitulé de séquence, c'est que la séquence précédente est terminée, on l'ajoute donc à notre liste finale.
                    act_line = ''
                else :
                    act_line += seq
        sequence_list.append(act_line) 
        return sequence_list[1:len(sequence_list)]


