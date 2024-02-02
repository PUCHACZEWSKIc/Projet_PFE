#Ce fichier permet d'extraire les séquences d'un fichier fasta sous forme d'une liste contrenant des string, chaque string étant une séquence.
#Ce procédé efface les noms des séquences.
import os

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
        self.__sequence_list = sequence_list[1:len(sequence_list)]
        return sequence_list[1:len(sequence_list)]

    def sequences_packages(self, pack_size) : #Cette fonction permet de récupérer les séquences et de les dispatcher dans un dossier "sequences_packages" en des fichiers contenant "pack_size" sequences
        if not os.path.exists("sequences_packages"):
            os.makedirs("sequences_packages") #On vérifie si le dossier d'acceuil existe, et si non, on le créé
        nbr_seq = 1
        nbr_pack = 0
        while nbr_seq-1 != len(self.__sequence_list) : #Cette boucle permet de créer plusieurs fichiers, tant qu'on a pas atteint la fin de la liste.
            title = "sequences_packages/pack_" + str(nbr_pack) + ".fasta" #On définit le titre du fichier actuel
            with open(title, 'a') as fileOut : 
                while nbr_seq%pack_size != 0 : #Tant qu'on a pas remplis le fichier avec le nombre exact de séquences, on continue
                    fileOut.write(">seq n°" + str(nbr_seq) + "\n" + self.__sequence_list[nbr_seq-1] + '\n\n')
                    nbr_seq += 1
                    if nbr_seq == len(self.__sequence_list) : #Si on tombe sur la fin de seq_list à ce moment là, on arrête tout.
                        break 
                with open(title, 'a') as fileOut : #Pour éviter de sauter une séquence lorsque le fichier est plein, on y rajoute la dernière séquence (la dernière séquence a un rapport nbr_seq%pack_size de 0 )
                    fileOut.write(">seq n°" + str(nbr_seq) + "\n" + self.__sequence_list[nbr_seq-1] + '\n\n') 
            nbr_seq += 1
            nbr_pack += 1

file = Fasta_extract('sequences_20240202_4425212.fasta') #Test sur un fichier dont le compte de séquence n'est pas rond (237)
file.sequences()
file.sequences_packages(100)

#Pensez à supprimer les fichiers dans le dossier recevant entre les deux test.

file = Fasta_extract('sequences_20240202_4425212.fasta') #Test sur un fichier dont le compte de séquence est rond (200)
file.sequences()
file.sequences_packages(100)