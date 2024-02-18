#Ce fichier permet d'extraire les séquences d'un fichier fasta sous forme d'une liste contrenant des string, chaque string étant une séquence.
#Ce procédé efface les noms des séquences.
import os
import argparse
import glob

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
        self.sequence_list = sequence_list[1:len(sequence_list)]
        return sequence_list[1:len(sequence_list)]

    def sequences_packages(self, pack_size = 1000) : #Cette fonction permet de récupérer les séquences et de les dispatcher dans un dossier "sequences_packages" en des fichiers contenant "pack_size" sequences
        if not os.path.exists("sequences_packages"):
            os.makedirs("sequences_packages") #On vérifie si le dossier d'acceuil existe, et si non, on le créé
        nbr_seq = 1
        nbr_pack = 0
        while nbr_seq-1 != len(self.sequence_list) : #Cette boucle permet de créer plusieurs fichiers, tant qu'on a pas atteint la fin de la liste.
            title = "sequences_packages/pack_" + str(nbr_pack) + ".fasta" #On définit le titre du fichier actuel
            with open(title, 'a') as fileOut : 
                while nbr_seq%pack_size != 0 : #Tant qu'on a pas remplis le fichier avec le nombre exact de séquences, on continue
                    fileOut.write(">seq n°" + str(nbr_seq) + "\n" + self.sequence_list[nbr_seq-1] + '\n\n')
                    nbr_seq += 1
                    if nbr_seq == len(self.sequence_list) : #Si on tombe sur la fin de seq_list à ce moment là, on arrête tout.
                        break 
                with open(title, 'a') as fileOut : #Pour éviter de sauter une séquence lorsque le fichier est plein, on y rajoute la dernière séquence (la dernière séquence a un rapport nbr_seq%pack_size de 0 )
                    fileOut.write(">seq n°" + str(nbr_seq) + "\n" + self.sequence_list[nbr_seq-1] + '\n\n') 
            nbr_seq += 1
            nbr_pack += 1

    def add_sequences(self) : 
        list_fichier = glob.glob('sequences_packages/pack_*.fasta') #Afin de pouvoir ajouter des séquences, on doit savoir combien de fichier on a déjà
        nom_fichier_base = list_fichier[0] #On va aussi regarder, dans le premier fichier la taille des paquets (/!\ si le paquet 0 n'est pas complet, le nombre de séquences dedans sera pris comme référence)
        with open(nom_fichier_base) as read :
            pack_size = 0
            for l in read : 
                if l[0] == '>' : 
                    pack_size += 1
        nbr_seq = 1 #À partir d'ici, le code est quasi identique à celui de sequences_packages.
        nbr_pack = len(list_fichier) #Avecle nombre de fichier on définit le nouveau point de départ de nos n° de fichier
        while nbr_seq-1 != len(self.sequence_list) : #Cette boucle permet de créer plusieurs fichiers, tant qu'on a pas atteint la fin de la liste.
            title = "sequences_packages/pack_" + str(nbr_pack) + ".fasta" #On définit le titre du fichier actuel
            with open(title, 'a') as fileOut : 
                while nbr_seq%pack_size != 0 : #Tant qu'on a pas remplis le fichier avec le nombre exact de séquences, on continue
                    fileOut.write(">seq n°" + str(nbr_seq) + "\n" + self.sequence_list[nbr_seq-1] + '\n\n')
                    nbr_seq += 1
                    if nbr_seq == len(self.sequence_list) : #Si on tombe sur la fin de seq_list à ce moment là, on arrête tout.
                        break 
                with open(title, 'a') as fileOut : #Pour éviter de sauter une séquence lorsque le fichier est plein, on y rajoute la dernière séquence (la dernière séquence a un rapport nbr_seq%pack_size de 0 )
                    fileOut.write(">seq n°" + str(nbr_seq) + "\n" + self.sequence_list[nbr_seq-1] + '\n\n') 
            nbr_seq += 1
            nbr_pack += 1

#file = Fasta_extract('sequences_20240202_4425212.fasta') #Test sur un fichier dont le compte de séquence n'est pas rond (237)
#file.sequences()
#file.sequences_packages(100)

#Pensez à supprimer les fichiers dans le dossier recevant entre les deux test.

#file = Fasta_extract('sequences_20240202_4425212.fasta') #Test sur un fichier dont le compte de séquence est rond (200)
#file.sequences()
#file.sequences_packages(100)
            
if __name__ == "__main__"  : 
    arg_manager = argparse.ArgumentParser()
    arg_manager.add_argument('-add_new_sequences', '-n', type = str)
    arg_manager.add_argument('-input', '-i', type = str)
    arg_manager.add_argument('-pack_size', '-ps', type = int)
    args = arg_manager.parse_args()
    if args.add_new_sequences == "y": 
        file = Fasta_extract(args.input)
        file.sequences()
        file.add_sequences()
    else : 
        file = Fasta_extract(args.input)
        file.sequences()
        file.sequences_packages(args.pack_size)