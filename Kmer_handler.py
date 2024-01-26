from Extract_fasta import Fasta_extract 

#Cette classe permet d'extraire des kmer depuis une liste contenant les str de séquences. 

class kmer_handler : 

    def __init__(self, data) : 
        #La classe prends en entrée simplement la liste
        self.__data = data

    def __kmer(self, seq, k, inter) : 
        #Cette fonction privée permet, pour une séquence, de récupérer des kmers de la taille k, prenant les kmers avec un intervale de inter.
        #Les kmers sont stocks dans un dictionnaire dans lequel ils sont la clé, et la valeur est chacun de leur position.
        dict_kmer = {}
        for carac in range(0, len(seq), inter) :
            kmer = seq[carac:carac+k]
            if kmer not in dict_kmer : 
                dict_kmer[kmer] = [carac]
            else :
                dict_kmer[kmer].append(carac) 
        return dict_kmer

    def kmer_splicer(self, k, inter) : 
        #Cette fonction permet de créer une liste de dictionnaire, chaque dictionnaire contenant tout les kmers détectés dans une séquence, selon les paramètres k et inter (voir doc __kmer)
        #kmer_splicer utilise la fonction privé __kmer sur chaque séquence de la liste de séquence donnée en condition à la classe
        list_dict_kmer = []
        for sequence in self.__data :
            list_dict_kmer.append(self.__kmer(sequence, k, inter))
        return list_dict_kmer


extractor = Fasta_extract('sequences_SARS-CoV-2_ech.fasta')

data = extractor.sequences()

test = kmer_handler(data)

kmers = test.kmer_splicer(300, 2)

print(len(kmers[1]))
