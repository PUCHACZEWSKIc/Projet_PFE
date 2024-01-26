from Extract_fasta import Fasta_extract 


class kmer_handler : 

    def __init__(self, data) : 
        self.__data = data

    def kmer(self, seq, k, inter) : 
        dict_kmer = {}
        for carac in range(0, len(seq), inter) :
            kmer = seq[carac:carac+k]
            if kmer not in dict_kmer : 
                dict_kmer[kmer] = [carac]
            else :
                dict_kmer[kmer].append(carac) 
        return dict_kmer

    def kmer_splicer(self, k, inter) : 
        list_dict_kmer = []
        for sequence in self.__data :
            list_dict_kmer.append(self.kmer(sequence, k, inter))
        return list_dict_kmer


extractor = Fasta_extract('sequences_SARS-CoV-2_ech.fasta')

data = extractor.sequences()

test = kmer_handler(data)

kmers = test.kmer_splicer(300, 2)

print(len(kmers[1]))
