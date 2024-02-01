from Extract_fasta import Fasta_extract

import numpy as np
import matplotlib.pyplot as plt


class SequenceAnalyser:
	def __init__(self, consensus_seq, list_sequences):
		self.consensus_seq = consensus_seq
		self.list_sequences = list_sequences

	@staticmethod
	def kmers_from_sequence(sequence, k):
		kmers = []
		for i in range(0, len(sequence)):
			kmer = sequence[i:i + k]
			if len(kmer) == k:
				kmers.append(kmer)
		return set(kmers)

	@staticmethod
	def taux_kmers_communs(kmers_sequence, kmers_consensus):
		kmers_communs = kmers_sequence.intersection(kmers_consensus)
		if len(kmers_sequence) == 0:
			return 0
		pourcentage = (len(kmers_communs) / len(kmers_sequence)) * 100
		return pourcentage

	def get_list_pourcentage_communs(self, list_k_sizes):
		k_pourc_list = []
		for k_size in list_k_sizes:
			print(f"K-size: {k_size}")
			kmers_consensus = self.kmers_from_sequence(self.consensus_seq, k_size)
			seq_pourc_list = []
			for sequence in self.list_sequences:
				kmers_sequence = self.kmers_from_sequence(sequence, k_size)
				pourcentage = self.taux_kmers_communs(kmers_sequence, kmers_consensus)
				seq_pourc_list.append(pourcentage)
			k_pourc_list.append(seq_pourc_list)
		return k_pourc_list

	def build_graph_pourcentage_communs(self, pourcentages_list, list_k_sizes, title):
		moyennes = [np.mean(p) for p in pourcentages_list]
		ecarts_types = [np.std(p) for p in pourcentages_list]

		plt.plot(list_k_sizes, moyennes)
		plt.fill_between(list_k_sizes, np.array(moyennes) - np.array(ecarts_types),
		                 np.array(moyennes) + np.array(ecarts_types), color='blue', alpha=0.2)

		plt.title("Moyenne des Pourcentages d'Identité par Taille de k-mer avec Écart Type")
		plt.xlabel("Taille de k-mer")
		plt.ylabel("Moyenne des Pourcentages d'Identité")

		plt.savefig(title)
		plt.close()


if __name__ == "__main__":
	consensus_seq = Fasta_extract("consensus.fasta").sequences()[0]
	list_sequences = Fasta_extract("sequences_SARS-CoV-2_ech.fasta").sequences()
	analyser = SequenceAnalyser(consensus_seq, list_sequences)
	list_k_sizes = [i for i in range(10, 1010, 10)]
	pourcentages = analyser.get_list_pourcentage_communs(list_k_sizes)
	analyser.build_graph_pourcentage_communs(pourcentages, list_k_sizes, "graphique_kmers_tous.png")
