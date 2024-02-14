import numpy as np


class Aligner:

	COMPOSITION = ["A", "C", "G", "T", "-"]

	def __init__(self, consensus, composition, nb_seq, kmer=200):
		self.kmer_size = kmer
		self.nb_sequences = nb_seq
		self.consensus = consensus
		self.composition = np.array(composition, dtype=int)
		self.kmer_consensus = self.kmers_from_sequence(self.consensus, self.kmer_size)

	@staticmethod
	def kmers_from_sequence(sequence, k, step=1):
		kmers = []
		for i in range(0, len(sequence), step):
			kmer = sequence[i:i + k]
			if len(kmer) == k:
				kmers.append(kmer)
		return set(kmers)

	def add_aligned_sequence(self, sequence):
		for position, nucleotide in enumerate(sequence):
			if nucleotide in self.COMPOSITION:  # Les N ne sont pas pris en compte
				self.composition[position, self.COMPOSITION.index(nucleotide)] += 1
		self.nb_sequences += 1

	def update_consensus(self):
		max_indices = self.composition.argmax(axis=1)
		self.consensus = ''.join(self.COMPOSITION[index] for index in max_indices)

	def prevalence(self, position):
		column = self.composition[position]
		total = column.sum()
		if total == 0:
			return None
		prevalence = (column.max() / total) * 100
		return prevalence

	def occupancy(self, position):
		column = self.composition[position]
		print(column[-1])
		nb_nuc_total = column.sum() - column[-1]
		occupancy = (nb_nuc_total / self.nb_sequences) * 100
		return occupancy

	def add_gap(self, position):
		gap_column = np.zeros((self.composition.shape[0]), dtype=int)
		gap_column[-1] = 1
		self.composition = np.insert(self.composition, position, gap_column, axis=0)

	def align_new_sequence(self, new_seq):
		kmer_new_seq = self.kmers_from_sequence(new_seq, self.kmer_size)
		common_kmers = self.kmer_consensus.intersection(kmer_new_seq)

		for kmer in common_kmers:
			start_consensus = self.consensus.find(kmer)
			start_new_seq = new_seq.find(kmer)


if __name__ == "__main__":
	test_cons = "ATG-C"
	test_comp = [[2, 0, 0, 0, 0],
	             [0, 0, 0, 2, 0],
	             [0, 0, 2, 0, 0],
	             [0, 0, 0, 0, 2],
	             [0, 2, 0, 0, 0]]
	test = Aligner(test_cons, test_comp, 2, kmer=2)
	print(test.consensus)
	print(test.composition)
	test.add_aligned_sequence("ACG-C")
	test.add_aligned_sequence("ACG-N")
	test.update_consensus()
	print(test.consensus)
	print(test.composition)
	test.add_gap(2)
	test.update_consensus()
	print(test.consensus)
	print(test.composition)
	print(test.prevalence(1))
	print(test.occupancy(4))
	print(test.kmer_consensus)
