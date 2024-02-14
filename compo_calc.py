import argparse
import numpy as np
from Bio import SeqIO


class CompositionCalculator:
	"""
    Calcule la composition en nucléotides pour chaque position d'un alignement de séquences et génère une séquence
    consensus.

    Attributes
    ----------
    consensus : str
        La séquence consensus calculée à partir de l'alignement.
    composition : numpy.ndarray
        Matrice de composition des nucléotides à chaque position de l'alignement.
    """

	COMPOSITION = ["A", "C", "G", "T", "N", "-"]

	def __init__(self):
		"""
		Initialise l'analyseur de composition avec des attributs par défaut.
		"""
		self.consensus = ""
		self.composition = None

	def add_sequence(self, sequence):
		"""
        Ajoute une séquence à la matrice de composition globale.

        Parameters
        ----------
        sequence : str
            Une séquence nucléotidique à ajouter à la composition globale.
        """
		if self.composition is None:
			self.composition = np.zeros((len(sequence), len(self.COMPOSITION)))
		for position, nucleotide in enumerate(sequence):
			if nucleotide.upper() in self.COMPOSITION:
				# Les lettres caractérisant plusieurs nucléotides ne sont pas pris en compte
				self.composition[position, self.COMPOSITION.index(nucleotide.upper())] += 1

	def compute_consensus(self):
		"""
		Calcule la séquence consensus à partir de la matrice de composition actuelle.
		"""
		self.consensus = ''.join([self.COMPOSITION[col.argmax()] for col in self.composition])

	def export_csv(self, csv_outfile):
		"""
        Exporte la matrice de composition transposée dans un fichier CSV.

        Parameters
        ----------
        csv_outfile : str
            Chemin du fichier CSV pour exporter la composition.
        """
		with open(csv_outfile, 'ab') as file_out:  # Mode append binaire
			np.savetxt(file_out, self.composition.T, delimiter=',', fmt='%d')   # Format entier décimal pour les chiffres

	def append_fasta(self, fasta_infile, fasta_outfile):
		"""
        Ajoute la séquence consensus à un fichier FASTA existant.

        Parameters
        ----------
        fasta_infile : str
            Nom du fichier FASTA d'entrée utilisé pour générer la séquence consensus.
        fasta_outfile : str
            Chemin du fichier FASTA pour append la séquence consensus.
        """
		with open(fasta_outfile, 'a') as file_out:
			file_out.write(f">{fasta_infile}\n{self.consensus}\n")


def main(fasta_infile, csv_outfile, fasta_outfile):
	"""
    Fonction principale pour analyser un fichier FASTA et produire la matrice de composition et la séquence consensus.

    Parameters
    ----------
    fasta_infile : str
        Chemin du fichier multifasta d'entrée.
    csv_outfile : str
        Chemin du fichier CSV de sortie pour la composition.
    fasta_outfile : str
        Chemin du fichier FASTA de sortie pour la séquence consensus.
    """
	analyser = CompositionCalculator()

	for record in SeqIO.parse(fasta_infile, "fasta"):
		analyser.add_sequence(str(record.seq))

	analyser.compute_consensus()
	analyser.export_csv(csv_outfile)
	analyser.append_fasta(fasta_infile, fasta_outfile)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Produit une matrice de composition et une séquence consensus à "
	                                             "partir d'un alignement multiple")
	parser.add_argument("fasta_infile", type=str, help="Le fichier multifasta d'entrée")
	parser.add_argument("csv_outfile", type=str, help="Le fichier CSV de sortie pour la composition")
	parser.add_argument("fasta_outfile", type=str,
	                    help="Le fichier FASTA de sortie pour la séquence consensus")

	args = parser.parse_args()

	main(args.fasta_infile, args.csv_outfile, args.fasta_outfile)

