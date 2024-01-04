from __future__ import division
from Bio.Data import IUPACData
import numpy as np
#from numpy import *
import argparse

#Karan

positive_pKs = {'Nterm': 7.5, 'K': 10.0, 'R': 12.0, 'H': 5.98}
negative_pKs = {'Cterm': 3.55, 'D': 4.05, 'E': 4.45, 'C': 9.0, 'Y': 10.0}
pKcterminal = {'D': 4.55, 'E': 4.75}
pKnterminal = {'A': 7.59, 'M': 7.0, 'S': 6.93, 'P': 8.36, 'T': 6.82, 'V': 7.44,
               'E': 7.7}
charged_aas = ('K', 'R', 'H', 'D', 'E', 'C', 'Y')

class Charge_Calculator(object):

    def __init__(self, ProteinSequence):
        """Initialize the class."""
        self.sequence = ProteinSequence
        self.amino_acids_content = self.count_amino_acids(self.sequence)
        self.charged_aas_content = self._select_charged(self.amino_acids_content)

    def count_amino_acids(self, seq):
        """Count standard amino acids, return a dict.

        Counts the number times each amino acid is in the protein
        sequence. Returns a dictionary {AminoAcid:Number}.

        The return value is cached in self.amino_acids_content.
        It is not recalculated upon subsequent calls.
        """
        #if self.amino_acids_content is None:
        prot_dic = dict((k, 0) for k in IUPACData.protein_letters)
        for aa in prot_dic:
            prot_dic[aa] = seq.count(aa)

            #self.amino_acids_content = prot_dic

        return prot_dic

    # This function creates a dictionary with the contents of each charged aa,
    # plus Cterm and Nterm.
    def _select_charged(self, AminoAcidsContent):
        charged = {}
        for aa in charged_aas:
            charged[aa] = float(AminoAcidsContent[aa])
        charged['Nterm'] = 1.0
        charged['Cterm'] = 1.0
        return charged

    # This function calculates the total charge of the protein at a given pH.
    def _chargeR(self, pH, pos_pKs, neg_pKs):
        PositiveCharge = 0.0
        for aa, pK in pos_pKs.items():
            CR = 10 ** (pK - pH)
            partial_charge = CR / (CR + 1.0)
            PositiveCharge += self.charged_aas_content[aa] * partial_charge

        NegativeCharge = 0.0
        for aa, pK in neg_pKs.items():
            CR = 10 ** (pH - pK)
            partial_charge = CR / (CR + 1.0)
            NegativeCharge += self.charged_aas_content[aa] * partial_charge

        return PositiveCharge - NegativeCharge

    def charge_values(self, pH_x):

        pos_pKs = dict(positive_pKs)
        neg_pKs = dict(negative_pKs)
        nterm = self.sequence[0]
        cterm = self.sequence[-1]
        if nterm in pKnterminal:
            pos_pKs['Nterm'] = pKnterminal[nterm]
        if cterm in pKcterminal:
            neg_pKs['Cterm'] = pKcterminal[cterm]

        charge_y = []


        for pH in pH_x:
            charge = self._chargeR(pH, pos_pKs, neg_pKs)
            charge_y.append(charge)

        return charge_y

    def charge_curve(self, name):

        pH_x = np.arange(1,14,0.1)

        charge_y = self.charge_values(pH_x)

        np.savetxt(str(name)+'.dat', np.column_stack((pH_x, charge_y)), fmt='%f')


def main():
    parser = argparse.ArgumentParser(description='Calculate the charge curve of a protein sequence.')
    parser.add_argument('sequence', type=str, help='Protein sequence to calculate the charge curve for')
    parser.add_argument('-o', '--output', type=str, default='charge_curve', help='Output file name without extension')
    args = parser.parse_args()

    protein_sequence = args.sequence
    output_file = args.output

    # Ensure the sequence only contains valid amino acids
    valid_aas = set(IUPACData.protein_letters)
    if not set(protein_sequence).issubset(valid_aas):
        raise ValueError("Invalid protein sequence: contains non-standard amino acids.")

    # Create an instance of Charge_Calculator and calculate the charge curve
    charge_calculator = Charge_Calculator(protein_sequence)
    charge_calculator.charge_curve(output_file)
    print(f"Charge curve data saved to {output_file}.dat")

if __name__ == "__main__":
    main()
