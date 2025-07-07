"""
Author: Yesse Monkou
Date: July 7th 2025

This script mutates given SNPs onto a given sequence and checks the 
effect that it has on that sequence. This program identifies wether 
the mutation is: symomynous, nonsense or missense.

Access this script through command line using the following code:
python variant_predictor.py --fasta gene.fasta --variants variants.csv --output results.csv

- gene.fasta can contain as many sequences as you would like
- variant.csv must only contain: position,reference_base,alternative_base
"""

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import argparse


def read_fasta(file_path):
    """Reads in fasta file and checks if there is only 1 sequence

    Args:
        file_path (str): Path of fasta file

    Raises:
        ValueError: There are multiple sequences inside the fasta file

    Returns:
        str: The sequence
    """
    # Load in fasta file
    records = list(SeqIO.parse(file_path, "fasta"))
    record_dict = {}
    for record in records:
        record_dict[str(record.id)] = str(record.seq)
    
    return record_dict


def mutate_classify(sequence, variants, output, id):
    # Getting the original amino acid sequence for comparison
    original_aa_seq = str(Seq(sequence).translate())
    
    # Creating a list for mutated_sequence so the data is adjustable
    mutated_sequence = list(sequence)
    
    # Loops through every row in the variants dataframe
    for row in variants.itertuples(index=False):
        # Checks if the reference base from the variants dataframe is actually in the sequence
        if mutated_sequence[int(row.position)-1] == row.reference_base:
            # Mutates base
            mutated_sequence[int(row.position)-1] = row.alternative_base
            
            # Filling the output-dictionary
            output["fasta_header"].append(id)
            output["position"].append(row.position)
            output["ref_base"].append(row.reference_base)
            output["alt_base"].append(row.alternative_base)
            codon_index = (int(row.position)-1) // 3
            codon_start = codon_index * 3
            codon_end = codon_start + 3
            output["original_codon"].append(sequence[codon_start:codon_end])
            output["mutated_codon"].append("".join(mutated_sequence[codon_start:codon_end]))
            output["original_aa"].append(str(Seq(output["original_codon"][-1]).translate()))
            output["mutated_aa"].append(str(Seq(output["mutated_codon"][-1]).translate()))
            
            # For comparison to original aa sequence
            mutated_aa_seq = str(Seq("".join(mutated_sequence)).translate())
            
            # Checks 
            if mutated_aa_seq == original_aa_seq:
                output["mutation_type"].append("Synonynous")
            elif len(mutated_aa_seq) < len (original_aa_seq) or output["mutated_aa"][-1] == "*":
                output["mutation_type"].append("Nonsense")
            else:
                output["mutation_type"].append("Missense")
        
        # Reset mutation        
        mutated_sequence = list(sequence)

    return output

def main(fasta_path, variants_path, output_path):
    sequences = read_fasta(fasta_path)
    variants = pd.read_csv(variants_path)
    variants['reference_base'] = variants['reference_base'].str.upper()
    variants['alternative_base'] = variants['alternative_base'].str.upper()

    
    # Checks for any missing values
    if variants.isnull().values.any():
        raise ValueError("Your --variants file contains a missing value.")

    general_output = {"fasta_header": [],
              "position" : [],
              "ref_base": [],
              "alt_base": [],
              "original_codon": [],
              "mutated_codon": [],
              "original_aa": [],
              "mutated_aa": [],
              "mutation_type": []}
    for id in sequences:
        sequence = sequences[id].upper()
        general_output = mutate_classify(sequence, variants, general_output, id)
    
    # Write to CSV
    pd.DataFrame(general_output).to_csv(output_path, index=False)
    

if __name__ == "__main__":
    # To parse arguments
    parser = argparse.ArgumentParser(description="SNP Effect Predictor for Protein-Coding Genes")
    parser.add_argument('--fasta', required=True, help='Path to input FASTA file')
    parser.add_argument('--variants', required=True, help='Path to SNP variants CSV file)')
    parser.add_argument('--output', required=True, help='Path to output results CSV file')

    args = parser.parse_args()

    fasta_path = args.fasta
    variants_path = args.variants
    output_path = args.output
    
    main(fasta_path, variants_path, output_path)
