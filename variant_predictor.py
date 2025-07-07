from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import argparse


def read_fasta(file_path):
    records = list(SeqIO.parse("gene.fasta", "fasta"))
    if records[0].id != records[-1].id:
        raise ValueError("Multiple sequences detected. This tool supports single-gene FASTA files only.")
    
    return str(records[0].seq)


def mutate_classify(sequence, variants):
    original_aa_seq = str(Seq(sequence).translate())
    mutated_sequence = list(sequence)
    output = {"position" : [],
              "ref_base": [],
              "alt_base": [],
              "original_codon": [],
              "mutated_codon": [],
              "original_aa": [],
              "mutated_aa": [],
              "mutation_type": []}
    
    for row in variants.itertuples(index=False):
        if mutated_sequence[int(row.position)-1] == row.reference_base:
            mutated_sequence[int(row.position)-1] = row.alternative_base
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
            mutated_aa_seq = str(Seq("".join(mutated_sequence)).translate())
            
            if mutated_aa_seq == original_aa_seq:
                output["mutation_type"].append("Synomynous")
            elif len(mutated_aa_seq) < len (original_aa_seq) or output["mutated_aa"][-1] == "*":
                output["mutation_type"].append("Nonsense")
            else:
                output["mutation_type"].append("Missense")
                
        mutated_sequence = list(sequence)

    return pd.DataFrame(output)

def main():
    parser = argparse.ArgumentParser(description="SNP Effect Predictor for Protein-Coding Genes")
    parser.add_argument('--fasta', required=True, help='Path to input FASTA file (SEQUENCE MUST BE IN CAPITAL)')
    parser.add_argument('--variants', required=True, help='Path to SNP variants CSV file (BASES MUST BE IN CAPITAL)')
    parser.add_argument('--output', required=True, help='Path to output results CSV file')

    args = parser.parse_args()

    fasta_path = args.fasta
    variants_path = args.variants
    output_path = args.output

    sequence = read_fasta(fasta_path)
    variants = pd.read_csv(variants_path)

    output = mutate_classify(sequence, variants)
    output.to_csv(output_path, index=False)
    

if __name__ == "__main__":
    main()
    