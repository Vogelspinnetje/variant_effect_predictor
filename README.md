# 🧬 SNP Effect Predictor for Protein-Coding Genes

## Goal:
Build a Python program that predicts the effect of single nucleotide polymorphisms (SNPs) on protein-coding DNA sequences. The program identifies whether the mutation is:
- Synonymous (no amino acid change)
- Missense (amino acid change)
- Nonsense (introduces stop codon)

## Test Dataset and Example Query:
### FASTA File (coding DNA sequence):
File: `gene.fasta`
```
>EXAMPLE_GENE
ATGGAGGAGCTGGTGGTGACGTGGACCTGAAG
```

### Variant File (CSV):
File: `variants.csv`
```
position,reference_base,alternative_base
4,G,A
15,G,A
31,C,G
```

### Expected Outcome Example:
```
position,ref_base,alt_base,original_codon,mutated_codon,original_aa,mutated_aa,mutation_type
4,G,A,GAG,AAG,E,K,Missense
15,G,A,GTG,ATG,V,M,Missense
31,C,G,TGA,TGG,*,W,Nonsense Rescue
```

## Program Description:
Your program should:
- ✅ Read a coding DNA sequence from a FASTA file
- ✅ Read a list of SNPs from a CSV file
- ✅ Simulate each SNP in the DNA sequence
- ✅ Translate both the original and mutated DNA to protein
- ✅ Compare the proteins to classify the mutation (Synonymous, Missense, Nonsense)
- ✅ Output a report (CSV or console) summarizing each SNP's effect

# Bonus (optional):
- Validate SNPs (correct position, matching reference base)
- Handle errors gracefully

## Deliverables:
- Python script: variant_predictor.py
- Output file with results
- README.md explaining:
    > What the program does
    > How to run it
    > Example command:
      `python variant_predictor.py --fasta gene.fasta --variants variants.csv --output results.csv`
    > Clean, readable code with comments
