# Genetic-Variant-Classifier

A Python-based tool for classifying genetic variants into **non-coding, synonymous, and non-synonymous** categories. It processes **VCF, GFF, and FASTA** files, applies quality filters, and outputs classified variants along with a visualization of variant distribution.

## Features

1.Parses **VCF files** to extract genomic variants  
2.Uses **GFF annotations** to determine coding regions  
3.Reads **FASTA sequences** to analyze protein-coding effects  
4.Classifies variants into **non-coding, synonymous, and non-synonymous**  
5.Filters out low-quality variants (QUAL < 20)  
6.Logs errors and saves results in structured output files  
7.Generates a **visualization plot** for variant distribution  

## Running the Script:

python3 GenVarClassifier.py --vcf "path to variants.vcf" --gff "path to annotations.gff" --fasta "path to genome.fasta"
