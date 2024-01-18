import argparse
import os

def parse_eccDNA_arguments():
    parser = argparse.ArgumentParser(description='CircSim eccDNA Utility Functions')
    parser.add_argument('-t', '--total_eccDNA', type=int, default=100, help='Total number of eccDNA to simulate')
    parser.add_argument('-l', '--length_min', type=int, default=300, help='Minimum length of eccDNA')
    parser.add_argument('-L', '--length_max', type=int, default=2000, help='Maximum length of eccDNA')
    parser.add_argument('-g', '--genome_fasta', type=str, required=True, help='Path to the genome FASTA file')
    parser.add_argument('-o', '--output_bed', type=str, default=os.getcwd(), help='Path to the output BED file for eccDNA')
    return parser.parse_args()

def parse_reads_arguments():
    parser = argparse.ArgumentParser(description='CircSim Reads Utility Functions')
    parser.add_argument('-s', '--sequence', type=str, default='short', help='Coverage for reads simulation')
    parser.add_argument('-c', '--coverage', type=float, default=30, help='Coverage for reads simulation')
    parser.add_argument('-r', '--reads_length', type=int, default=150, help='Length of simulated reads')
    parser.add_argument('-i', '--insert_length', type=int, default=500, help='Insert length for reads simulation')
    parser.add_argument('-g', '--genome_fasta', type=str, required=True, help='Path to the genome FASTA file')
    parser.add_argument('-b', '--input_bed', type=str, required=True, help='Path to the input BED file for eccDNA')
    parser.add_argument('-o', '--output_fastq', type=str, default=os.getcwd(), help='Path to the output FASTQ file for reads')
    return parser.parse_args()