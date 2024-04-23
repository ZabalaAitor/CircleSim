import argparse
import os

def parse_coordinates_arguments(parser):
    parser.add_argument('-t', '--type', type=str, required=True, choices=['DNA', 'RNA'], help='Type of circle (DNA or RNA)')
    parser.add_argument('-n', '--number', type=int, default=100, help='Number of circles to simulate')
    parser.add_argument('-d', '--distribution', type=str, choices=['uniform', 'lognormal'], default='uniform', help='Distribution of circle length (uniform or lognormal) [default = uniform]')
    parser.add_argument('-l', '--length_min', type=int, default=300, help='Minimum length of eccDNA [default = 300]')
    parser.add_argument('-L', '--length_max', type=int, default=2000, help='Maximum length of eccDNA [default = 2000]')
    parser.add_argument('-m', '--mean', type=int, default=1000, help='Mean length of eccDNA [default = 1000]')
    parser.add_argument('-sd', '--sd', type=int, default=1, help='Standard deviation length of eccDNA [default = 1]')
    parser.add_argument('-g', '--genome_fasta', type=str, help='Path to the genome FASTA file')
    parser.add_argument('-r', '--transcript_bed', type=str, help='Path to the trasncript BED file')
    parser.add_argument('-o', '--output_bed', type=str, default=os.getcwd(), help='Path to the output BED file for circles')

def parse_reads_arguments(parser):
    parser.add_argument('-t', '--type', type=str, required=True, choices=['DNA', 'RNA', 'linear'], help='Type of circle (DNA or RNA)')
    parser.add_argument('-s', '--sequence', type=str, default='short', help='Coverage for reads simulation')
    parser.add_argument('-c', '--coverage', type=float, default=30, help='Coverage for reads simulation')
    parser.add_argument('-r', '--reads_length', type=int, default=150, help='Length of simulated reads')
    parser.add_argument('-i', '--insert_length', type=int, default=500, help='Insert length for reads simulation')
    parser.add_argument('-a', '--alpha', type=int, default=0.5, help='Alpha value for beta distribution [default = 0.5]')
    parser.add_argument('-v', '--beta', type=int, default=0.5, help='Beta value for beta distribution [default = 0.5]')
    parser.add_argument('-g', '--genome_fasta', type=str, help='Path to the genome FASTA file')
    parser.add_argument('-b', '--input_bed', type=str, required=True, help='Path to the input BED file for eccDNA')
    parser.add_argument('-o', '--output_fastq', type=str, default=os.getcwd(), help='Path to the output FASTQ file for reads')

def parse_join_arguments(parser):
    parser.add_argument("-c", "--circle_fastq", type=str, help="Path to the input fastq file for circle", required=True)
    parser.add_argument("-l", "--linear_fastq", type=str, help="Path to the second BED file", required=True)
    parser.add_argument("-o", "--output_fastq", type=str, default=os.getcwd(), help="Path to the output FASTQ file for reads")

def set_default_genome_fasta(args):
    default_genome_fasta = '/home/unidad/Descargas/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz'

    if os.path.exists(default_genome_fasta):
        args.genome_fasta = default_genome_fasta
    else:
        print("Default genome FASTA file not found locally. Downloading from GitHub...")
        subprocess.run(['wget', default_genome_fasta])
        args.genome_fasta = default_genome_fasta
