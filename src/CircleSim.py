import subprocess
import argparse
import os
import sys
from coordinates import SimulateCoordinates 
from reads import SimulateReads
from reads import read_bed_file
from join import Join
from utils import parse_coordinates_arguments, parse_reads_arguments, parse_join_arguments


class CircleSim:

    def __init__(self):
        self.cs_version = "1.0"

        self.parser = argparse.ArgumentParser(
            description='CircleSim',
            usage=f'''CircleSim <subprograms> [options]

version = {self.cs_version}

CircleSim

Commands:

    coordinates  Simulate coordinates
    reads   Simulate reads 
    join   Join fastq files 
'''
        )
        subparsers = self.parser.add_subparsers(dest='subcommand')
        self.coordinates = subparsers.add_parser(
            name='coordinates',
            description='Simulate coordinates',
            prog='CircleSim coordinates',
            usage='''CircleSim coordinates [options]'''
        )
        self.coordinates.add_argument('-t', '--type', type=str, required=True, choices=['DNA', 'RNA'], help='Type of circle (DNA or RNA)')
        self.coordinates.add_argument('-n', '--number', type=int, default=100, help='Number of circles to simulate')
        self.coordinates.add_argument('-d', '--distribution', type=str, choices=['uniform', 'lognormal'], default='uniform', help='Distribution of circle length (uniform or lognormal [default = uniform])')
        self.coordinates.add_argument('-l', '--length_min', type=int, default=300, help='Minimum length of circle [default = 300]')
        self.coordinates.add_argument('-L', '--length_max', type=int, default=2000, help='Maximum length of circle [default = 2000]')
        self.coordinates.add_argument('-m', '--mean', type=int, default=1000, help='Mean length of circle [default = 1000]')
        self.coordinates.add_argument('-sd', '--sd', type=int, default=1, help='Standard deviation length of circle [default = 1]')
        self.coordinates.add_argument('-g', '--genome_fasta', type=str, help='Path to the genome FASTA file')
        self.coordinates.add_argument('-r', '--transcript_bed', type=str, help='Path to the trasncript BED file')
        self.coordinates.add_argument('-o', '--output_bed', type=str, default=os.getcwd(), help='Path to the output BED file for circles')

        self.reads = subparsers.add_parser(
            name='reads',
            description='Simulate circle reads',
            prog='CircleSim reads',
            usage='''CircleSim reads [options]'''
        )
        self.reads.add_argument('-t', '--type', type=str, required=True, choices=['DNA', 'RNA', 'linear'], help='Type of circle (DNA or RNA)')
        self.reads.add_argument('-s', '--sequence', type=str, default='short', help='Coverage for reads simulation')
        self.reads.add_argument('-c', '--coverage', type=float, default=30, help='Coverage for reads simulation [default = 30]')
        self.reads.add_argument('-r', '--reads_length', type=int, default=150, help='Length of simulated reads [default = 150]')
        self.reads.add_argument('-i', '--insert_length', type=int, default=500, help='Insert length for reads simulation [default = 500]')
        self.reads.add_argument('-a', '--alpha', type=int, default=0.5, help='Alpha value for beta distribution [default = 0.5]')
        self.reads.add_argument('-v', '--beta', type=int, default=0.5, help='Beta value for beta distribution [default = 0.5]')
        self.reads.add_argument('-g', '--genome_fasta', type=str, required=True, help='Path to the genome FASTA file')
        self.reads.add_argument('-b', '--input_bed', type=str, help='Path to the input BED file for eccDNA')
        self.reads.add_argument('-o', '--output_fastq', type=str, default=os.getcwd(), help='Path to the output FASTQ file for reads')

        self.join = subparsers.add_parser(
            name='join',
            description='Join fastq files',
            prog='CircleSim join',
            usage='''CircleSim join [options]'''
        )
        self.join.add_argument('-c', '--circle_fastq', type=str, help='Path to the input fastq file for circle')
        self.join.add_argument('-l', '--linear_fastq', type=str, help='Path to the input fastq file for linear')
        self.join.add_argument('-o', '--output_fastq', type=str, default=os.getcwd(), help='Path to the output FASTQ file for reads')

    def args_circles(self):
        return self.parser.parse_args()

    def args_reads(self):
        return self.parser.parse_args()

    def set_default_genome_fasta(args):
        default_genome_fasta = ''
        # Set default genome FASTA file based on the circle type
        default_genome_fasta = '/home/unidad/Descargas/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz'
    
        # Check if default genome FASTA file exists locally
        if os.path.exists(default_genome_fasta):
            args.genome_fasta = default_genome_fasta
        else:
            print("Default genome FASTA file not found locally. Downloading from GitHub...")
            # Download default genome FASTA file
            subprocess.run(['wget', default_genome_fasta])
            args.genome_fasta = default_genome_fasta

    def main(self):
        args = self.parser.parse_args()

        if args.subcommand == 'coordinates':
            # Check if genome FASTA file is not provided
            if not args.genome_fasta:
                set_default_genome_fasta(args)

            circles_instance = SimulateCoordinates(self.args_circles())
            circles_instance.run()
        elif args.subcommand == 'reads':
            # Check if genome FASTA file is not provided
            if not args.genome_fasta:
                set_default_genome_fasta(args)
   
            circle_bed = read_bed_file(args.input_bed)
            reads_instance = SimulateReads(args, circle_bed)
            reads_instance.simulate_reads()
        elif args.subcommand == 'join':
            join_instance = Join(args)
            join_instance.join()

if __name__ == "__main__":
    circlesim = CircleSim()
    circlesim.main()

