import argparse
import os
import sys
import time
from simulate_eccDNA import SimulateEccDNA
from simulate_reads import SimulateReads

class CircSim:

    def __init__(self):
        self.cs_version = "1.0"  # Define cs_version here or retrieve it from your code

        self.parser = argparse.ArgumentParser(
            description='CircSim',
            usage=f'''CircSim <subprograms> [options]

version = {self.cs_version}

CircSim

Commands:

    eccDNA  Simulate eccDNA
    reads   Simulate eccDNA reads
'''
        )
        subparsers = self.parser.add_subparsers(dest='subcommand')

        self.eccDNA = subparsers.add_parser(
            name='eccDNA',
            description='Simulate eccDNA',
            prog='CircSim eccDNA',
            usage='''CircSim eccDNA [options]'''
        )

        self.reads = subparsers.add_parser(
            name='reads',
            description='Simulate eccDNA reads',
            prog='CircSim reads',
            usage='''CircSim reads [options]'''
        )

        if len(sys.argv) <= 1:
            self.parser.print_help()
            time.sleep(0.01)
            sys.stderr.write("\nNo argument given to CircSim")
            sys.exit(0)

        elif sys.argv[1] == 'eccDNA':
            self.subprogram = self.args_eccDNA()
            self.args = self.subprogram.parse_args(sys.argv[2:])

            eccdna_instance = eccDNA(
                total_eccDNA=self.args.total_eccDNA,
                length_min=self.args.length_min,
                length_max=self.args.length_max,
                genome_fasta=self.args.genome_fasta,
                output_bed=self.args.output_bed
            )
            eccdna_instance.run()

        elif sys.argv[1] == 'reads':
            self.subprogram = self.args_reads()
            self.args = self.subprogram.parse_args(sys.argv[2:])
            
            reads_instance = reads(
                sequence=self.args.sequence,
                coverage=self.args.coverage,
                reads_length=self.args.reads_length,
                insert_length=self.args.insert_length,
                genome_fasta=self.args.genome_fasta,
                input_bed=self.args.input_bed,
                output_fastq=self.args.output_fastq
            )
            reads_instance.run()

    def args_eccDNA(self):
        parser = self.eccDNA
        parser._action_groups.pop()
        parser.add_argument('-t', '--total_eccDNA', type=int, default=100, help='Total number of eccDNA to simulate')
        parser.add_argument('-l', '--length_min', type=int, default=300, help='Minimum length of eccDNA')
        parser.add_argument('-L', '--length_max', type=int, default=2000, help='Maximum length of eccDNA')
        parser.add_argument('-g', '--genome_fasta', type=str, required=True, help='Path to the genome FASTA file')
        parser.add_argument('-o', '--output_bed', type=str, default=os.getcwd(), help='Path to the output BED file for eccDNA')
        return parser

    def args_reads(self):
        parser = self.reads
        parser._action_groups.pop()
        parser.add_argument('-s', '--sequence', type=str, default='short', help='Coverage for reads simulation')
        parser.add_argument('-c', '--coverage', type=float, default=30, help='Coverage for reads simulation')
        parser.add_argument('-r', '--reads_length', type=int, default=150, help='Length of simulated reads')
        parser.add_argument('-i', '--insert_length', type=int, default=500, help='Insert length for reads simulation')
        parser.add_argument('-g', '--genome_fasta', type=str, required=True, help='Path to the genome FASTA file')
        parser.add_argument('-b', '--input_bed', type=str, required=True, help='Path to the input BED file for eccDNA')
        parser.add_argument('-o', '--output_fastq', type=str, default=os.getcwd(), help='Path to the output FASTQ file for reads')
        return parser

if __name__ == "__main__":
    circsim = CircSim()
