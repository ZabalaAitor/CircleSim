import subprocess
import os

# Define the path for the CircleSim.py script
CIRCLESIM_PATH = '/data/CircleSim/src/CircleSim.py'

# Define the path for the genome fasta file
GENOME_FASTA_PATH = '/data/CircleSim/database/genome.fasta'

# Define the path for the transcript BED file
TRANSCRIPT_BED_PATH = '/data/CircleSim/database/Homo_sapiens.GRCh38.cdna.all.short.bed'

# Define the base path for output files
OUTPUT_BASE_PATH = '/data/CircleSim/data'

def main():
    coverages = [5, 7, 10, 15, 20, 30, 50, 70, 100]

    for circle_type in ['DNA', 'RNA']:
        # Step 1: Generate coordinates for circular and linear DNA/RNA once
        generate_coordinates(circle_type, circular=True)
        generate_coordinates(circle_type, circular=False)
        
        for coverage in coverages:
            # Step 2: Generate reads for circular and linear DNA/RNA
            generate_reads(circle_type, coverage, circular=True)
            generate_reads(circle_type, coverage, circular=False)
            
            # Step 3: Join the R1 circular with the linear and the R2 circular with the linear
            join_reads(circle_type, coverage)

def generate_coordinates(circle_type, circular):
    min_length = 300 if circular else 501
    file_suffix = "circular" if circular else "linear"
    output_path = os.path.join(OUTPUT_BASE_PATH, f'{file_suffix}{circle_type}.bed')
    print(f"Generating {file_suffix} coordinates for {circle_type}")
    subprocess.run([
        'python', CIRCLESIM_PATH, 'coordinates',
        '-t', circle_type,
        '-n', '1000',  # Example number, adjust as needed
        '-d', 'lognormal',
        '-l', str(min_length),
        '-L', '10000',
        '-m', '1000',
        '-sd', '1',
        '-g', GENOME_FASTA_PATH,
        '-r', TRANSCRIPT_BED_PATH,
        '-o', output_path
    ])

def generate_reads(circle_type, coverage, circular):
    file_suffix = "circular" if circular else "linear"
    reads_type = "linear" if not circular else circle_type
    coordinates_path = os.path.join(OUTPUT_BASE_PATH, f'{file_suffix}{circle_type}.bed')
    output_base_path = os.path.join(OUTPUT_BASE_PATH, f'{file_suffix}{circle_type}_cov{coverage}_')
    print(f"Generating {file_suffix} reads for {circle_type} with coverage {coverage}")
    subprocess.run([
        'python', CIRCLESIM_PATH, 'reads',
        '-t', reads_type,
        '-c', str(coverage),
        '-r', '150',
        '-i', '500',
        '-a', '0.5',
        '-v', '0.5',
        '-g', GENOME_FASTA_PATH,
        '-b', coordinates_path,
        '-o', output_base_path
    ])

def join_reads(circle_type, coverage):
    circular_r1_path = os.path.join(OUTPUT_BASE_PATH, f'circular{circle_type}_cov{coverage}_R1.fastq')
    linear_r1_path = os.path.join(OUTPUT_BASE_PATH, f'linear{circle_type}_cov{coverage}_R1.fastq')
    output_r1_path = os.path.join(OUTPUT_BASE_PATH, f'{circle_type}_cov{coverage}_R1.fastq')

    circular_r2_path = os.path.join(OUTPUT_BASE_PATH, f'circular{circle_type}_cov{coverage}_R2.fastq')
    linear_r2_path = os.path.join(OUTPUT_BASE_PATH, f'linear{circle_type}_cov{coverage}_R2.fastq')
    output_r2_path = os.path.join(OUTPUT_BASE_PATH, f'{circle_type}_cov{coverage}_R2.fastq')

    print(f"Joining reads for {circle_type} with coverage {coverage}")
    subprocess.run([
        'python', CIRCLESIM_PATH, 'join',
        '-c', circular_r1_path,
        '-l', linear_r1_path,
        '-o', output_r1_path
    ])
    subprocess.run([
        'python', CIRCLESIM_PATH, 'join',
        '-c', circular_r2_path,
        '-l', linear_r2_path,
        '-o', output_r2_path
    ])

if __name__ == "__main__":
    main()
