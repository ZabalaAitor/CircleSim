import subprocess
import os

# Define the path for the CircleSim.py script
CIRCLESIM_PATH = 'src/CircleSim.py'

# Define the path for the genome fasta file
GENOME_FASTA_PATH = 'database/genome.fasta'

# Define the path for the transcript BED file
TRANSCRIPT_BED_PATH = 'database/Homo_sapiens.GRCh38.cdna.all.short.bed'

# Define the base path for output files
OUTPUT_BASE_PATH = 'data/circDNA/short'

def main():
    coverages = [5, 7, 10, 15, 20, 30, 50, 70, 100]

    for circle_type in ['DNA']:
        for molecule_type in ['circular', 'linear']:
            generate_coordinates(circle_type, molecule_type)

    for circle_type in ['DNA']:
        for coverage in coverages:
            generate_reads(circle_type, coverage, 'circular')
            generate_reads(circle_type, coverage, 'linear')
            join_reads(circle_type, coverage)

def generate_coordinates(circle_type, molecule_type):
    min_length = 175 if molecule_type == 'circular' else 501
    max_length = 425 if molecule_type == 'circular' else 1000
    mean_length = 300 if molecule_type == 'circular' else 750
    output_path = os.path.join(OUTPUT_BASE_PATH, f'{molecule_type}{circle_type}.bed')
    print(f"Generating {molecule_type} coordinates for {circle_type}")
    subprocess.run([
        'python', CIRCLESIM_PATH, 'coordinates',
        '-t', circle_type,
        '-T', molecule_type,
        '-n', '1000', 
        '-d', 'lognormal',
        '-l', str(min_length),
        '-L', str(max_length),
        '-m', str(mean_length),
        '-sd', '1',
        '-g', GENOME_FASTA_PATH,
        '-r', TRANSCRIPT_BED_PATH,
        '-o', output_path
    ])

def generate_reads(circle_type, coverage, molecule_type):
    coordinates_path = os.path.join(OUTPUT_BASE_PATH, f'{molecule_type}{circle_type}.bed')
    output_base_path = os.path.join(OUTPUT_BASE_PATH, f'{molecule_type}{circle_type}_cov{coverage}_')
    print(f"Generating {molecule_type} reads for {circle_type} with coverage {coverage}")
    subprocess.run([
        'python', CIRCLESIM_PATH, 'reads',
        '-t', circle_type,
        '-T', molecule_type,
        '-c', str(coverage),
        '-r', '150',
        '-i', '500',
        '-a', '0.5',
        '-v', '0.5',
        '-g', GENOME_FASTA_PATH,
        '-b', coordinates_path,
        '-o', output_base_path,
        '--mutation',
        '--save_unmutated'
    ])

def join_reads(circle_type, coverage):
    # Paths for mutated reads
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

    # Paths for unmutated reads
    circular_unmutated_r1_path = os.path.join(OUTPUT_BASE_PATH, f'circular{circle_type}_cov{coverage}_unmutated_R1.fastq')
    linear_unmutated_r1_path = os.path.join(OUTPUT_BASE_PATH, f'linear{circle_type}_cov{coverage}_unmutated_R1.fastq')
    output_unmutated_r1_path = os.path.join(OUTPUT_BASE_PATH, f'{circle_type}_cov{coverage}_unmutated_R1.fastq')

    circular_unmutated_r2_path = os.path.join(OUTPUT_BASE_PATH, f'circular{circle_type}_cov{coverage}_unmutated_R2.fastq')
    linear_unmutated_r2_path = os.path.join(OUTPUT_BASE_PATH, f'linear{circle_type}_cov{coverage}_unmutated_R2.fastq')
    output_unmutated_r2_path = os.path.join(OUTPUT_BASE_PATH, f'{circle_type}_cov{coverage}_unmutated_R2.fastq')

    if os.path.exists(circular_unmutated_r1_path) and os.path.exists(linear_unmutated_r1_path):
        print(f"Joining unmutated reads for {circle_type} with coverage {coverage}")
        subprocess.run([
            'python', CIRCLESIM_PATH, 'join',
            '-c', circular_unmutated_r1_path,
            '-l', linear_unmutated_r1_path,
            '-o', output_unmutated_r1_path
        ])
        subprocess.run([
            'python', CIRCLESIM_PATH, 'join',
            '-c', circular_unmutated_r2_path,
            '-l', linear_unmutated_r2_path,
            '-o', output_unmutated_r2_path
        ])

if __name__ == "__main__":
    main()
