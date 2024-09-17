# CircleSim

**CircleSim** is a command-line tool designed to simulate circular DNA and RNA sequences and generate corresponding sequencing reads. This tool can be helpful for researchers working with circular DNA (like eccDNA) or RNA, providing an easy way to simulate the coordinates of these molecules and generate sequencing files for further analysis.

## Features

- **Simulate Coordinates**: Generates genomic coordinates for circular DNA or RNA molecules with customizable parameters.
- **Simulate Reads**: Creates simulated sequencing reads for the generated circular sequences, including options to introduce mutations and sequencing errors.
- **Merge FASTQ Files**: Merges FASTQ files from different sequencing simulations into a single output.

## Installation

To install CircleSim in a Conda environment, follow these steps:

```bash
git clone https://github.com/yourusername/CircleSim.git
cd CircleSim
conda create -n circlesim python=3.8
conda activate circlesim
pip install -r requirements.txt
```

## Data Requirements

CircleSim requires a reference genome in FASTA format for simulating circular DNA or RNA sequences. You can download a suitable genome file and place it in a designated folder called database. If a reference genome is not specified in the scripts, CircleSim will automatically look for the genome file in this folder. For circRNA simulations, circles can also be created from transcriptome regions using the ```database/Homo_sapiens.GRCh38.cdna.all.short.bed``` file.

## Using CircleSim

**Simulate Coordinates:**

```bash
python CircleSim.py coordinates -t DNA -T circular -n 100 -l 500 -L 1500 -o circular_coordinates.bed
python CircleSim.py coordinates -t DNA -T linear -n 100 -l 500 -L 1500 -o linear_coordinates.bed
```

**Simulate reads:**

```bash
python CircleSim.py reads -b circular_coordinates.bed -o circular_reads.fastq
python CircleSim.py reads -b linear_coordinates.bed -o linear_reads.fastq
```

**Merge FASTQ Files:**

```bash
python CircleSim.py  -c circular_reads.fastq -l linear_reads.fastq -o reads.fastq
```
## Parameters

### 1. Simulating Coordinates

The following parameters are used to generate genomic coordinates for circular DNA or RNA molecules:

- `-t`, `--type`: Type of molecule to simulate. Options: `DNA`, `RNA`. Default is `DNA`.
- `-T`, `--molecule`: Specifies whether the molecule is `linear` or `circular`. Default is `circular`.
- `-n`, `--number`: Number of circular sequences to simulate. Default is `100`.
- `-d`, `--distribution`: Distribution of the circle lengths. Options: `uniform`, `lognormal`. Default is `uniform`.
- `-l`, `--length_min`: Minimum length of the circular DNA/RNA. Default is `300`.
- `-L`, `--length_max`: Maximum length of the circular DNA/RNA. Default is `10,000`.
- `-m`, `--mean`: Mean length of the circles when using lognormal distribution. Default is `1,000`.
- `-sd`, `--sd`: Standard deviation for the lognormal distribution of circle lengths. Default is `1`.
- `-s`, `--split`: Specifies a sequence (e.g., `AGGT`) to enable splitting of circular sequences at specific points.
- `-g`, `--genome_fasta`: Path to the reference genome FASTA file.
- `-r`, `--transcript_bed`: Path to the transcript BED file for circRNA simulations.
- `-o`, `--output_bed`: Path for saving the output BED file with simulated coordinates. Default is the current directory.

### 2. Simulating Reads

These parameters are used to generate sequencing reads for the simulated circles:

- `-t`, `--type`: Type of circle. Options: `DNA` or `RNA` **Required**
- `-T`, `--molecule`: Specifies the molecule as either `linear` or `circular`. Default is `circular`.
- `-s`, `--sequence`: Coverage type for read simulation. Default is `short`.
- `-c`, `--coverage`: Sequencing coverage. Default is `30x`.
- `-r`, `--reads_length`: Length of the simulated reads. Default is `150` base pairs.
- `-i`, `--insert_length`: Insert length for paired-end reads. Default is `500`.
- `-a`, `--alpha`: Alpha parameter for the beta distribution, which influences mutation rates. Default is `0.5`.
- `-v`, `--beta`: Beta parameter for the beta distribution. Default is `0.5`.
- `--mutation`: Enables mutation simulation.
- `--save_unmutated`: Option to save unmutated reads.
- `--mutation_rate`: Specifies the mutation rate for the simulated reads. Default is `0.01`.
- `--error_rate`: Specifies the sequencing error rate. Default is `0.001`.
- `-g`, `--genome_fasta`: Path to the genome FASTA file.
- `-b`, `--input_bed`: Path to the input BED file with coordinates for the circular sequences. **Required**
- `-o`, `--output_fastq`: Path for saving the simulated FASTQ reads. Default is the current directory.

### 3. Merging FASTQ Files

Use these parameters to merge FASTQ files from different simulations:

- `-c`, `--circle_fastq`: Path to the FASTQ file for circular sequences. **Required**
- `-l`, `--linear_fastq`: Path to the FASTQ file for linear sequences. **Required**
- `-o`, `--output_fastq`: Path for saving the merged FASTQ file. Default is the current directory.

## License

CircleSim is freely available under the [MIT license](./LICENSE).

## Acknowledgements

CircleSim is developed by Aitor Zabala, Iñigo Prada-Luengo, Alex Martínez Ascensión and David Otaegui at Biogipuzkoa Health Research Institute.