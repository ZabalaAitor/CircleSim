import os

class Join:
    def __init__(self, args):
        self.circle_fastq = args.circle_fastq
        self.linear_fastq = args.linear_fastq
        self.output_fastq = args.output_fastq

    # Read two input FASTQ files (circular and linear) and combines them into one output FASTQ.
    def join(self):
        try:
            with open(self.circle_fastq, 'r') as f1, open(self.linear_fastq, 'r') as f2:
                fastq_data = f1.read().strip() + '\n' + f2.read().strip()
            with open(self.output_fastq, 'w') as output:
                output.write(fastq_data)
        except FileNotFoundError as e:
            print(f"Error: {e}. One or both of the input FASTQ files not found.")
        except Exception as e:
            print(f"An error occurred: {e}")

def main():
    args = parse_arguments()
    fastq_joiner = Join(args)
    fastq_joiner.join()

if __name__ == "__main__":
    main()
