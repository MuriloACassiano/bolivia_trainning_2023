import gzip
import argparse

def separate_r1_r2(input_filename, output_r1, output_r2):
    with gzip.open(input_filename, 'rt') as f:
        with gzip.open(output_r1, 'wt') as r1_out, gzip.open(output_r2, 'wt') as r2_out:
            while True:
                # Each FASTQ entry consists of 4 lines
                header = f.readline().strip()
                seq = f.readline().strip()
                plus = f.readline().strip()
                qual = f.readline().strip()

                # Check for end of file
                if not header:
                    break

                if ".1 " in header:
                    # This is a R1 read
                    header = header.replace(".1 "," ")
                    r1_out.write(header + '\n')
                    r1_out.write(seq + '\n')
                    plus = plus.replace(".1 "," ")
                    r1_out.write(plus + '\n')
                    r1_out.write(qual + '\n')
                elif ".2 " in header:
                    # This is a R2 read
                    header = header.replace(".2 "," ")
                    r2_out.write(header + '\n')
                    r2_out.write(seq + '\n')
                    plus = plus.replace(".2 "," ")
                    r2_out.write(plus + '\n')
                    r2_out.write(qual + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Separate R1 and R2 reads from a single fastq.gz file.")
    parser.add_argument('input_file', type=str, help='Path to the input fastq.gz file.')
    parser.add_argument('output_r1', type=str, help='Path to the output R1 fastq.gz file.')
    parser.add_argument('output_r2', type=str, help='Path to the output R2 fastq.gz file.')

    args = parser.parse_args()

    separate_r1_r2(args.input_file, args.output_r1, args.output_r2)
