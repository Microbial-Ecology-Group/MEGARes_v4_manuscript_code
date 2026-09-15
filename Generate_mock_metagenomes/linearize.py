#!/usr/bin/env python3

import sys

def linearize_fasta(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        first_header = True
        
        for line in infile:
            line = line.strip()
            if not line:
                continue  # Skip blank lines
                
            if line.startswith('>'):
                # Add a newline before the next header (except the first one)
                if not first_header:
                    outfile.write('\n')
                outfile.write(line + '\n')
                first_header = False
            else:
                # Concatenate sequence lines without newlines
                outfile.write(line)
        
        # Write the final trailing newline at the end of the file
        if not first_header:
            outfile.write('\n')

if __name__ == '__main__':
    # Usage: python script.py input.fasta output.fasta
    linearize_fasta(sys.argv[1], sys.argv[2])

