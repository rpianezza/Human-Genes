import argparse
import subprocess
import os

# Parse command-line arguments
parser = argparse.ArgumentParser(description="")
parser.add_argument("blasted", help="Path to fasta file with the sequences to BLAST to the library")
parser.add_argument("library", help="Path to fasta library file")
parser.add_argument("output", help="Path to the output file")
args = parser.parse_args()

'''

'''
def write_single_fasta(blasted, out, i):
    with open(blasted, 'r') as blasted:
        with open(out, 'w') as output_file:
            for counter, line in enumerate(blasted):
                if counter == (i*2)-2 or counter == (i*2)-1:
                    output_file.write(line)
                elif counter > (i*2)-1:
                    break

def count_rows(file):
    row_count = 0
    with open(file, 'r') as file:
        for line in file:
            row_count += 1
        file.seek(0)
    return row_count

tot_rows = count_rows(args.blasted)
rows = int(tot_rows/2)+1

def filter_blast(file, out):
    with open(file, 'r') as file:
        with open(out, 'w') as output_file:
            for line in file:
                cells = line.strip().split('\t')
                seq1 = cells[0]
                seq2 = cells[1]
                if seq1==seq2 or float(cells[12])<1000:
                    continue
                else:
                    cells[0]=seq1
                    cells[1]=seq2
                    new_line = '\t'.join(cells) + "\n"
                    output_file.write(new_line)

def calculate_similarity(file, out):
    with open(file, 'r') as file:
        with open(out, 'w') as output_file:
            geneA = "gene"
            geneB = "gene"
            matches = []
            for line in file:
                cells = line.strip().split('\t')
                pid = float(cells[2])
                match_length = int(cells[4])
                match_start = int(cells[7])
                match_end = int(cells[8])
                if (cells[1] != geneB) or (cells[0] != geneA):
                    if geneA!="gene":
                        homology = round((equal_bases/gene_length*100),2)
                        equal_bases = round(equal_bases)
                        gene_length = round(gene_length)
                        new_line = geneA + '\t' + geneB + '\t' + str(homology) + '\t' + str(equal_bases) + '\t' + str(gene_length) + '\n'
                        output_file.write(new_line)
                    geneB = cells[1]
                    matches = []
                    hit = (match_start, match_end)
                    matches.append(hit)
                    equal_bases = match_length*(pid/100)
                    if cells[0] != geneA:
                        print("Finding homology of " + str(geneA))
                        geneA = cells[0]
                        gene_length = int(cells[3])
                else:
                    # Keep first gene and second gene, add homology to previous
                    for segment in matches:
                        #print("Evaluating match (" + str(match_start) + ", " + str(match_end) + ")")
                        #print("On " + str(segment))
                        if (match_start >= segment[0]) and (match_end <= segment[1]):
                            #print("Nested or same, skipping")
                            #print("\n")
                            continue
                        elif ((match_start >= segment[0]) and (match_start < segment[1])) and (match_end > segment[1]):
                             match_start = segment[1]
                             #print("Moving match start to: " + str(match_start))
                        elif (match_start < segment[0]) and ((match_end < segment[1]) and (match_end > segment[0])):
                             match_end = segment[0]
                             #print("Moving match end to: " + str(match_end))
                        else:
                            continue
                            #print("Non overlapping")
                        #print("\n")
                    hit = (match_start, match_end)
                    matches.append(hit)
                    to_add = (match_end-match_start)*(pid/100)
                    equal_bases += to_add
                            

out_file=args.output
out_folder_path=out_file.strip().split('/')[:-1]
out_folder="/".join(out_folder_path)

merge_filtered_blast = f"cat {out_folder}/*filtered.tsv >> {out_folder}/merged.blast"
remove_single_blasts = f"rm {out_folder}/*filtered.tsv"

for row in range(1,rows):
    write_single_fasta(args.blasted, args.output+str(row)+".fasta", row)
    process = subprocess.Popen(["bash", "/Volumes/Temp1/human-genes/blast-matrix/blast-matrix-v2.0.sh", args.output+str(row)+".fasta", args.library, args.output+str(row)+".tsv"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    os.remove(args.output+str(row)+".fasta")
    filter_blast(args.output+str(row)+".tsv", args.output+str(row)+"filtered.tsv")
    os.remove(args.output+str(row)+".tsv")
    subprocess.run(merge_filtered_blast, shell=True)
    subprocess.run(remove_single_blasts, shell=True)

calculate_similarity(out_folder+"/merged.blast", out_folder+"/similarity.blast")

