import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('exons', metavar='exons', type=str, help='path to the file with exon coverages')
parser.add_argument('out', metavar='out', type=str, help='path to output file with gene coverages')
args = parser.parse_args()

gene="void"
gene_start=int
gene_end=int
mean_coverage=0
gene_length=0

with open(args.exons, 'r') as exons, open(args.out, 'w') as out:
    for line in exons:
        line = line.strip()
        cells = line.split('\t')
        exon = cells[0]
        coverage = float(cells[1])
        temp = str(cells[2])
        if (gene!=temp) and (gene!="void"):
            gene_coverage = round((mean_coverage/gene_length),2)
            row = chromosome+":"+str(gene_start)+"-"+str(gene_end)+"\t"+str(gene)+"\t"+str(gene_coverage)+"\t"+str(gene_length)+"\n"
            out.write(row)
            gene = cells[2]
            chromosome = exon.split(":")[0]
            start = int(exon.split(":")[1].split("-")[0])
            end = int(exon.split(":")[1].split("-")[1])
            gene_length = 0
            mean_coverage = 0
            gene_start = start
            gene_end = end
            length = end-start
            cov_exon = coverage*length
            mean_coverage = mean_coverage+float(cov_exon)
            gene_length = gene_length+length
        else:
            chromosome = exon.split(":")[0]
            gene = cells[2]
            start = int(exon.split(":")[1].split("-")[0])
            if gene_length==0:
                gene_start=start
            end = int(exon.split(":")[1].split("-")[1])
            length = end-start
            cov_exon = coverage*length
            mean_coverage = mean_coverage+float(cov_exon)
            gene_length = gene_length+length
            gene_end = end