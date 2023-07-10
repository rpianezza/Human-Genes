import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('scg', metavar='scg', type=str, help='path to the file with scg coverages')
parser.add_argument('pcg', metavar='pcg', type=str, help='path to the file with pcg coverages')
parser.add_argument('out', metavar='out', type=str, help='path to output file with normalized coverages')
args = parser.parse_args()

def scg_coverage(file_path):
    with open(file_path, 'r') as file:
        s = 0
        count = 0
        for line in file:
            line = line.strip()
            cells = line.split('\t')
            chromo = cells[0].split(':')
            if chromo[0]!="chrX":
                s += float(cells[1])
                count += 1
        scg_cov = s/count
        return scg_cov


def pcg_coverage(file_path, scg_cov, output):
    with open(file_path, 'r') as file, open(output, 'w') as out:
        out.write("position"+'\t'+"gene"+'\t'+"raw_coverage"+'\t'+"copynumber"+'\t'+"scg_mean"+'\n')
        for line in file:
            line = line.strip()
            cells = line.split('\t')
            old_cov = round((float(cells[2])),2)
            norm_cov = round((float(cells[2])/float(scg_cov)),2)
            new_row = cells[0]+'\t'+cells[1]+'\t'+str(old_cov)+'\t'+str(norm_cov)+'\t'+str(scg_cov)+'\n'
            out.write(new_row)

normalizer = scg_coverage(args.scg)
print("Mean SCG coverage (normalizer): "+str(round((normalizer),2)))
pcg_coverage(args.pcg, normalizer, args.out)