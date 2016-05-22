__author__ = 'jlu96'

import csv
import sys


def get_genes(segment):
    if segment[-4:] in {'loss', 'gain'}:
        alt_type = segment[-4:]
        genes = segment[:-4].split('_')
    else:
        alt_type = 'somatic'
        genes = segment.split('_')

    return set(genes)

def get_parser():

    import argparse

    description = 'Find and replace every instance of one file with another'
    parser = argparse.ArgumentParser(description=description)

    # General parameters
    parser.add_argument('-c', '--COSMICFILE', default = "COSMICCensus.txt",
                        help='File name with conversions.')

    parser.add_argument('-ic', '--input_IDs', type=str, nargs='+', default='Name')
    parser.add_argument('-oc', '--output_IDs', default="", type=str, nargs='+', help="Must have same length as input ids")
    parser.add_argument('-i', '--annotate_file', required=True)
    parser.add_argument('-o', '--output_file', default=None)
    parser.add_argument('-ip', '--insertion_point', type=int, default=None, help="the index at which to insert the new output fieldnames")

    return parser


def run(args):

    # Load the parameters

    # Gene list
    COSMICFILE = args.COSMICFILE

    with open(COSMICFILE, 'rU') as f:
        cosmic_gene_list = set(l.rstrip().split()[0] for l in f if not l.startswith("#") and l)



    annotate_file = args.annotate_file

    output_file = args.output_file
    if not output_file:
        output_file = annotate_file + '_COSMIC.txt'

    genecols = args.input_IDs
    new_cols = args.output_IDs
    if not new_cols:
        new_cols = ["COSMICGenes_" + genecol for genecol in genecols]
    if len(new_cols) != len(genecols):
        print "Number of output IDs must be same as number of input IDs"
        raise ValueError

    records = []

    with open(annotate_file, 'rU') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        header = reader.fieldnames
        for row in reader:
            for i in range(len(genecols)):
                genecol = genecols[i]
                new_col = new_cols[i]
                segment = row[genecol]
                genes = get_genes(segment)
                cosmic_genes = genes.intersection(cosmic_gene_list)

                if cosmic_genes:
                    if segment[-4:] in {'loss', 'gain'}:
                        row[new_col] = '_'.join(list(cosmic_genes)) + segment[-4:]
                    else:
                        row[new_col] = '_'.join(list(cosmic_genes))
                else:
                    if segment[-4:] in {'loss', 'gain'}:
                        row[new_col] = 'NA' + segment[-4:]
                    else:
                        row[new_col] = 'NA'
            records.append(row)




    insertion_point = args.insertion_point
    new_header = header


    for new_col in reversed(new_cols):
        if insertion_point == None:
            new_header.append(new_col)
        else:
            new_header.insert(insertion_point, new_col)

    with open(output_file, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, delimiter='\t', fieldnames=new_header)
        writer.writeheader()
        for row in records:
            writer.writerow(row)

def main():
    run(get_parser().parse_args(sys.argv[1:]))


if __name__ == '__main__':
    main()