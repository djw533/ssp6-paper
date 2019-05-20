#! /usr/bin python

import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import translate
from Bio.Seq import reverse_complement





def parseArgs():
    """Parse command line arguments"""

    import argparse

    try:
        parser = argparse.ArgumentParser(
            description='Get proteome of gff file')

        parser.add_argument('-g',
                '--gff',
                action='store',
                help='Input gff file')
        parser.add_argument('-t',
                            '--type',
                            action = 'store',
	                    help="Type of output desired - 'nuc' or 'prot'")
        parser.add_argument('-o',
			    '--output',
			    action='store',
			    help='Output file - include .faa')

    except:
        print "An exception occurred with argument parsing. Check your provided options."
        traceback.print_exc()

    return parser.parse_args()




def _gff_to_fasta(gff_file, fasta_file, type):
    '''convert a gff file with the appended FASTA to
    gff_file = input gff file
    fasta_file = output file
    output: a protein coding FASTA file '''
    out = open("tmp_fasta.fa", "w")
    contigs = {}
    with open(gff_file) as f:
        fasta = False
        for line in f:
            if fasta:
                out.write(line)
                continue
            if line.startswith("##FASTA"):
                fasta = True
                continue
            if line.startswith("#"):
                continue
            toks = line.strip().split("\t")
            if toks[2] != "CDS":
                continue
            name = toks[-1].split("ID=")[1].split(';')[0]
            if toks[0] not in contigs:
                contigs[toks[0]] = []
            contigs[toks[0]].append({"name": name, "start": int(toks[3])-1, "stop": int(toks[4]), "strand": toks[6]})
    out.close()

    ## read the contigs and save the final fasta file
    out = open(fasta_file, "w")
    with open("tmp_fasta.fa") as handle:
        for values in SimpleFastaParser(handle):
            curr_contig = values[0].split()[0]
            print contigs.keys()
            if curr_contig not in contigs: # no CDSs in this contig
                continue
            print 'still continuing'
            for cds in contigs[curr_contig]:
                out.write(">" + cds["name"] + "\n")
                seq = values[1][cds["start"]:cds["stop"]]
                if cds["strand"] == "-":
                    seq = reverse_complement(seq)
                if type == 'prot':
                    out.write(translate(seq) + "\n")
                if type == 'nuc':
                    out.write(seq + "\n")
    out.close()
    os.remove("tmp_fasta.fa")
    return


def main():

    args = parseArgs()

    _gff_to_fasta(args.gff, args.output, args.type)



if __name__ == '__main__':
    main()
