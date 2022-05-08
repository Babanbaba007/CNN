import csv
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from itertools import product
from collections import Counter


def csv2fasta(filename):
    # initializing the titles and rows list
    rows = []

    # reading csv file
    with open(filename, 'r') as csvfile:
        # creating a csv reader object
        csvreader = csv.reader(csvfile)
        fields = next(csvreader)

        # extracting each data row one by one
        for row in csvreader:
            rows.append(row)


    IDtaxon = []
    seqList = []

    for row in rows:
        temp = []
        temp.append(row[0])  # seq_id
        temp.extend(row[4:7])  # taxon from order to genus
        IDtaxon.append(temp)

        seq = row[-1] + '\n'
        seqList.append(seq)

    # create new records with fasta format
    new_record = ''
    for i in range(len(IDtaxon)):
        # replace taxon containing space with '_'   ex) E.col -> E_coli
        converter = lambda x: x.replace(' ', '_')
        new_taxon = list(map(converter, IDtaxon[i]))

        new_record += '>' + (' '.join(new_taxon)) + '\n' + (''.join(seqList[i]))

    with open(filename.split('.')[0] + '.fasta', "w") as f:  # rewrite file content
        f.write(new_record)
        f.close()



def fasta2matrix(filename,k):
    ALPHABET = "ACGT"
    kmers = [''.join(chars) for chars in product(*(k*(ALPHABET,)))]


    matrice = []
    for record in SeqIO.parse(filename, "fasta"):
        seq = str(record.seq)

        kmer_content = []
        for i in range(0, len(seq) -k+1):
            kmer_content.append(seq[i:i+k])

        counts = Counter(kmer_content)

        kmer_dict = {}
        for kmer in kmers:
            kmer_dict[kmer] = counts[kmer]

        temp = []
        temp.append(record.id)
        temp.append(kmer_dict)
        matrice.append(temp)



    with open(filename.split('.')[0] + '_matrix.fasta', "w") as f:  # rewrite file content
        first_row = 'seq_id' + ','
        for kmer in kmers:
            first_row += kmer + ','
        first_row = first_row[:-1] + '\n'
        f.write(first_row)

        for matrix in matrice:
            row = matrix[0] + ','
            for kmer in matrix[1]:
                row += str(matrix[1][kmer]) + ','
            row = row[:-1] + '\n'
            f.write(row)



def fasta2taxon(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        f.close()

    contents = ''
    for line in lines:
        if line[0] == '>':
            contents += line

    with open(filename.split('.')[0] + '_taxonomy.txt', 'w') as f:
        f.write(contents)
        f.close()


def matrix2input(matrixfile, taxonfile):
    filename = matrixfile.split('.')[0][:-6]
    taxonlevels = ['Order', 'Family', 'Genus']

    for taxonlevel in taxonlevels:
        outputFile = open(filename + "input_{}.txt".format(taxonlevel[0]), 'w')
        matrix = list(open(matrixfile, 'r'))
        records = open(taxonfile, 'r')

        outputFile.write(matrix[0])
        matrix = matrix[1:]

        i = 0
        for seq in records:
            elements = seq.split(" ")
            seqID_tax = elements[0][1:]
            taxon_index = taxonlevels.index(taxonlevel)
            taxon = elements[taxon_index+1]
            seqID_mat = matrix[i].split(",")[0]

            if seqID_mat == seqID_tax:
                outputFile.write(matrix[i].split('\n')[0] + ',' + taxon + '\n')
                i += 1
        outputFile.close()





#######################################################################################################################
if __name__ == "__main__":
    csvfile = sys.argv[1]
    k = int(sys.argv[2])

    csv2fasta(csvfile)
    fastafile = csvfile.split('.')[0] + '.fasta'

    fasta2matrix(fastafile,k)
    matrixfile = fastafile.split('.')[0] + '_matrix' + '.fasta'

    fasta2taxon(fastafile)
    taxonfile = fastafile.split('.')[0] + '_taxonomy' + '.txt'

    finalfile = matrix2input(matrixfile,taxonfile)

