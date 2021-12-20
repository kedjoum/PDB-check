import datetime as date
import gzip
import os
import shutil
from Bio import SeqIO
import wget as wg
import requests

workDirectory = os.getcwd()


# we write a fonction to dowload the file if the date is Wednessday

def downloadfile(url, filename):
    actual_execute_date = date.datetime.today().weekday()
    if actual_execute_date == 2:
        if filename == 'current_author.idx':
            filepath = workDirectory + '/last_author.idx'
            if os.path.exists(filepath):
                os.remove(filepath)  # we remove the last file named last_author.idx if it  exist
                os.renames("current_author.idx", "last_author.idx")  # current file become the last file
            wg.download(url, workDirectory + '/' + filename)  # we download the new current author file

        elif filename == 'entry_type.txt' or filename == 'pdb_seqres.txt.gz':
            os.remove(filename)  # remove last entry_type or last PDB_seqres  and download the new one
            wg.download(url, workDirectory + '/' + filename)


# we write a function to decompress a file
def unziprnafile(zipfilename, unzipfilename):
    with gzip.open(zipfilename, 'rb') as f_in:
        with open(unzipfilename, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


# We write a founction to get ID
def getidcode(filename):
    with open(filename, "r") as datei:  # open, read and close a file
        content = datei.readlines()  # read file  line to line and save them into a list
        idcodes = []
        for element in content[5:]:
            idcode, author = element.split(" ;")
            idcodes.append(idcode.lower())  # get Id of the new file
    return idcodes


# We write function to compare RNA-sequence with Rfam
def rfamsearch(fileName):
    url = "https://rfam.org/search/batch"
    payload = {'email': 'phangou10@gmail.com'}  # put the email in which you would like to receive the results
    files = [
        ('batchSeq', (fileName, open(fileName, 'rb'), 'application/octet-stream'))
    ]
    # move caracter like "'" , "*" and X in RNA file because Rfam does not recognise them.
    sequences = []
    with open(fileName) as f:
        for sequence in f:
            sequence = sequence.replace("*", '')
            if sequence.startswith(">") and "(5'" in sequence:
                sequence.replace("'", '')
            elif "X" in sequence:
                sequence = sequence.replace("X", '')
            sequences.append(sequence)
    with open(fileName, "w") as datei2:
        datei2.writelines(sequences)
    # send the resquest to Rfam and verifie if you have errors
    resp = requests.request("POST", url, data=payload, files=files)
    if resp.ok:
        return resp
    else:
        raise TypeError("check if the fasta file contains special characters")


# We dowload the PDB-file
author_idurl = 'ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/author.idx'
downloadfile(author_idurl, 'current_author.idx')

# We extract the IDcode from the file

current_idcode = getidcode('current_author.idx')
last_idcode = getidcode('last_author.idx')

# We extract new IDcode
new_idcode = list(set(current_idcode) - set(last_idcode))

# We download file that content IDcode associated with Nucleid acid or Protein
id_contenturl = 'ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt'
downloadfile(id_contenturl, 'entry_type.txt')
# We extract New IDcode  containing nucleic acids
id_type = {}

with open('entry_type.txt', 'r') as entry_type:
    content_entry = entry_type.readlines()
    for line in content_entry:
        IDCODE_CONTENT = line.split()[0:2]
        if IDCODE_CONTENT[1] != 'prot':
            id_type[IDCODE_CONTENT[0]] = IDCODE_CONTENT[1]

new_idcode_content = {key: id_type[key] for key in id_type.keys() & set(new_idcode)}
print(new_idcode_content.keys())
# We dowload and unzip file with protein- and Nucleid acid sequences

pdb_sequrl = 'ftp://ftp.pdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz'

downloadfile(pdb_sequrl, 'pdb_seqres.txt.gz')

unziprnafile('pdb_seqres.txt.gz', 'pdb_seqres.txt')

# using new_Idcode_content we sort new idcode that content just RNA sequences
# by filtering just RNA-sequences that content U

new_rna = 'new_rna.fa'
new_rna_seq = []

for seq_record in SeqIO.parse('pdb_seqres.txt', 'fasta'):
    header = seq_record.name
    description = seq_record.description
    seq = seq_record.seq
    for accession_id in new_idcode_content.keys():
        if accession_id == header[0:4] and 'U' in seq:  # because we just want RNA
            description = ">" + str(description) + "\n"
            seq = str(seq) + "\n"
            lines = [description, seq]
            new_rna_seq.append(lines)

file = open("new_rna.fa", "w")
for elt in new_rna_seq:
    file.writelines(elt)
file.close()

# Now comparing the new RNA with the RNA in Rfam using curl

response = rfamsearch("new_rna.txt")
print(response.status_code)
print(response.text)
