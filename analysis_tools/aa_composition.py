from Bio import SeqIO
import glob
import time

def count_occurence(f, amino_acids_to_track='ACDEFGHIKLMNPQRSTVWY'):
    amino_acids_per_train = 0
    amino_acids = dict.fromkeys(amino_acids_to_track, 0)
    # or amino_acids = { k: 0 for k in amino_acids_to_track}

    for record in SeqIO.parse(f, "fasta"):
        for char in record.seq:
            if char in amino_acids:
                amino_acids_per_train += 1
                amino_acids[char] += 1

    percentages = { k : (v * 100.0) / amino_acids_per_train for k, v in amino_acids.items()}
    return percentages

def write_output(out_filename, from_file, percentages):
    with open(out_filename, 'a+') as out_file:
        for k, v in percentages.items():
            out_file.write('{} is {}% in {}\n'.format(k, v, from_file))

if __name__ == '__main__':
    current_time = time.strftime("%d.%m.%y %H:%M", time.localtime())
    out_file= 'test#{}.txt'.format(current_time)
    files = glob.glob('*.faa')
    for f in files:
        percentages = count_occurence(f)
        write_output(out_file, f, percentages)
