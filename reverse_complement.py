import sys

def complement(seq):
    nuc_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}
    new_seq = ""
    for nuc in seq:
        new_seq = new_seq + nuc_dict[nuc]
    return new_seq


def reverse(seq):
    return seq[::-1]

if __name__ == '__main__':
    sequence = sys.argv[1]
    print(reverse(complement(sequence)))
