# cmd: cd <DIR>; python PA.py

# from pairwise_alignment_analysis import PA_limitedIndel_unite_mm as PA_script
# from data_processing import arrange_data
# from utilities import fetchSequenceFromGenome
# from defs import *

def OT_bulge_context(seq, ignore_PAM=False):
    if ignore_PAM:
        seq = seq[:-3]
    indices = bulge_context(seq, ignore_PAM)
    if len(indices) == 2 and ((0 not in indices) or (len(seq) - 1 not in indices)):
        return seq[indices[0]], seq[indices[1]]


def PAM_type(seq):
    seq = seq[-2:]
    # NGG, NAG, NGA, NTG
    '''
    NAA
    NAC
    NAG
    NAT
    NCA
    NCC
    NCG
    NCT
    NGA
    NGC
    NGG
    NGT
    NTA
    NTC
    NTG
    NTT
    '''
    lst_of_PAMs = []
    second_position, third_position = "ACGT", "ACGT"
    for i in range(4):
        for j in range(4):
            lst_of_PAMs.append(second_position[i] + third_position[j])
    return lst_of_PAMs


def cnt_RNA_bulge(aligned_sgRNA, aligned_offtarget):
    return count_bulges_in_seq(aligned_offtarget)


def cnt_DNA_bulge(aligned_sgRNA, aligned_offtarget):
    return count_bulges_in_seq(aligned_sgRNA)


def count_bulges(seq1, seq2):
    cnt = 0
    ending = len(seq1) - 3
    for i in range(ending):
        cnt += seq1[i] == "-" or seq2[i] == "-"
    return cnt


def count_consecutive_inconsistencies(aligned_sgRNA, aligned_offtarget):
    cnt = 0
    current_cnt = 0

    for i in range(len(aligned_offtarget) - 3):
        if aligned_offtarget[i] != aligned_sgRNA[i]:
            current_cnt += 1
        else:
            cnt += current_cnt > 0
            current_cnt = 0
    return cnt


# nuc_upstream to bulges?                                                   X
# nuc_downstream to bulges?                                                 X
# what is in front of the bulge?                                            V V
# max length of gap                                                         V V
# if more than one bulge, what the distance between them?                   V V
# location of bulges in the protospacer                                     V V
# what is in the opposite strand? (maybe the strength of link is weak)      V V
# bps frequency, GC content (on unaligned)                                  V V
# PAM type                                                                  V
# shape (on unaligned)                                                      V
# enthalpy                                                                  V V
# are there mismatches near the bulges?                                     V V