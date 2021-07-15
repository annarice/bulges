import sys
import pandas as pd
import numpy as np
import re
import pickle
import random

MAX_ALLOWED_GAPS = 3
PA_SCORE_THRES = 14.75  # as 95% of the data
MATCH_SCORE = 1
MISMATCH_PENALTY = 0
GAP_PENALTY = -1.25

DNA_PAIRS_THERMODYNAMICS = {"AA": 9.1, "AT": 8.6, "TA": 6.0, "CA": 5.8, "GT": 6.5, "CT": 7.8, "GA": 5.6, "CG": 11.9,
                            "GC": 11.1, "GG": 11.0, "TT": 9.1, "TG": 5.8, "AC": 6.5, "AG": 7.8, "TC": 5.6,
                            "CC": 11.0}  # Breslauer et al.


def align_pair(seqA, seqB, match_score, mismatch_score, gap_score, gaps_allowed, non_gapped_5p_len=0):
    """
    returns the global alignment of seqA and seqB
    non_gapped_5p_len sets the length from left that cannot have gaps
    """

    # Initialization
    alignmentA = ""
    alignmentB = ""
    scoring_matrix = init_matrix(seqA, seqB, match_score, mismatch_score, gap_score, gaps_allowed, non_gapped_5p_len)
    j = len(seqA)
    i = len(seqB)
    c = gaps_allowed
    score = 0

    # scoring_matrix: lines represent seqB, columns represent seqA
    # meaning: decreasing i = proceeding with seqB
    #          holding j = putting gaps in AlignmentA

    while i > 0 or j > 0:
        if i > 0:
            charB = seqB[i - 1]
        if j > 0:
            charA = seqA[j - 1]

        if i > 0 and j > 0:
            sigma = two_chars_score(charA, charB, match_score, mismatch_score)
            diag_score = scoring_matrix[i - 1, j - 1, c] + sigma
            if c != 0:
                up_score = scoring_matrix[i - 1, j, c - 1] + gap_score
                left_score = scoring_matrix[i, j - 1, c - 1] + gap_score
            else:
                up_score = left_score = float("-inf")
            max_score = max(diag_score, up_score, left_score)

        else:  # have to initiate arguments
            diag_score = up_score = left_score = 0
            max_score = -1

        # check in which direction to head
        if diag_score == max_score:
            # diagonal - both sequences are aligned at position, no gap
            alignmentA = seqA[j - 1] + alignmentA
            alignmentB = seqB[i - 1] + alignmentB
            i -= 1
            j -= 1
            score += sigma

        elif j == 0 or up_score == max_score:
            # up - gap in seqA
            # base case: j==0 , seqA is completed (adding gaps to beginning), seqB not yet
            alignmentA = "-" + alignmentA
            alignmentB = seqB[i - 1] + alignmentB
            i -= 1
            c -= 1
            score += gap_score

        elif i == 0 or left_score == max_score:
            # left - gap in seqB
            # base case: i==0 , seqB id completed (adding gaps to beginning), seqA not yet
            alignmentA = seqA[j - 1] + alignmentA
            alignmentB = "-" + alignmentB
            j -= 1
            c -= 1
            score += gap_score

    return (alignmentA, alignmentB, score)


def init_matrix(seqA, seqB, match_score, mismatch_score, gap_score, gaps_allowed, non_gapped_5p_len):
    """
    initiates the matrix according to the global alignment function
    """

    scoring_matrix = np.full((len(seqB) + 1, len(seqA) + 1, gaps_allowed + 1), fill_value=float("-inf"), dtype=float)
    for c in range(0, gaps_allowed + 1):
        scoring_matrix[0, 0, c] = 0
        for l in range(max(non_gapped_5p_len, 1), c + 1):
            scoring_matrix[l, 0, c] = scoring_matrix[l - 1, 0, c] + gap_score
            scoring_matrix[0, l, c] = scoring_matrix[0, l - 1, c] + gap_score

    for i in range(1, len(seqB) + 1):
        for j in range(1, len(seqA) + 1):
            for c in range(0, gaps_allowed + 1):
                match = scoring_matrix[i - 1, j - 1, c] + two_chars_score(seqB[i - 1], seqA[j - 1], match_score,
                                                                          mismatch_score)
                if c != 0 and i >= non_gapped_5p_len and j >= non_gapped_5p_len:
                    delete = scoring_matrix[i - 1, j, c - 1] + gap_score
                    insert = scoring_matrix[i, j - 1, c - 1] + gap_score
                else:
                    delete = insert = float("-inf")

                scoring_matrix[i, j, c] = max(match, insert, delete)
    return scoring_matrix


def two_chars_score(C1, C2, match_score, mismatch_score):
    """
    returns the score of 2 chars comparison
    """
    if C1 == C2:
        return match_score
    else:
        return mismatch_score


def calculate_pa_score(seq1, seq2, match=MATCH_SCORE, mm=MISMATCH_PENALTY, gap=GAP_PENALTY):  # assumes PAM exists
    l = min(len(seq1), len(seq2))
    seq1 = seq1[-l:]
    seq2 = seq2[-l:]
    assert (len(seq1) == len(seq2))
    score = 0
    for i in range(l - 3):
        ch1 = seq1[i]
        ch2 = seq2[i]
        if ch1 == "-" or ch2 == "-":
            score += gap
        elif ch1 == ch2:
            score += match
        else:
            score += mm

    return score


def assert_reading_frame(seq, lst):
    # adjust reading frame to input length
    # sequences are between 20 and 25
    # 26 - 7 , 25 - 6, 24 - 5, 23 - 4, 22 - 3, 21 - 2, 20 - 1
    curr_length = 26
    while len(seq) < curr_length:
        lst.remove(max(lst))
        curr_length -= 1
    return lst


def best_pa(offtarget, sgRNA, match=MATCH_SCORE, mm=MISMATCH_PENALTY, gap=GAP_PENALTY):
    """
    :param offtarget: offtarget class object
    :param offtarget if iit's a string, then contains PAM and 3 nucleotides to the left
    :return:
    """
    extended_offtarget_seq = offtarget
    max_score = float("-inf")

    reading_frames_ind = [0, 6, 1, 5, 2, 4, 3]
    reading_frames_ind = assert_reading_frame(offtarget, reading_frames_ind)

    for i in reading_frames_ind:  # because starting at 3 is the original - we'd prefer that
        current_dna = extended_offtarget_seq[i:-3]
        (alnA, alnB, score) = \
            align_pair(seqA=sgRNA[:-3], seqB=current_dna, match_score=match, mismatch_score=mm, gap_score=gap,
                       gaps_allowed=MAX_ALLOWED_GAPS)  # regular pa
        # mm_start_score=0, mm_end_score=mm, gap_score=gap, gaps_allowed=MAX_ALLOWED_GAPS)
        if re.search("^\-", alnA) is None and score >= max_score:
            # the target can begin a '-', it means that the last nt of the sg is not paired
            # the sg cannot begin with a '-' - in that case, a better alignment would be found (shorter DNA).
            #   However, if we first found this target and then another with a different score- we'd prefer the other
            (alignmentA, alignmentB, max_score) = (alnA, alnB, score)
    aligned_sgRNA = alignmentA + sgRNA[-3:]
    aligned_offtarget = alignmentB + extended_offtarget_seq[-3:]

    return aligned_sgRNA, aligned_offtarget, max_score


def count_matches(aligned_sgRNA, aligned_offtarget, ignore_PAM=False):
    cnt = 0
    if ignore_PAM:
        ending = len(aligned_offtarget) - 3
    else:
        ending = len(aligned_offtarget)

    for i in range(ending):
        cnt += int(aligned_offtarget[i] == aligned_sgRNA[i])
    return cnt


def count_mismatches(aligned_sgRNA, aligned_offtarget, ignore_PAM=False):
    cnt = 0
    if ignore_PAM:
        ending = len(aligned_offtarget) - 3
    else:
        ending = len(aligned_offtarget)

    for i in range(ending):
        cnt += int(aligned_offtarget[i] != aligned_sgRNA[i] and aligned_offtarget[i] != "-" and aligned_sgRNA[i] != "-")
    return cnt


def count_bulges_in_seq(seq, ignore_PAM=False):
    if ignore_PAM:
        seq = seq[:-3]
    return seq.count("-")


def base_freq(seq, ignore_PAM=False):
    d = {"A": 0, "C": 0, "G": 0, "T": 0}
    if ignore_PAM:
        seq = seq[:-3]
    for base in d.keys():
        d[base] = seq.count(base) / len(seq) * 100
    return list(d.values())


def GC_content(seq, ignore_PAM=False):
    d = {"C": 0, "G": 0}
    if ignore_PAM:
        seq = seq[:-3]
    for base in d.keys():
        d[base] = seq.count(base)
    return sum(d.values()) / len(seq) * 100


def max_gap_length(seq, ignore_PAM=False):
    if ignore_PAM:
        seq = seq[:-3]
    # return max([len(match.group(0)) for match in re.finditer("-+", seq)]) # one-liner
    max_len = 0
    for match in re.finditer("-+", seq):
        if max_len < len(match.group(0)):
            max_len = len(match.group(0))
    return max_len


def distance_between_bulges(seq, ignore_PAM=False):
    if ignore_PAM:
        seq = seq[:-3]
    max_dist = 0
    for match in re.finditer("-+([ACGT]+)", seq):
        if max_dist < len(match.group(1)):
            max_dist = len(match.group(1))
    return max_dist


def bulge_location_in_seq(seq, ignore_PAM=False):
    locations = []
    if ignore_PAM:
        seq = seq[:-3]
    for match in re.finditer('-', seq):
        locations.append(match.start())
    locations += ['NA'] * (3 - len(locations))
    return locations


def bps_opposite_to_bulge(seq1, seq2, ignore_PAM=False):
    """
    seq1: the sequence with the bulges
    seq2: the sequence to extract the opposite bp from
    return: three opposite bulges bps; if there are less than three - return NA
    """
    if ignore_PAM:
        seq1 = seq1[:-3]
        seq2 = seq2[:-3]
    lst_of_bulge_locations = bulge_location_in_seq(seq1, ignore_PAM)
    opposite_bps = "000"
    for i in range(3):  # max number of bulges
        if lst_of_bulge_locations[i] is not "NA":
            opposite_bps = opposite_bps[:i] + seq2[lst_of_bulge_locations[i]] + opposite_bps[i + 1:]
    return opposite_bps
    # return [seq2[ind] for ind in lst_of_bulge_locations if ind is not "NA"] # one-liner version


def up_down_context(seq, ignore_PAM=False):
    if ignore_PAM:
        seq = seq[:-3]
    bulge_loc = bulge_location_in_seq(seq, ignore_PAM)
    for location in bulge_loc:
        if location is not "NA":
            if location == 0:  # in 5' end
                return [1]
            elif location == len(seq) - 1:  # in 3' end
                return [len(seq) - 1]
            else:
                return [location - 1, location + 1]  # in the middle


def mismatches_near_bulge(seq1, seq2, ignore_PAM=False):
    if ignore_PAM:
        seq1 = seq1[:-3]
        seq2 = seq2[:-3]
    context = up_down_context(seq1, ignore_PAM)

    for ind in context:
        if seq1[ind] != seq2[ind]:
            return 1
    return 0


def bulge_context(seq, target_context=None, ignore_PAM=False):
    if ignore_PAM:
        seq = seq[:-3]
    context = up_down_context(seq, ignore_PAM)

    try:
        if len(context) == 2:  # bulge in the middle of sequence
            return seq[context[0]], seq[context[1]]
        elif context[0] == 1:  # 5'
            if target_context is not None and isinstance(target_context, str):  # extracting for OT, retrieve from row[target_context]
                return target_context[target_context.find(seq) - 1], seq[1]
            else:
                return "NA", seq[1]  # extracting for sgRNA, retrieve only downstream
        else:  # 3'
            if target_context is not None and isinstance(target_context, str):  # extracting for OT, retrieve from row[target_context]
                return seq[len(seq) - 2], target_context[target_context.find(seq) + len(seq) + 1]
            else:
                return "NA", seq[1]  # extracting for sgRNA, retrieve only downstream
    except Exception as e:
        print(e)
        print("****" + seq)


def calculate_enthalpy(seq, ignore_PAM=False):
    if ignore_PAM:
        seq = seq[:-3]
    score = 0
    for i in range(1, len(seq)):
        current_couple = seq[i - 1:i + 1]
        score += DNA_PAIRS_THERMODYNAMICS[current_couple]
    return score


def calculate_dna_shape(seq, d, ignore_PAM=False):
    dna_shape_vars = ["MGW", "roll-1", "roll+1", "roll", "ProT", "HelT-1", "HelT+1", "HelT"]
    tmp = {"MGW": [], "roll-1": [], "roll+1": [], "roll": [], "ProT": [], "HelT-1": [], "HelT+1": [], "HelT": []}
    if ignore_PAM:
        seq = seq[:-3]
    m = len(seq)
    for i in range(2, m - 2):
        pentamer = seq[i - 2:i + 3]
        for var in dna_shape_vars:
            tmp[var].append(d[pentamer][var])
    results_dict = {key: round(np.mean(value), 2) for key, value in tmp.items()}
    return results_dict.values()


def complement(seq):
    nuc_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}
    new_seq = ""
    for nuc in seq:
        new_seq = new_seq + nuc_dict[nuc]
    return new_seq


def reverse(seq):
    return seq[::-1]


def pairwise_alignment_features(RNA, DNA, ignore_PAM=False):
    if ignore_PAM:
        DNA = DNA[:-3]
        RNA = RNA[:-3]
    no_PA = calculate_pa_score(DNA, RNA)
    (aligned_sgRNA, aligned_OT, PA_score) = best_pa(DNA, RNA)
    (aligned_sgRNA_rev, aligned_OT_rev, PA_score_rev) = best_pa(reverse(complement(DNA)), RNA)

    # check if reverse complement gives a better alignment
    if PA_score_rev > PA_score:
        aligned_sgRNA, aligned_OT, PA_score = aligned_sgRNA_rev, aligned_OT_rev, PA_score_rev
    matches = count_matches(aligned_sgRNA, aligned_OT)
    mismatches = count_mismatches(aligned_sgRNA, aligned_OT)

    return no_PA, aligned_sgRNA, aligned_OT, PA_score, matches, mismatches


def sequence_features(RNA, DNA, ignore_PAM=False):
    if ignore_PAM:
        DNA = DNA[:-3]
        RNA = RNA[:-3]
    base_freq_OT, base_freq_sgRNA = base_freq(DNA), base_freq(RNA)
    GC_OT, GC_sgRNA = GC_content(DNA), GC_content(RNA)

    return base_freq_OT, base_freq_sgRNA, GC_OT, GC_sgRNA


def bulges_features(aligned_seq1, aligned_seq2, bulge_count, target_context=None, ignore_PAM=False):
    if ignore_PAM:
        aligned_seq1 = aligned_seq1[:-3]
        aligned_seq2 = aligned_seq2[:-3]
    opposite_bps = bps_opposite_to_bulge(aligned_seq1, aligned_seq2)
    max_gap = max_gap_length(aligned_seq1)
    bulges_locations = bulge_location_in_seq(aligned_seq1)
    adjacent_mismatches = mismatches_near_bulge(aligned_seq1, aligned_seq2)
    up_context, down_context = bulge_context(aligned_seq1, target_context)
    bulges_distance = "NA"
    if bulge_count > 1:
        bulges_distance = distance_between_bulges(aligned_seq1)

    return opposite_bps, max_gap, bulges_locations, adjacent_mismatches, bulges_distance, up_context, down_context



def rigidity_features(RNA, DNA, ignore_PAM=False):
    if ignore_PAM:
        DNA = DNA[:-3]
        RNA = RNA[:-3]
    sgRNA_enthalpy, off_target_enthalpy = calculate_enthalpy(RNA), calculate_enthalpy(DNA)
    [sgRNA_shape_MGW, sgRNA_shape_roll_minus, sgRNA_shape_roll_plus, sgRNA_shape_roll, sgRNA_shape_prot,
     sgRNA_shape_helt_minus, sgRNA_shape_helt_plus, sgRNA_shape_helt] = calculate_dna_shape(RNA, DNAshape_pickle_dict)
    [OT_shape_MGW, OT_shape_roll_minus, OT_shape_roll_plus, OT_shape_roll, OT_shape_prot, OT_shape_helt_minus,
     OT_shape_helt_plus, OT_shape_helt] = calculate_dna_shape(DNA, DNAshape_pickle_dict)

    return sgRNA_enthalpy, off_target_enthalpy, \
           sgRNA_shape_MGW, sgRNA_shape_roll_minus, sgRNA_shape_roll_plus, sgRNA_shape_roll, sgRNA_shape_prot, sgRNA_shape_helt_minus, sgRNA_shape_helt_plus, sgRNA_shape_helt, \
           OT_shape_MGW, OT_shape_roll_minus, OT_shape_roll_plus, OT_shape_roll, OT_shape_prot, OT_shape_helt_minus, OT_shape_helt_plus, OT_shape_helt


def features_calculation(row):
    if DATASET == "sql":
        off_target = row["target_sequence"]
        sgRNA = row["grna_target_sequence"]
        target_context = row["target_context"]
    else:
        off_target = row["offtarget_sequence"]
        sgRNA = row["target"]
        target_context = None

    if off_target == "#NAME?" or sgRNA == "#NAME?":
        return

    # Ns are replaced with a random nucleotide
    #sgRNA = sgRNA.replace("N", random.choice(["A", "G", "C", "T"]))
    # Ns are replaced with the corresponding nucleotide
    sgRNA = sgRNA.replace("N", off_target[-3])
    #off_target = off_target.replace("N", random.choice(["A", "G", "C", "T"]))


    no_PA, aligned_sgRNA, aligned_OT, PA_score, matches, mismatches = pairwise_alignment_features(sgRNA, off_target)
    base_freq_OT, base_freq_sgRNA, GC_OT, GC_sgRNA = sequence_features(sgRNA, off_target)
    sgRNA_bulges, OT_bulges = count_bulges_in_seq(aligned_sgRNA), count_bulges_in_seq(aligned_OT)

    sgRNA_opposite_bps, OT_opposite_bps, sgRNA_max_gap, OT_max_gap, sgRNA_bulges_distance, OT_bulges_distances, sgRNA_bulges_locations, OT_bulges_locations, \
    sgRNA_adjacent_mismatches, OT_adjacent_mismatches, sgRNA_up_context, sgRNA_down_context, OT_up_context, OT_down_context, \
    sgRNA_enthalpy, off_target_enthalpy, \
    sgRNA_shape_MGW, sgRNA_shape_roll_minus, sgRNA_shape_roll_plus, sgRNA_shape_roll, sgRNA_shape_prot, sgRNA_shape_helt_minus, sgRNA_shape_helt_plus, sgRNA_shape_helt, \
    OT_shape_MGW, OT_shape_roll_minus, OT_shape_roll_plus, OT_shape_roll, OT_shape_prot, OT_shape_helt_minus, OT_shape_helt_plus, OT_shape_helt = \
        return_NAs(32)

    if sgRNA_bulges > 0:
        sgRNA_opposite_bps, sgRNA_max_gap, sgRNA_bulges_locations, sgRNA_adjacent_mismatches, sgRNA_bulges_distance, sgRNA_up_context, sgRNA_down_context = \
            bulges_features(aligned_sgRNA, aligned_OT, sgRNA_bulges)
    if OT_bulges > 0:
        OT_opposite_bps, OT_max_gap, OT_bulges_locations, OT_adjacent_mismatches, OT_bulges_distance, OT_up_context, OT_down_context = \
            bulges_features(aligned_OT, aligned_sgRNA, OT_bulges, target_context)

    sgRNA_enthalpy, off_target_enthalpy, \
    sgRNA_shape_MGW, sgRNA_shape_roll_minus, sgRNA_shape_roll_plus, sgRNA_shape_roll, sgRNA_shape_prot, sgRNA_shape_helt_minus, sgRNA_shape_helt_plus, sgRNA_shape_helt, \
    OT_shape_MGW, OT_shape_roll_minus, OT_shape_roll_plus, OT_shape_roll, OT_shape_prot, OT_shape_helt_minus, OT_shape_helt_plus, OT_shape_helt = \
        rigidity_features(sgRNA, off_target)

    return aligned_sgRNA, aligned_OT, no_PA, PA_score, matches, mismatches, \
           OT_bulges, sgRNA_bulges, sgRNA_opposite_bps, OT_opposite_bps, sgRNA_max_gap, OT_max_gap, \
           sgRNA_bulges_distance, OT_bulges_distances, sgRNA_bulges_locations, OT_bulges_locations, \
           sgRNA_adjacent_mismatches, OT_adjacent_mismatches, sgRNA_up_context, sgRNA_down_context, OT_up_context, OT_down_context, \
           sgRNA_enthalpy, off_target_enthalpy, \
           sgRNA_shape_MGW, sgRNA_shape_roll_minus, sgRNA_shape_roll_plus, sgRNA_shape_roll, sgRNA_shape_prot, sgRNA_shape_helt_minus, sgRNA_shape_helt_plus, sgRNA_shape_helt, \
           OT_shape_MGW, OT_shape_roll_minus, OT_shape_roll_plus, OT_shape_roll, OT_shape_prot, OT_shape_helt_minus, OT_shape_helt_plus, OT_shape_helt, \
           base_freq_sgRNA[0], base_freq_sgRNA[1], base_freq_sgRNA[2], base_freq_sgRNA[3], \
           base_freq_OT[0], base_freq_OT[1], base_freq_OT[2], base_freq_OT[3], GC_sgRNA, GC_OT


def return_NAs(n):
    return ("NA",) * n


if __name__ == '__main__':
    '''
    (aligned_sgRNA, aligned_OT, PA_score) = best_pa("TGGACGCATAAAGATGAGACGCTGG", "GACGCATAAAGATGAGACGCTGG")  # DNA, RNA
    no_PA, aligned_sgRNA, aligned_OT, PA_score, matches, mismatches = pairwise_alignment_features("CAGCAGGGCTGGGTGAGAAAG", "GAGCAGGGCTGGGGAGAAGG")
    print("Final alignment:")
    print(aligned_sgRNA)
    print(aligned_OT)
    print(matches, mismatches)
    exit()
    '''

    in_file = sys.argv[1]  # input file
    out_file = sys.argv[2]  # output file
    pickle_dir = sys.argv[3]  # data structure containing DNA-shape data located in ~/crispr-il/data/pickles/DNAshape.pkl
    DATASET = sys.argv[4]  # "sql" or "CHANGE"

    lst_of_columns = ["aligned_sgRNA", "aligned_OT", "score_without_PA", "PA_score", "total_matches",
                      "total_mismatches",
                      "RNA_bulges", "DNA_bulges", "RNA_opposite_bps", "OT_opposite_bps", "sgRNA_max_gap",
                      "OT_max_gap",
                      "sgRNA_bulges_distance", "OT_bulges_distances", "sgRNA_bulges_locations",
                      "OT_bulges_locations",
                      "sgRNA_adjacent_mismatches", "OT_adjacent_mismatches",
                      "sgRNA_up_context", "sgRNA_down_context", "OT_up_context", "OT_down_context",
                      "sgRNA_enthalpy", "off_target_enthalpy",
                      "sgRNA_shape_MGW", "sgRNA_shape_roll_minus", "sgRNA_shape_roll_plus", "sgRNA_shape_roll",
                      "sgRNA_shape_prot", "sgRNA_shape_helt_minus", "sgRNA_shape_helt_plus", "sgRNA_shape_helt",
                      "OT_shape_MGW", "OT_shape_roll_minus", "OT_shape_roll_plus", "OT_shape_roll", "OT_shape_prot",
                      "OT_shape_helt_minus", "OT_shape_helt_plus", "OT_shape_helt",
                      "A_freq_sgRNA", "C_freq_sgRNA", "G_freq_sgRNA", "T_freq_sgRNA",
                      "A_freq_OT", "C_freq_OT", "G_freq_OT", "T_freq_OT", "GC_sgRNA", "GC_OT"]

    data = pd.read_csv(in_file)

    if DATASET == "sql":
        data = data[~data.target_sequence.str.contains("-")]  # remove strings with "-"
        data = data[~data.grna_target_sequence.str.contains("-")]
        final_lst_of_columns = ["target_sequence", "grna_target_sequence", "cleavage_freq"] + lst_of_columns
    else:  # CHANGE-seq
        data = data[~data.target.str.contains("-")]  # remove strings with "-"
        data = data[~data.offtarget_sequence.str.contains("-")]
        final_lst_of_columns = data.columns.tolist() + lst_of_columns

    DNAshape_pickle_dict = pickle.load(open(pickle_dir, "rb"))
    # pickle load of dna_shape.
    # pickle created by the script /features_eng/process_features.py function create_dna_shape_dict.
    # the data is in /groups/itay_mayrose/annarice/crispr-il/data/pickles/DNAshape.pkl

    data[lst_of_columns] = data.apply(features_calculation, axis=1, result_type="expand")

    data.to_csv(out_file, index=False, columns=final_lst_of_columns)
