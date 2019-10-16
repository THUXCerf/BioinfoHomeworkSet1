"""
Global Alignment.py
Alignment

Created by Lukas Yu on 2019/10/15. All Rights Reserved
"""

import sys


# Define the function to compare the two string and return the scoring
def scoring(char1, char2):
    if char1 == char2:
        return score_alpha
    elif char1 == '-' or char2 == '-':
        return score_gamma
    else:
        return score_beta


# Initializing the scoring matrix
def intial(size):
    s_matrix = []
    for x in range(size[0]):
        s_matrix.append([])
        for y in range(size[1]):
            s_matrix[-1].append(0)
    return s_matrix


# Define the output function
def align(align1, align2):
    # Reverse the two alignment sequence
    align_alpha = align1[::-1]
    align_beta = align2[::-1]

    print("The alignment result is")
    print("Seq1:  ", align_alpha)
    print("Seq2:  ", align_beta)


# Define the function to fill scoring matrix, return the maximum globing alignment score and show the traceback result
def traceback(seq1, seq2):
    m, n = len(seq1), len(seq2)
    score = intial((m + 1, n + 1))

    # Fill progress
    for i in range(0, m + 1):                   # Initialing the column 0
        score[i][0] = score_gamma * i
    for j in range(0, n + 1):                   # Initialing the row 0
        score[0][j] = score_gamma * j
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1][j - 1] + scoring(seq1[i - 1], seq2[j - 1])
            insert_up = score[i - 1][j] + score_gamma
            insert_left = score[i][j - 1] + score_gamma
            score[i][j] = max(match, insert_up, insert_left)

    print("The global alignment score is %d." % (score[m][n]))

    # Traceback Progress
    align1, align2 = '', ''
    i, j = m, n

    # End traceback to the top or the left edge, according to the Needleman-Wunsch Algorithm
    while i > 0 and j > 0:
        score_current = score[i][j]
        score_diag = score[i-1][j-1]
        score_up = score[i][j-1]
        score_left = score[i-1][j]

        if score_current == score_diag + scoring(seq1[i-1], seq2[j-1]):
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1
        elif score_current == score_left + score_gamma:
            align1 += seq1[i-1]
            align2 += '-'
            i -= 1
        elif score_current == score_up + score_gamma:
            align1 += '-'
            align2 += seq2[j-1]
            j -= 1

    # Finish traceback up to the top left point (1,1)
    while i > 0:
        align1 += seq1[i - 1]
        align2 += '-'
        i -= 1
    while j > 0:
        align1 += '-'
        align2 += seq2[j - 1]
        j -= 1

    align(align1, align2)


# Define the alignment score of matching, gap penalty abd mismatch penalty
score_alpha = 1             # Matching score
score_beta = -1             # Mismatching penalty
score_gamma = -1            # Gap penalty

seq_alpha = input("Please input the first sequence: ")
seq_beta = input("Please input the second sequence: ")
traceback(seq_alpha, seq_beta)
