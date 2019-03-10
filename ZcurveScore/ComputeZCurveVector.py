from __future__ import division

import sys

BASES = 'ACGT'

def read_z_curve_score_params(paramFileName) :
    # Read numbers from the file, ignoring lines that start with #, 
    #   and return the list of numbers. There should be 190 of them.
    # The first 189 represent a vector and the last is an offset.
    import os
    params = []
    for line in open(os.path.abspath(os.path.expanduser(paramFileName))) :
        if line[0] == '#' :
            continue
        params += map(float, line.split())
    return params
 

def compute_z_curve_score_from_zvectors(zCurveVec, zCurveParams) :
    # Compute the score of the given vectors
    #print zCurveVec
    #print zCurveParams[:6]

    return sum(t * z for t, z in zip(zCurveParams[:6], zCurveVec)) + zCurveParams[6]


def compute_z_curve_score(seq, zCurveParams) :
    # Compute the score of the given sequence of As, Ts, Cs, and Gs (length 5 or more).
    zCurveVec = compute_z_curve_vector(seq)
    return sum(t * z for t, z in zip(zCurveParams[:189], zCurveVec)) + zCurveParams[189]

def compute_z_curve_vector(seq) :
    # Given a sequence of As, Ts, Cs, and Gs, return the 189-vector of phase specific 1, 2, and 3-mer
    #   frequencies described in Feng Gao and Chun-Ting Zhang, "Comparison of various algorithms for 
    #   recognizing short coding sequences of human genes" BIOINFORMATICS Vol. 20 no. 5 2004, 
    #   pages 673-681 DOI: 10.1093/bioinformatics/btg467
    # Sequence length must be at least 5.
    
    nn = len(seq)
    if nn < 5 :
        print >>sys.stderr, 'compute_z_curve_vector: sequence too short', len(seq)
        seq="ATGTGA"
    
    freq1mer = [{}, {}, {}] # freq1mer[0]['A'] is frequency of 'A' at 1st codon position, etc.
    freq2mer = [{}, {}, {}] # freq2mer[0]['AA'] is frequency of 'AA' starting at 1st codon position, etc.
    freq3mer = [{}, {}, {}] # freq3mer[0]['AAA'] is frequency of 'AAA' starting at 1st codon position, etc.
    
    for codonPos in [0, 1, 2] :
        for base1 in BASES :
            freq1mer[codonPos][base1] = 0
            for base2 in BASES :
                freq2mer[codonPos][base1 + base2] = 0
                for base3 in BASES :
                    freq3mer[codonPos][base1 + base2 + base3] = 0
    num1mers = [(nn + 2) // 3, (nn + 1) // 3, (nn) // 3 ]
    num2mers = [(nn + 1) // 3, (nn) // 3, (nn - 1) // 3]
    num3mers = [(nn) // 3, (nn - 1) // 3, (nn - 2) // 3]
    for ii in range(nn) :
        codonPos = ii % 3
        freq1mer[codonPos][seq[ii]] += 1 / num1mers[codonPos]
        if ii < nn - 1 :
            freq2mer[codonPos][seq[ii : ii + 2]] += 1 / num2mers[codonPos]
        if ii < nn - 2 :
            freq3mer[codonPos][seq[ii : ii + 3]] += 1 / num3mers[codonPos]
            
    xx = [] # [codonPos]
    yy = [] # [codonPos]
    zz = [] # [codonPos]
    xB = {} # [base][codonPos]
    yB = {} # [base][codonPos]
    zB = {} # [base][codonPos]
    xBB = {} # [base + base][codonPos]        
    yBB = {} # [base + base][codonPos]        
    zBB = {} # [base + base][codonPos]        
    for base1 in BASES :
        xB[base1] = []
        yB[base1] = []
        zB[base1] = []
        for base2 in BASES :
            xBB[base1 + base2] = []
            yBB[base1 + base2] = []
            zBB[base1 + base2] = []
    for codonPos in range(3) :
        f1 = freq1mer[codonPos]
        xx.append(f1['A'] + f1['G'] - f1['C'] - f1['T'])
        yy.append(f1['A'] + f1['C'] - f1['G'] - f1['T'])
        zz.append(f1['A'] + f1['T'] - f1['G'] - f1['C'])
        for base1 in BASES :
            f2 = freq2mer[codonPos]
            xB[base1].append(f2[base1 + 'A'] + f2[base1 + 'G'] - f2[base1 + 'C'] - f2[base1 + 'T'])
            yB[base1].append(f2[base1 + 'A'] + f2[base1 + 'C'] - f2[base1 + 'G'] - f2[base1 + 'T'])
            zB[base1].append(f2[base1 + 'A'] + f2[base1 + 'T'] - f2[base1 + 'G'] - f2[base1 + 'C'])
            for base2 in BASES :
                base12 = base1 + base2
                f3 = freq3mer[codonPos]
                xBB[base12].append(f3[base12 + 'A'] + f3[base12 + 'G'] - f3[base12 + 'C'] - f3[base12 + 'T'])
                yBB[base12].append(f3[base12 + 'A'] + f3[base12 + 'C'] - f3[base12 + 'G'] - f3[base12 + 'T'])
                zBB[base12].append(f3[base12 + 'A'] + f3[base12 + 'T'] - f3[base12 + 'G'] - f3[base12 + 'C'])
    result = []
    for codonPos in [0, 1, 2] :
        result += [xx[codonPos], yy[codonPos], zz[codonPos]]
    for base1 in BASES :
        for codonPos in [0, 1, 2] :
            result += [xB[base1][codonPos], yB[base1][codonPos], zB[base1][codonPos]]
    for base1 in BASES :
        for base2 in BASES :
            for codonPos in [0, 1, 2] :
                result += [xBB[base1 + base2][codonPos],
                           yBB[base1 + base2][codonPos],
                           zBB[base1 + base2][codonPos]]
    return result
