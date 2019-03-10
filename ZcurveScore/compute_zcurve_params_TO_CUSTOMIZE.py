#!/usr/bin/env python
from __future__ import division
import sys
import rpy

###### INFO #######
##
##  1. Requires INPUT to be modified for file "INPUT_FOR_LDA" file
##  2. Requires NUM_VARS to be modified for the true number of VARS
##  So: RUN first --> cat $0 | sed "s/INPUT/$i/g" | sed "s/NUM_VARS/XX/g" > $0.tmp
##
##################

# Find a linear discriminant (affine actually) that distinguishes between 1, 2, and 3-mer 
# frequencies of coding and non-coding regions.

log_msg = print_time = sys.stderr.write
list_dotp = lambda vec1, vec2 : sum(v1 * v2 for v1, v2 in zip(vec1, vec2))

def compute_Z_curve_parameters(trainingDataFileName, outputThetaFileName, priorForPositive = 0.5, minSeqLength = 5) :
	#### INPUT FILE FORMAT #####
	#	The input file is a set of lines each of which consists of 195 tab-separated fields:
	# 	Field 0 is 1 for coding or 0 for non-coding.
	#	Fields 1 - 4 are chromosome, start, end, strand, which are ignored so you can leave them blank.
	#	Field 5 is the length of the original sequence of bases, which is there only so it can exclude sequences less than 5 bases (but of course you can do what you want with this -- it needn't be the actual number of bases if you don't want it to do this filtering).
	#	Fields 6-194 are the NUM_VARS z curve fields, computed by compute_z_curve_vector that I sent you earlier.
	##### HOW IT WORKS #######
    # Read the training data and invoke R function lda to calculate the theta vector for the linear discriminant. 
    # Then use the training data to compute an offset so that the Minimum Average Error (MAE)
    #     of the discriminant  f(vect) = dot(theta, vect) + offset   is at score threshold 0.
    # Print the minimum average error to stderr.
    # Write the parameters to outputThetaFileName out in the form of comment strings preceded by '#'
    #   then the NUM_VARS coordinates of the theta vector (separated by \n), then the offset 
    # Use priorForPositive as the prior for the linear discriminant, and as a weight when calculating the MAE.
    # Ignore training sequences shorter than minSeqLength nt (assumed to be the 5th field fo the training data).
    # Return (Fraction of negatives that are false at MAE point, Fraction of positives that are false at MAE point)
    log_msg('Using prior %f.' % priorForPositive)
    trainingVectors = []
    trainingClasses = []
    numEx = [0, 0]  # Num neg and pos training examples
    for count, line in enumerate(open(trainingDataFileName)) :
        words = line.split()
        length = int(words[4])
        if length < minSeqLength :
            continue
        cls = int(words[0])
        numEx[cls] += 1
        trainingClasses.append(cls)
        trainingVectors += [float(w) for w in words[6:]]
    log_msg('Using %d negative and %d positive examples.' % tuple(numEx))
    print_time('Creating matrix.')
    trainingMatrix = rpy.r.matrix(trainingVectors, ncol = NUM_VARS, byrow = True)
    rpy.r.require('MASS')
    print_time('Invoking lda.')    
    ldaResult = rpy.r.lda(trainingMatrix, trainingClasses, prior = [1 - priorForPositive, priorForPositive])
    theta = [vv[0] for vv in ldaResult['scaling']]  # ldaResult['scaling'] is a list of 1-element vectors
    print_time('Computing MAE.')
    scoreClassPairs = []
    for ii, cls in enumerate(trainingClasses) :
        trainingVect = trainingVectors[NUM_VARS * ii : NUM_VARS * (ii + 1)]
        score = list_dotp(trainingVect, theta)
        scoreClassPairs.append((score, cls))
    scoreClassPairs.sort()
    bestScore = None
    numSoFar = [0, 0]
    for score, cls in scoreClassPairs :
        numSoFar[cls] += 1
        fracFN = 1 - numSoFar[0] / numEx[0]
        fracFP = numSoFar[1] / numEx[1]
        aveErr = fracFN * (1 - priorForPositive) + fracFP * priorForPositive
        if bestScore == None or aveErr < bestAveErr :
            bestAveErr = aveErr
            bestScore = score
            fracFNAtBest = fracFN
            fracFPAtBest = fracFP
    log_msg('Minimum Average Error with prior %f is %f, at threshold %f.' % (priorForPositive, bestAveErr, bestScore))
    log_msg('Fraction of negatives that are false at the this threshold: %f' % fracFNAtBest)
    log_msg('Fraction of positives that are false at the this threshold: %f' % fracFPAtBest)    
    outFile = open(outputThetaFileName, 'w')
    print >>outFile, '# Z curve parameters for prior %g computed using training examples in' % priorForPositive
    print >>outFile, '#     ' + trainingDataFileName
    print >>outFile, '#     excluding examples shorter than %d nt.' % minSeqLength
    print >>outFile, '# Minimum Weighted Average Error is %f.' % bestAveErr
    print >>outFile, '# Fraction of negatives that are false at MAE point: %f' % fracFNAtBest
    print >>outFile, '# Fraction of positives that are false at MAE point: %f' % fracFPAtBest
    print >>outFile, '# First NUM_VARS lines are the coordinates of theta, last line is offset.'
    print >>outFile, '# Z curve score is    dot(theta, vect) + offset'
    for component in theta :
        print >>outFile, component 
    print >>outFile, -bestScore
    return fracFNAtBest, fracFPAtBest
    
    

compute_Z_curve_parameters("INPUT","INPUT.Params")
