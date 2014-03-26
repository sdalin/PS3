__author__ = 'Simona'

import random, sys, math

# -------------------------
# "Sequence" class that represents sequences, with
# some fields and methods helpful for using this class in a Gibbs sampler.
# -------------------------
class Sequence:
    # fields in this class
    # -------------------------
    seqName = ""  # name of this sequence (e.g. gene name)
    sequence = ""  # nucleotide sequence
    siteScores = []  # site odds ratios =P(motif | model)/P(motif | background)
    motif = -1  # current position of motif according to Gibbs sampler

    # methods in this class
    # -------------------------
    def __init__(self, name, seq):
        # initializes a new instance of a Sequence object, and initializes our
        # belief of where the motif is to be uniform across all possible motif
        # positions
        self.seqName = name
        self.sequence = seq

    def getMotif(self, *pos):
        # returns the motif of length motifWidth can either specify a position
        # (e.g. getMotif(4) returns the motif at position 4) or if no position
        # in specified it will return the motif at self.motif (e.g. getMotif()
        # returns self.sequence[self.motif : self.motif + motifWidth])
        if pos == ():
            idx = self.motif
        else:
            idx = pos[0]
        if idx < 0 or idx > len(self.sequence) - motifWidth:
            print "Error - tried to access motif beginning at", idx, \
                ", index out of bounds."
            sys.exit(1)
        else:
            return self.sequence[idx: idx + motifWidth]


# END OF SEQUENCE CLASS

def readFastaFile(filename):
    # INPUT: the name of a fasta file to read in OUTPUT: a list of Sequence
    # objects each corresponding to a sequence in the .fa file, initialized
    # according to the __init__() function under class Sequence
    try:
        fastaLines = open(filename).readlines()
    except IOError:
        # file "filename" cannot be opened, print an error message and exit
        print "Error - could not open", filename, + \
            ", please check that you entered the correct file."
        sys.exit(1)
    seqList = []
    curSeq = ''
    curName = None
    for line in fastaLines:
        if '>' in line:
            if curName != None:
                seqList.append(Sequence(curName, curSeq))
                curName = None
                curSeq = ''
            # remove first char '>' and leading/trailing whitespace
            curName = line.strip()[1:]
        else:
            curSeq = curSeq + line.strip().upper()
    if curName != None:
        seqList.append(Sequence(curName, curSeq))
    return seqList

def findSimpleBackgroundModel(sequences):
    # Finds background model assuming simple 0-order model (e.g. each position
    # independent) INPUTS:
    #	sequences - a list of Sequence objects
    # OUTPUTS:
    #	background - a dictionary mapping the four nucleotides to their
    #	frequencies across all sequences

    background = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    # -------------------------
    # PUT YOUR CODE HERE

    totalnts=0
    for sequenceObject in sequences:
        sequenceStr = sequenceObject.sequence
        for nt in sequenceStr:
            totalnts += 1
            if nt == 'A':
                background['A'] += 1
            elif nt == 'C':
                background['C'] += 1
            elif nt == 'G':
                background['G'] += 1
            elif nt == 'T':
                background['T'] += 1
            else:
                print "findSimpleBackground is broken!"

    for nt in background:
        background[nt] = float(background[nt])/float(totalnts)

    #print "findSimpleBackgroundModel() not yet implemented!"

    # -------------------------
    return(background)

def run():  # this is the main function

    # Get file name and motifWidth from command line
    if len(sys.argv) <= 1:
        print "Error-please specify a fasta file and motif width as inputs"
        sys.exit(1)
    else:
        fastaName = sys.argv[1]

    # Read in sequences (see readFastaFile)
    sequences = readFastaFile(fastaName)


    background = findSimpleBackgroundModel(sequences)

    print background

run()