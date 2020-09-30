import numpy as np
from Bio.SubsMat import MatrixInfo
from time import *

def read(filename):
    """ Reads in input from file. Returns (line1, line2) where line1 and line2
        are strings to be aligned """
    with open (filename, "r") as file:
        line1 = file.readline().strip()
        line2 = file.readline().strip()
        
    return line1, line2


def rev(s):
    """ Returned a reversed string """
    return s[::-1]

def maxind(s):
    """ Find the index of the last position where the maximum
        of a list occurs """
    
    maxim = -99999
    ind = -1    
    for i in range(len(s)):
        if s[i] >= maxim:
            maxim = s[i]
            ind = i   
    return ind


def getscore(c1, c2):
    """ Take in two characters and return corresonding match/mismatch
        score according to the blosum62 scoring matrix """
    score = MatrixInfo.pam250
    
    if score.get( (c1,c2), "NA" ) != "NA":
        return score[(c1,c2)]
    else:
        return score[(c2,c1)]
#    if c1==c2: return 1
#    else: return 0
    
    
def GlobalScore(s1, s2, pointers=False):
    """ Computes GlobalAlignment score matrix in linear space.
        Setting optional argument "pointers" to True also computers
        a column on backtracking pointers.
        
        We are pretty confident that this function works fine.
    """
    
    indel = -5
    n = len(s2)
    m = len(s1)
    
    score = np.full( (n+1,2), 0, dtype=np.int_)
    bT = np.full( (n+1,1), -1, dtype=np.int_)
    bT[0,0] = 0

    for i in range(1,n+1):
        score[i,0] = score[i-1,0] + indel
        
    for j in range(1,m+1):
        
        if j>1:
            score[:,0] = score[:,1]
        score[0,1] = score[0,0] + indel
        
        for i in range(1,n+1):
            hor = score[i,0] + indel
            diag = score[i-1,0] + getscore(s2[i-1],s1[j-1])
            vert = score[i-1,1] + indel
            incoming = [hor, diag, vert]       
            score[i,1] = max(incoming)
            
            if pointers == True and j==m: 
                bT[i,0] = np.argmax(incoming)
                
    #print(score)
    if len(s1)==0:
        return score[:,0], bT # return last column of score matrix and backtracking column matrix
    else:
        return score[:,-1], bT

def middle(s1, s2):
    """
    Returns
        middleNode: a tuple representing the coordinates of the middle node in a graph
        middleEdge: a number (0,1,2) representing the edge coming out of the middlenode
        
    """
    
    m = len(s1)
    n = len(s2)
    mid = m//2
    
    forward = GlobalScore( s1[:mid], s2)
    backward = GlobalScore( rev(s1[mid:]), rev(s2), pointers=True )
    
    fromSource = forward[0]
    toSink = backward[0]
    
    #print(fromSource)
    #print(toSink)
    toSink = np.flip(toSink)
    #print(toSink)
    
    length = fromSource + toSink
#    print(fromSource)
#    print(toSink)
#    print(length)
    
    bT  = np.flip(backward[1]).flatten() # flatten 2D array to a list   
    
    #print(length)
    #rowMidNod = maxind(length)
    #print(length)
    rowMidNod = np.argmax(length)
    #print(rowMidNod)
    #print(np.argmax(length))
    
    middleEdge = bT[rowMidNod]
    midCoord = ( rowMidNod, mid  ) # coordinates of middle node
    
    return midCoord, middleEdge

"""
GAGCAATTGCGG
ACTTAACTG
"""

def LinearSpaceAlignment(s1, s2):
       
    if len(s1) == 0:
        return ("-" * len(s2), s2)
    if len(s2) == 0:
        return (s1, "-" * len(s1))
        
    
    midCoord, middleEdge = middle(s1, s2)
    row, col = midCoord
#    print(s1,s2)
#    print(midCoord)
#    print(middleEdge)
    
    aligned1, aligned2 = LinearSpaceAlignment( s1[:col], s2[:row] )
    
    if middleEdge == 0:
        aligned1 += s1[col]
        aligned2 += '-'
        suffix1, suffix2 = LinearSpaceAlignment( s1[col+1:], s2[row:] )
               
    elif middleEdge == 1:
        aligned1 += s1[col]
        aligned2 += s2[row]
        suffix1, suffix2 = LinearSpaceAlignment( s1[col+1:], s2[row+1:] )
        
    elif middleEdge == 2:
        aligned1 += '-'
        aligned2 += s2[row]
        suffix1, suffix2 = LinearSpaceAlignment( s1[col:], s2[row+1:] )
  

    return (aligned1+suffix1, aligned2+suffix2)


def checkscore(s1,s2):
    """
    Given two aligned strings, iterate through the length of the strings
    and computer score according to getscore function.
    """
    score = 0
    indel = -5
    
    for i in range(len(s1)):
        if s1[i]=="-" or s2[i] =="-":
            score += indel
        else:
            score += getscore(s1[i], s2[i])
        
    return score

if __name__ == "__main__":
    
    data = read("rosalind_ba5e.txt")
    s1 = data[0]
    s2 = data[1]
 
    
#    GlobalScore(s1,s2)
#    midCoord, midEdge = middle(s1,s2)
#    print(midCoord, midEdge)
    
    aligned1, aligned2 = LinearSpaceAlignment(s1,s2)
    print(aligned1)
    print(aligned2)
    print(checkscore(aligned1,aligned2))