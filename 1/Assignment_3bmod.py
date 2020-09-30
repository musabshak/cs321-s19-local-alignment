import numpy as np
from Bio.SubsMat import MatrixInfo
from time import *

def read(filename):
    """ Reads in input from file. Returns (line1, line2) where line1 and line2
        are strings to be aligned"""
    with open (filename, "r") as file:
        line1 = file.readline().strip()
        line2 = file.readline().strip()
        
    return line1, line2

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
    
def backAlign(s1, s2, bT):
    """ Given a backtracking matrix, aligns two strings appropriately
        and returns the aligned strings in a tuple """
    
    aligned1 = ""
    aligned2 = ""
        
    i = len(s2)
    j = len(s1)
    
    #print(i,j)
    #print(aligned1,aligned2)
    
    while ( not(i == 0 and j == 0) ):
        if bT[i,j] == 1:
            i -= 1 
            j -= 1
            aligned1 = s1[j] + aligned1
            aligned2 = s2[i] + aligned2
            
        elif bT[i,j] == 0:
            j -= 1
            aligned1 = s1[j] + aligned1
            aligned2 = "-" + aligned2
            
        elif bT[i,j] == 2:
            i -= 1
            aligned1 = "-" + aligned1
            aligned2 = s2[i] + aligned2
            
        #print(i,j)
        #print(aligned1,aligned2)
        
            
    return aligned1, aligned2        



def GlobalAlignment(s1, s2):
    """ Implements GlobalAlignment """
    
    indel = -5
    n = len(s2)
    m = len(s1)
    
    # matrix that keeps track of the score for each graph node
    score = np.full((n+1,m+1), 0, dtype=np.int_)
    
    # backtracking matrix
    bT = np.full((n+1,m+1), -1, dtype=np.int_) # 0 horizontal; 1 diagonal; 2 vertical

    for i in range(1,n+1):
        score[i,0] = score[i-1,0] + indel
        bT[i,0] = 2 # in the first column, all incoming edges are vertical
        
    for j in range(1,m+1):
        score[0,j] = score[0,j-1] + indel
        bT[0,j] = 0 # in the first row, all incoming edges are horizontal
        
    for i in range(1,n+1):
        for j in range(1,m+1):
            incoming = [score[i,j-1] + indel, score[i-1,j-1] + getscore(s2[i-1],s1[j-1]), score[i-1,j] + indel] # arrange according to numbers assigned to edges
            maxEdge = np.argmax(incoming)
            #maxEdge = maxind(incoming)
            bT[i,j] = maxEdge            
            score[i,j] = incoming[maxEdge]

    #print(score)
    #print(bT)
    return backAlign(s1,s2,bT), score[n,m]
    #return score

def checkscore(s1,s2):
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
    
    a = GlobalAlignment(s1,s2)
#    b = GlobalAlignment(s1[::-1],s2[::-1])
#    c = np.flipud(b)

#    print(a)
#    print(c)
#    print(a[:,3]+c[:,4])

    
    print(str(a[1]))
    print(a[0][0])
    print(a[0][1])
#    
    print(checkscore(a[0][0],a[0][1]))
    
        
        