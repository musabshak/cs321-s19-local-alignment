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


def getscore(c1, c2, score):
    """ Take in two characters and return corresonding match/mismatch
        score according to the blosum62 scoring matrix """
    
    if score.get( (c1,c2), "NA" ) != "NA":
        return score[(c1,c2)]
    else:
        return score[(c2,c1)]

#    if c1==c2: return 1
#    else: return 0
    
    
def backAlign(s1, s2, bT, i, j):
    """ Given a backtracking matrix, aligns two strings appropriately
        and returns the aligned strings in a tuple """
    
    aligned1 = ""
    aligned2 = ""
       
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
            
        elif bT[i,j] == 3:
            i = 0
            j = 0
                        
        #print(i,j)
        #print(aligned1,aligned2)
        
    return aligned1, aligned2        



def LocalAlignment(s1, s2, score):
    """ Implements Local Alignment """
    
    n = len(s2)
    m = len(s1)
    indel=-5
    
    # matrix that contains the diagonal weights for matches/mismatches
    weightMat = np.full((n,m), -1, dtype=np.int_)
    
    # populate the weight matrix according to the blosum62 scoring matrix
    for i in range(len(weightMat)):
        for j in range(len(weightMat[0])):
            weightMat[i,j] = getscore(s2[i],s1[j], score)
    
    # matrix that keeps track of the score for each graph node
    S = np.full((n+1,m+1), 0, dtype=np.int_)
    
    # backtracking matrix
    bT = np.full((n+1,m+1), -1, dtype=np.int_) # 0 horizontal; 1 diagonal; 2 vertical

    for i in range(1,n+1):
        incoming = [-9999,-9999,S[i-1,0] + indel,0] # vertical
        maxEdge = np.argmax(incoming)
        S[i,0] = incoming[maxEdge]
        bT[i,0] = maxEdge
        
    for j in range(1,m+1):
        incoming = [S[0,j-1] +indel,-9999,-9999,0] # horizontal
        maxEdge = np.argmax(incoming)   
        S[0,j] = incoming[maxEdge]
        bT[0,j] = maxEdge
        
    for i in range(1,n+1):
        for j in range(1,m+1):
            # [horizontal, diagonal, vertical, 0]
            incoming = [S[i,j-1] +indel, S[i-1,j-1] + weightMat[i-1,j-1], S[i-1,j] +indel, 0] # taxi edge at index 3
            maxEdge = np.argmax(incoming)
            bT[i,j] = maxEdge            
            S[i,j] = incoming[maxEdge]
    
    ind = np.unravel_index(np.argmax(S, axis=None), S.shape) # start backtracking at highest score in S matrix
    i = ind[0]
    j = ind[1]
    
    print(S)
    
    return backAlign(s1,s2,bT,i,j), S[i,j]


    
if __name__ == "__main__":
    
    data = read("rosalind_ba5e.txt")
    s1 = data[0]
    s2 = data[1]
    
    
    score = MatrixInfo.pam250
    t = time()
    a = LocalAlignment(s1,s2, score)
    print("time taken: " + str((time()-t )/60) )
    print(a[1])
    print(a[0][0])
    print(a[0][1])
