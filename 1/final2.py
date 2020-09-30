import numpy as np
from Bio.SubsMat import MatrixInfo
from time import *
import final as final
import Assignment_3bmod as a3bmod

def read(filename):
    """ Reads in input from file. Returns (line1, line2) where line1 and line2
        are strings to be aligned """
    with open (filename, "r") as file:
        line1 = file.readline().strip()
        line2 = file.readline().strip()
        
    return line1, line2


def rev(s):
    return s[::-1]


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
          


def GlobalAlignment(s1, s2):
    """ Computes GlobalAlignment score matrix in linear space """
    
    indel = -5
    n = len(s2)
    m = len(s1)
    
    score = np.full( (n+1,2), 0, dtype=np.int_)
    
    maxAll = (-1, () )

    for i in range(1,n+1):
        vert = score[i-1,0] + indel
        taxi = 0
        incoming = [vert, taxi]
        score[i,0] = max(incoming)
        
    for j in range(1,m+1):
        if j>1: score[:,0] = score[:,1]
        hor = score[0,0] + indel
        taxi = 0
        incoming = [hor,taxi]
        #print(incoming)
        score[0,1] = max(incoming)
        for i in range(1,n+1):
            hor = score[i,0] + indel
            diag = score[i-1,0] + getscore(s2[i-1],s1[j-1])
            vert = score[i-1,1] + indel
            taxi = 0
            incoming = [hor, diag, vert, taxi]       
            score[i,1] = max(incoming)
            
            if max(incoming) >= maxAll[0]:
                maxAll = (max(incoming), (i,j) )
            # print(max(incoming))
    
        #print(score)

    #print(maxAll)
    return score, maxAll
                

def GlobalAlignment2(s1, s2, maxScore):
    
    indel = -5
    n = len(s2)
    m = len(s1)
    
    score = np.full( (n+1,2), 0, dtype=np.int_)
    
    maxAll = (-1, () )

    for i in range(1,n+1):
        vert = score[i-1,0] + indel
        taxi = 0
        incoming = [vert, taxi]
        score[i,0] = max(incoming)
        
    for j in range(1,m+1):
        if j>1: score[:,0] = score[:,1]
        hor = score[0,0] + indel
        taxi = 0
        incoming = [hor,taxi]
        #print(incoming)
        score[0,1] = max(incoming)
        #print(score)
        for i in range(1,n+1):
            hor = score[i,0] + indel
            diag = score[i-1,0] + getscore(s2[i-1],s1[j-1])
            vert = score[i-1,1] + indel
            taxi = 0
            incoming = [hor, diag, vert, taxi]       
            score[i,1] = max(incoming)
            
            if max(incoming) >= maxAll[0]:
                maxAll = (max(incoming), (i,j) )
                if maxAll[0] == maxScore:
                    startCoord = (i,j)
                    return startCoord

    

if __name__ == "__main__":
    import time
    data = read("rosalind_ba5e.txt")
    s1 = data[0]
    s2 = data[1]
    s = time.time()
    score, maxAll = GlobalAlignment(s1,s2)
    
    maxScore, endCoord = maxAll
    endrow,endcol = endCoord
    print(maxScore)
    
    startCoord = GlobalAlignment2(rev(s1[:endcol]), rev(s2[:endrow]), maxScore)
    
    
    startCoord = (len(s2[:endrow]) - startCoord[0], len(s1[:endcol]) - startCoord[1] )
    
    
    startrow, startcol = startCoord
    theend = final.LinearSpaceAlignment(s1[startcol:endcol], s2[startrow:endrow])
    #theend = a3bmod.GlobalAlignment(s1[startcol:endcol], s2[startrow:endrow])[0]
    print(theend[0])
    print(theend[1])
    
    print(final.checkscore(theend[0], theend[1]))
    print(time.time()-s)
