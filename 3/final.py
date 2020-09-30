import numpy as np
from time import time
from Bio.SubsMat import MatrixInfo
from vtml20 import vtml20
from middle import middle


SCORE = vtml20
INDEL = -15


def read(filename):
    """
    Reads in input from file.
    
    Args
        filename: A string.        
    Returns
        (line1, line2): A tuple where line1 and line2 are strings to be
                        aligned.
    """
    with open (filename, "r") as file:
        line1 = file.readline().strip()
        line2 = file.readline().strip()
        
    return line1, line2


def rev(string):
    """
    Reverses a string.
    
    Args
        string: String to be reversed.
    Returns
        Reversed string.
    """
    return string[::-1]


def getscore(char1, char2):
    """
    Takes in two characters and returns corresonding match/mismatch
    score according to the given SCORE matrix. The characters could be
    amino acids or DNA nucleotides.
    
    Args:
        char1, char2: Two characters.
    Returns
        Match/mismatch score of char1 and char2 according to the scoring matrix.
    """
    score = SCORE
    char1 = char1.upper()
    char2 = char2.upper()
    

    return score[char1][char2]
#    if char1==char2: return 1
#    else: return 0


def checkscore(s1,s2):
    """
    Given two aligned strings, iterate through the length of the strings
    and compute the alignment score according to the getscore function.
    """
    score = 0
    indel = INDEL
    
    for i in range(len(s1)):
        if s1[i]=="-" or s2[i]=="-":
            score += indel
        else:
            score += getscore(s1[i], s2[i])
        
    return score
    

def LinearSpaceAlignment(s1, s2):
    """
    """
    if len(s1) == 0:
        return ("-" * len(s2), s2)
    if len(s2) == 0:
        return (s1, "-" * len(s1))
        
    midCoord, middleEdge = middle(s1, s2)
    row, col = midCoord

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


def LocalScore(s1, s2, maxScore=0, reverse=False):
    
    indel = INDEL
    n = len(s2)
    m = len(s1)
    
    score = np.full( (n+1,2), 0, dtype=np.int_)
    
    maxAll = (-9999, () )

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
                
                if reverse == True:
                    if maxAll[0] == maxScore:
                        startCoord = (i,j)
                        return startCoord

    return maxAll          

    
def LocalAlign(s1, s2):
    maxAll = LocalScore(s1,s2)
    maxScore, endCoord = maxAll
    endrow, endcol = endCoord
    print(maxScore)
    
    startCoord = LocalScore( rev(s1[:endcol]), rev(s2[:endrow]), maxScore=maxScore, reverse=True)
    startCoord = (len(s2[:endrow]) - startCoord[0], len(s1[:endcol]) - startCoord[1] )
    
    startrow, startcol = startCoord
    
    aligned = LinearSpaceAlignment(s1[startcol:endcol], s2[startrow:endrow])
    
    return aligned



if __name__ == "__main__":
    
    #data = read("humanmouse.txt")
    data = read("rosalind_ba5e.txt")
    s1 = data[0]
    s2 = data[1]
 
    s = time()

    aligned1, aligned2 = LocalAlign(s1, s2)
    print(aligned1)
    print(aligned2)
    
    print(checkscore(aligned1, aligned2))
    
    print(time()-s)







