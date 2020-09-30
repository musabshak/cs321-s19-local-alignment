import numpy as np
from time import time
#from vtml20 import vtml20

#SCORE = vtml20
INDEL = -15
SCORE={('A', 'A'): 7, ('A', 'C'): -3, ('A', 'D'): -6, ('A', 'E'): -5, ('A', 'F'): -8, ('A', 'G'): -4, ('A', 'H'): -7, ('A', 'I'): -7, ('A', 'K'): -6, ('A', 'L'): -7, ('A', 'M'): -5, ('A', 'N'): -6, ('A', 'P'): -4, ('A', 'Q'): -5, ('A', 'R'): -7, ('A', 'S'): -2, ('A', 'T'): -3, ('A', 'V'): -3, ('A', 'W'): -9, ('A', 'Y'): -8, ('C', 'C'): 12, ('C', 'D'): -14, ('C', 'E'): -14, ('C', 'F'): -13, ('C', 'G'): -7, ('C', 'H'): -6, ('C', 'I'): -5, ('C', 'K'): -13, ('C', 'L'): -12, ('C', 'M'): -4, ('C', 'N'): -8, ('C', 'P'): -9, ('C', 'Q'): -13, ('C', 'R'): -7, ('C', 'S'): -3, ('C', 'T'): -5, ('C', 'V'): -3, ('C', 'W'): -15, ('C', 'Y'): -4, ('D', 'D'): 8, ('D', 'E'): -1, ('D', 'F'): -16, ('D', 'G'): -6, ('D', 'H'): -4, ('D', 'I'): -12, ('D', 'K'): -5, ('D', 'L'): -15, ('D', 'M'): -9, ('D', 'N'): -1, ('D', 'P'): -6, ('D', 'Q'): -4, ('D', 'R'): -12, ('D', 'S'): -5, ('D', 'T'): -6, ('D', 'V'): -9, ('D', 'W'): -10, ('D', 'Y'): -14, ('E', 'E'): 7, ('E', 'F'): -14, ('E', 'G'): -6, ('E', 'H'): -6, ('E', 'I'): -10, ('E', 'K'): -2, ('E', 'L'): -8, ('E', 'M'): -8, ('E', 'N'): -5, ('E', 'P'): -6, ('E', 'Q'): -1, ('E', 'R'): -10, ('E', 'S'): -5, ('E', 'T'): -6, ('E', 'V'): -7, ('E', 'W'): -16, ('E', 'Y'): -7, ('F', 'F'): 9, ('F', 'G'): -11, ('F', 'H'): -5, ('F', 'I'): -5, ('F', 'K'): -14, ('F', 'L'): -3, ('F', 'M'): -3, ('F', 'N'): -10, ('F', 'P'): -9, ('F', 'Q'): -8, ('F', 'R'): -10, ('F', 'S'): -7, ('F', 'T'): -8, ('F', 'V'): -6, ('F', 'W'): -3, ('F', 'Y'): 0, ('G', 'G'): 7, ('G', 'H'): -7, ('G', 'I'): -15, ('G', 'K'): -7, ('G', 'L'): -11, ('G', 'M'): -10, ('G', 'N'): -5, ('G', 'P'): -8, ('G', 'Q'): -8, ('G', 'R'): -7, ('G', 'S'): -4, ('G', 'T'): -8, ('G', 'V'): -10, ('G', 'W'): -9, ('G', 'Y'): -10, ('H', 'H'): 10, ('H', 'I'): -9, ('H', 'K'): -5, ('H', 'L'): -7, ('H', 'M'): -12, ('H', 'N'): -3, ('H', 'P'): -6, ('H', 'Q'): -2, ('H', 'R'): -3, ('H', 'S'): -5, ('H', 'T'): -5, ('H', 'V'): -8, ('H', 'W'): -6, ('H', 'Y'): -1, ('I', 'I'): 7, ('I', 'K'): -9, ('I', 'L'): -2, ('I', 'M'): -2, ('I', 'N'): -9, ('I', 'P'): -10, ('I', 'Q'): -9, ('I', 'R'): -8, ('I', 'S'): -9, ('I', 'T'): -5, ('I', 'V'): 1, ('I', 'W'): -6, ('I', 'Y'): -8, ('K', 'K'): 7, ('K', 'L'): -8, ('K', 'M'): -5, ('K', 'N'): -3, ('K', 'P'): -6, ('K', 'Q'): -2, ('K', 'R'): 0, ('K', 'S'): -5, ('K', 'T'): -4, ('K', 'V'): -8, ('K', 'W'): -9, ('K', 'Y'): -8, ('L', 'L'): 6, ('L', 'M'): 0, ('L', 'N'): -9, ('L', 'P'): -7, ('L', 'Q'): -6, ('L', 'R'): -8, ('L', 'S'): -8, ('L', 'T'): -7, ('L', 'V'): -3, ('L', 'W'): -6, ('L', 'Y'): -6, ('M', 'M'): 10, ('M', 'N'): -7, ('M', 'P'): -10, ('M', 'Q'): -4, ('M', 'R'): -6, ('M', 'S'): -8, ('M', 'T'): -4, ('M', 'V'): -3, ('M', 'W'): -13, ('M', 'Y'): -11, ('N', 'N'): 8, ('N', 'P'): -8, ('N', 'Q'): -4, ('N', 'R'): -5, ('N', 'S'): -2, ('N', 'T'): -4, ('N', 'V'): -9, ('N', 'W'): -10, ('N', 'Y'): -6, ('P', 'P'): 9, ('P', 'Q'): -5, ('P', 'R'): -7, ('P', 'S'): -4, ('P', 'T'): -6, ('P', 'V'): -7, ('P', 'W'): -9, ('P', 'Y'): -15, ('Q', 'Q'): 9, ('Q', 'R'): -2, ('Q', 'S'): -4, ('Q', 'T'): -5, ('Q', 'V'): -7, ('Q', 'W'): -15, ('Q', 'Y'): -12, ('R', 'R'): 8, ('R', 'S'): -6, ('R', 'T'): -6, ('R', 'V'): -9, ('R', 'W'): -8, ('R', 'Y'): -7, ('S', 'S'): 7, ('S', 'T'): -1, ('S', 'V'): -8, ('S', 'W'): -8, ('S', 'Y'): -6, ('T', 'T'): 8, ('T', 'V'): -4, ('T', 'W'): -15, ('T', 'Y'): -8, ('V', 'V'): 7, ('V', 'W'): -13, ('V', 'Y'): -8, ('W', 'W'): 12, ('W', 'Y'): -2, ('Y', 'Y'): 9}

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

    if score.get( (char1,char2), "NA" ) != "NA":
        return score[(char1,char2)]
    else:
        return score[(char2,char1)]
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
    
    
def GlobalScore(s1, s2, pointers=False):
    """
    Computes score matrix for global alignment in linear space.
    Setting optional argument "pointers" to True also computes
    a column of backtracking pointers.
    
    Args
        s1, s2: Two strings that are being aligned.
        pointers: (optional) Setting this to True ensures that backtracking
                pointers are stored for the last column of the score matrix.
    Returns
        Last column of the (n+1) x 2 score matrix. (In the case that the length of
        s1, the string on "top" of the matrix, is 0, the first column of the
        score matrix is returned).
        backTr: A (n+1) x 1 column matrix representing backtracking pointers.
                0, 1, 2 represent horizontal, diagonal, and vertical, respectively.
    """
    indel = INDEL
    n = len(s2)
    m = len(s1)
    
    score = np.full( (n+1,2), 0, dtype=np.int_)
    backTr = np.full( (n+1,1), -1, dtype=np.int_)
    backTr[0,0] = 0

    for i in range(1,n+1):
        score[i,0] = score[i-1,0] + indel
        
    for j in range(1,m+1):
        # On all except the first iteration, the score matrix will have to be updated;
        # the first column will need to take on the values of the second column
        if j>1: 
            score[:,0] = score[:,1]
        score[0,1] = score[0,0] + indel # Accounting for the top-right node in the score matrix
        
        for i in range(1,n+1):
            hor = score[i,0] + indel
            diag = score[i-1,0] + getscore(s2[i-1],s1[j-1])
            vert = score[i-1,1] + indel
            incoming = [hor, diag, vert]
            score[i,1] = max(incoming)
            
            if pointers == True and j==m: # j==m because only need pointers for last column
                backTr[i,0] = np.argmax(incoming)
        
    # In the case that the length of the string on "top" of the matrix is 0 the first column of the
    # score matrix is returned.
    if len(s1)==0:
        return score[:,0], backTr 
    else:
        return score[:,-1], backTr


def middle(s1, s2):
    """
    Finds the middle node and middle edge.
    
    Args
        s1, s2: Two strings to be aligned.
        
    Returns
        midCoord: A tuple representing the coordinates (indices) of the middle node in the alignment graph.
        middleEdge: A number (0,1,2) representing the edge coming out of the middle node.         
    """
    m = len(s1)
    n = len(s2)
    mid = m//2
    
    forward = GlobalScore( s1[:mid], s2)
    backward = GlobalScore( rev(s1[mid:]), rev(s2), pointers=True )
    
    fromSource = forward[0] # Get the last column of score matrix that GlobalScore returned
    toSink = backward[0]
    toSink = np.flip(toSink) # Reverse the toSink matrix
   
    length = fromSource + toSink
    rowMidNod = np.argmax(length)
    
    backTr  = np.flip(backward[1]).flatten() # Reverse, and flatten 2D backTr array to a 1D list   

    middleEdge = backTr[rowMidNod]
    midCoord = ( rowMidNod, mid  ) # Coordinates of middle node
    
    return midCoord, middleEdge


def LinearSpaceAlignment(s1, s2):
    """
    Implement linear space global alignment recursively.
    
    Args
        s1, s2: Two strings to be aligned.
    Returns
        Global alignment of the two strings.
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
    """
    Computes score matrix for local alignment in linear space.
    Setting optional argument "reverse" to True stops the function when
    the specified maxScore has been seen.
    
    Args
        s1, s2: Two strings that are being aligned.
        maxScore: (optional) The maximum score from running the LocalScore
                 in the forward direction.
        reverse: (optional) Is the algorithm being run on the reverse of the
                 two strings?
    Returns
        If reverse == False:
            maxAll: A tuple containing the maxScore and the coordinates in the
                    alignment matrix where the maxScore was reached. These coordinates
                    form the end coordinates of the local alignment.
        If reverse == True:
            startCoord: A tuple containing the start coordinates of the local
                        alignment, found by running LocalScore on the reverse
                        of the two strings and stopping it when maxScore has
                        been seen.
    """
    
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
    """
    Perform local alignment on two given sequnces.
    
    Args
        s1, s2: Two sequences to be aligned
        aligned: A tuple containing the two aligned sequences
    """
    maxAll = LocalScore(s1,s2)
    maxScore, endCoord = maxAll
    endrow, endcol = endCoord
    print(maxScore)
    
    # Run LocalScore on reverse of the strings and stop when the maxScore has been seen
    startCoord = LocalScore( rev(s1[:endcol]), rev(s2[:endrow]), maxScore=maxScore, reverse=True)
    
    # Need to translate startCoord with respect to the original alignment graph (the coordinates received
    # above are for the reverse alignment graph
    startCoord = (len(s2[:endrow]) - startCoord[0], len(s1[:endcol]) - startCoord[1] )
    
    startrow, startcol = startCoord
    
    # Run linear space global alignment between the start and end coordinates
    aligned = LinearSpaceAlignment(s1[startcol:endcol], s2[startrow:endrow])
    
    return aligned



if __name__ == "__main__":
    
    data = read("humanmouse1.txt")
    s1 = data[0]
    s2 = data[1]
 
    s = time()

    aligned1, aligned2 = LocalAlign(s1, s2)
    print(checkscore(aligned1, aligned2))
    print( (time()-s)/3600)
    
#    print(aligned1)
#    print(aligned2)
    with open("results1", "w") as file:
        file.write(aligned1)
        file.write('\n')
        file.write(aligned2)
    
    







