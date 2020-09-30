cimport numpy as np
ctypedef np.int_t DTYPE_t



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

cimport cython
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
    
    cdef int n,m,indel
    indel = INDEL
    n = len(s2)
    m = len(s1)
    
    cdef np.ndarray[DTYPE_t, ndim=2] score = np.full( (n+1,2), 0, dtype=np.int_)
    
    cdef np.ndarray[DTYPE_t, ndim=1] backTr = np.full( (n+1,1), -1, dtype=np.int_)
    backTr[0,0] = 0

    cdef int i, j
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
        middleNode: A tuple representing the coordinates (indices) of the middle node in the alignment graph.
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