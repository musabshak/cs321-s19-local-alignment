with open("alignedproteins1.txt", "r") as file:
    aligned = file.readline().strip().split('123')
    aligned1 = aligned[0]
    aligned2 = aligned[1]
    
#    with open("results1.txt", "w") as file1:
#        file1.write(aligned1)
#        file1.write('\n')
#        file1.write(aligned2)
 
INDEL = 0
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
#    score = SCORE
#    char1 = char1.upper()
#    char2 = char2.upper()
#
#    if score.get( (char1,char2), "NA" ) != "NA":
#        return score[(char1,char2)]
#    else:
#        return score[(char2,char1)]
    if char1==char2: return 1
    else: return 0


def checkscore(s1,s2):
    """
    Given two aligned strings, iterate through the length of the strings
    and compute the alignment score according to the getscore function.
    """
    score = 0
    indelnum = 0
    indel = INDEL
    
    for i in range(len(s1)):
        if s1[i]=="-" or s2[i]=="-":
            score += indel
            indelnum+=1
        else:
            score += getscore(s1[i], s2[i])
        
    return score, indelnum