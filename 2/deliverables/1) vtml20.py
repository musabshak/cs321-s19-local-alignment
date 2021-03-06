import numpy as np
import itertools

names = "ACDEFGHIKLMNPQRSTVWY"

matrix = np.genfromtxt("vtml20.txt", dtype=np.intc, skip_header=1, max_rows=20, usecols=range(1,21))

rows, cols = np.triu_indices(20)

vtml20 = {}
    
for i in range(len(rows)):
    pair = (names[rows[i]], names[cols[i]])
    vtml20[pair] = matrix[rows[i], cols[i]]

## vtml20 = {('A', 'A'): 7, ('A', 'C'): -3, ('A', 'D'): -6, ('A', 'E'): -5, ('A', 'F'): -8, ('A', 'G'): -4, ('A', 'H'): -7, ('A', 'I'): -7, ('A', 'K'): -6, ('A', 'L'): -7, ('A', 'M'): -5, ('A', 'N'): -6, ('A', 'P'): -4, ('A', 'Q'): -5, ('A', 'R'): -7, ('A', 'S'): -2, ('A', 'T'): -3, ('A', 'V'): -3, ('A', 'W'): -9, ('A', 'Y'): -8, ('C', 'C'): 12, ('C', 'D'): -14, ('C', 'E'): -14, ('C', 'F'): -13, ('C', 'G'): -7, ('C', 'H'): -6, ('C', 'I'): -5, ('C', 'K'): -13, ('C', 'L'): -12, ('C', 'M'): -4, ('C', 'N'): -8, ('C', 'P'): -9, ('C', 'Q'): -13, ('C', 'R'): -7, ('C', 'S'): -3, ('C', 'T'): -5, ('C', 'V'): -3, ('C', 'W'): -15, ('C', 'Y'): -4, ('D', 'D'): 8, ('D', 'E'): -1, ('D', 'F'): -16, ('D', 'G'): -6, ('D', 'H'): -4, ('D', 'I'): -12, ('D', 'K'): -5, ('D', 'L'): -15, ('D', 'M'): -9, ('D', 'N'): -1, ('D', 'P'): -6, ('D', 'Q'): -4, ('D', 'R'): -12, ('D', 'S'): -5, ('D', 'T'): -6, ('D', 'V'): -9, ('D', 'W'): -10, ('D', 'Y'): -14, ('E', 'E'): 7, ('E', 'F'): -14, ('E', 'G'): -6, ('E', 'H'): -6, ('E', 'I'): -10, ('E', 'K'): -2, ('E', 'L'): -8, ('E', 'M'): -8, ('E', 'N'): -5, ('E', 'P'): -6, ('E', 'Q'): -1, ('E', 'R'): -10, ('E', 'S'): -5, ('E', 'T'): -6, ('E', 'V'): -7, ('E', 'W'): -16, ('E', 'Y'): -7, ('F', 'F'): 9, ('F', 'G'): -11, ('F', 'H'): -5, ('F', 'I'): -5, ('F', 'K'): -14, ('F', 'L'): -3, ('F', 'M'): -3, ('F', 'N'): -10, ('F', 'P'): -9, ('F', 'Q'): -8, ('F', 'R'): -10, ('F', 'S'): -7, ('F', 'T'): -8, ('F', 'V'): -6, ('F', 'W'): -3, ('F', 'Y'): 0, ('G', 'G'): 7, ('G', 'H'): -7, ('G', 'I'): -15, ('G', 'K'): -7, ('G', 'L'): -11, ('G', 'M'): -10, ('G', 'N'): -5, ('G', 'P'): -8, ('G', 'Q'): -8, ('G', 'R'): -7, ('G', 'S'): -4, ('G', 'T'): -8, ('G', 'V'): -10, ('G', 'W'): -9, ('G', 'Y'): -10, ('H', 'H'): 10, ('H', 'I'): -9, ('H', 'K'): -5, ('H', 'L'): -7, ('H', 'M'): -12, ('H', 'N'): -3, ('H', 'P'): -6, ('H', 'Q'): -2, ('H', 'R'): -3, ('H', 'S'): -5, ('H', 'T'): -5, ('H', 'V'): -8, ('H', 'W'): -6, ('H', 'Y'): -1, ('I', 'I'): 7, ('I', 'K'): -9, ('I', 'L'): -2, ('I', 'M'): -2, ('I', 'N'): -9, ('I', 'P'): -10, ('I', 'Q'): -9, ('I', 'R'): -8, ('I', 'S'): -9, ('I', 'T'): -5, ('I', 'V'): 1, ('I', 'W'): -6, ('I', 'Y'): -8, ('K', 'K'): 7, ('K', 'L'): -8, ('K', 'M'): -5, ('K', 'N'): -3, ('K', 'P'): -6, ('K', 'Q'): -2, ('K', 'R'): 0, ('K', 'S'): -5, ('K', 'T'): -4, ('K', 'V'): -8, ('K', 'W'): -9, ('K', 'Y'): -8, ('L', 'L'): 6, ('L', 'M'): 0, ('L', 'N'): -9, ('L', 'P'): -7, ('L', 'Q'): -6, ('L', 'R'): -8, ('L', 'S'): -8, ('L', 'T'): -7, ('L', 'V'): -3, ('L', 'W'): -6, ('L', 'Y'): -6, ('M', 'M'): 10, ('M', 'N'): -7, ('M', 'P'): -10, ('M', 'Q'): -4, ('M', 'R'): -6, ('M', 'S'): -8, ('M', 'T'): -4, ('M', 'V'): -3, ('M', 'W'): -13, ('M', 'Y'): -11, ('N', 'N'): 8, ('N', 'P'): -8, ('N', 'Q'): -4, ('N', 'R'): -5, ('N', 'S'): -2, ('N', 'T'): -4, ('N', 'V'): -9, ('N', 'W'): -10, ('N', 'Y'): -6, ('P', 'P'): 9, ('P', 'Q'): -5, ('P', 'R'): -7, ('P', 'S'): -4, ('P', 'T'): -6, ('P', 'V'): -7, ('P', 'W'): -9, ('P', 'Y'): -15, ('Q', 'Q'): 9, ('Q', 'R'): -2, ('Q', 'S'): -4, ('Q', 'T'): -5, ('Q', 'V'): -7, ('Q', 'W'): -15, ('Q', 'Y'): -12, ('R', 'R'): 8, ('R', 'S'): -6, ('R', 'T'): -6, ('R', 'V'): -9, ('R', 'W'): -8, ('R', 'Y'): -7, ('S', 'S'): 7, ('S', 'T'): -1, ('S', 'V'): -8, ('S', 'W'): -8, ('S', 'Y'): -6, ('T', 'T'): 8, ('T', 'V'): -4, ('T', 'W'): -15, ('T', 'Y'): -8, ('V', 'V'): 7, ('V', 'W'): -13, ('V', 'Y'): -8, ('W', 'W'): 12, ('W', 'Y'): -2, ('Y', 'Y'): 9}
## 


