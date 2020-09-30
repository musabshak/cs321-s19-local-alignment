import numpy as np
import itertools

names = "ACDEFGHIKLMNPQRSTVWY"

matrix = np.genfromtxt("vtml20.txt", dtype=np.intc, skip_header=1, max_rows=20, usecols=range(1,21))

rows, cols = np.triu_indices(20)

vtml20 = {}
    
for i in range(len(rows)):
    pair = (names[rows[i]], names[cols[i]])
    vtml20[pair] = matrix[rows[i], cols[i]]




