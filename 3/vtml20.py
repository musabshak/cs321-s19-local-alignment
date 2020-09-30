import numpy as np
import pandas as pd

names = list("ACDEFGHIKLMNPQRSTVWY")

matrix = np.genfromtxt("vtml20.txt", dtype=np.intc, skip_header=1, max_rows=20, usecols=range(1,21))

vtml20 = pd.DataFrame(matrix, columns=names, index = names)