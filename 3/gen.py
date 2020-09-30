from random import choice
"https://www.ncbi.nlm.nih.gov/protein/"
aminoacids = "ARNDCQEGHILKMFPSTWYV"

def f(x):
    return (10**-5)*((x)**2) + (0.0178)*(x) - 10.765




with open("proteinlike.txt","w") as file:
    
    for i in range(2500):
        file.write(choice(aminoacids))
        
    file.write('\n')
    
    for i in range(2500):
        file.write(choice(aminoacids))
        
        