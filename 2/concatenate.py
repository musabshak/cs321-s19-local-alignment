def readtext(filename): 
    with open(filename, "r") as file:
        string = ''
        for line in file.readlines():
            array = line.strip().split()
            strng = ("".join(array[1::]))
            
            string+=strng

        return string

#def readtextdog(filename): 
#    with open(filename, "r") as file:
#        string = ''
#        for line in file.readlines():
#            array = line.strip().split()
#            strng = ("".join(array[::]))
#            
#            string+=strng
#
#        return string
    
str1 = readtext('human1.txt')
str2 = readtext('mouse1.txt')
str2=str2.lower()

#with open("humanmouse1.txt", "w") as file:
#    file.write(str1)
#    file.write('\n')
#    file.write(str2)

