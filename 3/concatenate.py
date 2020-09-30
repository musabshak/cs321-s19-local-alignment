def readtext(filename): 
    with open(filename, "r") as file:
        string = ''
        for line in file.readlines():
            array = line.strip().split()
            strng = ("".join(array[1::]))
            
            string+=strng

        return string
    
str1 = readtext('human.txt')
str2 = readtext('mouse.txt')

with open("humanmouse.txt", "w") as file:
    file.write(str1)
    file.write('\n')
    file.write(str2)

