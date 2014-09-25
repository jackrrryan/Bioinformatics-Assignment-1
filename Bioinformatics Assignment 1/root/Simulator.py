'''
Created on Sep 19, 2014

@author: Jack Ryan B00321033
'''

from random import randint

file = open('Assnt1_sampleinput.fna.txt', 'r')  #opening the FASTA file
template = file.read()  #reading the FASTA file to a string

#some variables used in all calculations
nucleotides = ['A','T','C','G']
antinucleotides = ['T','A','G','C']
preamble = template[0:72]   #This is all the text before the start of the sequence
template = template[73:]    #This is the actual sequence
templateList = []   #this is a copy of the template string to be used in determining gaps in sequencing
def makeTemplateList(template):
    for i in range(0,len(template)):
        templateList.append(template[i])

makeTemplateList(template)

def readSizePosition(sequenceType):
    if sequenceType == 1:
        readLength = 700  #700bp is length of pyrosequencing reads
        start = randint(0,len(template)-readLength) #read must end before last bp in sequence
        finish = start+readLength
    elif sequenceType == 2:
        readLength = randint(50,300)    #Illumina can have reads between 50-300 bp
        start = randint(0,(len(template)-readLength))
        finish = start+readLength
    else:
        readLength = randint(400,900)   #Sanger can have reads between 400-900 bp
        start = randint(0,(len(template)-readLength))
        finish = start+readLength
    return (start,finish)

def sequence(sequenceType, readNum, template, preamble, antinucleotides, nucleotides):
    
    readsfile = open('reads.txt', 'w')
    readsfile.write(preamble + '\n\nBegin sequencing:\n\n')

    print("# of reads: {}".format(readNum))

    counter = 0

    for i in range(0,readNum):
        antisense = ''
        sfpositions = readSizePosition(int(sequenceType))
        start = sfpositions[0]    # this is the start position of this specific read
        finish = sfpositions[1]   # this is the end position of this specific read
        templateList = findGaps(start, finish)
        for j in range(start,finish):
            templateN = template[j]
            if templateN == 'A':
                antisense = antisense + antinucleotides[0]
            elif templateN == 'T':
                antisense = antisense + antinucleotides[1]
            elif templateN == 'C':
                antisense = antisense + antinucleotides[2]
            else:
                antisense = antisense + antinucleotides[3]
        readsfile.write(">Read{}[{}:{}]\n{}\n\n".format(counter,start,finish,antisense))
        counter += 1

    readsfile.flush()
    readsfile.close()
    return templateList

def randsequence(sequenceType, readNum, template, preamble, antinucleotides, nucleotides):
    if sequenceType == 1:
        errorpercent = 1    #1/1000, or 99.9% for Pyrosequencing
    elif sequenceType == 2:
        errorpercent = 20   #20/1000, or 98% for Illumina
    else:
        errorpercent = 1    #1/1000, or 99.9% for Sanger

    readsfile = open('reads.txt', 'w')
    readsfile.write(preamble + '\n\nBegin sequencing:\n\n')

    counter = 0

    for i in range(0,readNum):
        antisense = ''
        sfpositions = readSizePosition(int(sequenceType))
        start = sfpositions[0]    # this is the start position of this specific read
        finish = sfpositions[1]   # this is the end position of this specific read
        templateList = findGaps(start, finish)
        for j in range(start,finish):
            templateN = template[j]
            if templateN == 'A':
                if randint(0,1000) <= errorpercent:
                    antisense = antisense + antinucleotides[randint(0,3)]
                else:
                    antisense = antisense + antinucleotides[0]
            elif templateN == 'T':
                if randint(0,1000) <= errorpercent:
                    antisense = antisense + antinucleotides[randint(0,3)]
                else:
                    antisense = antisense + antinucleotides[1]
            elif templateN == 'C':
                if randint(0,1000) <= errorpercent:
                    antisense = antisense + antinucleotides[randint(0,3)]
                else:
                    antisense = antisense + antinucleotides[2]
            else:
                if randint(0,1000) <= errorpercent:
                    antisense = antisense + antinucleotides[randint(0,3)]
                else:
                    antisense = antisense + antinucleotides[3]
        readsfile.write(">Read{}[{}:{}]\n{}\n\n".format(counter,start,finish,antisense))
        counter += 1

    readsfile.flush()
    readsfile.close()
    return templateList

def findGaps(start, finish):
    for i in range(start, finish):
        templateList[i] = "X"
    return templateList

def askErrors(sequenceType):
    e = ''
    #loop until user gives a proper 'y' or 'n' answer
    while e != 'y' or 'n':
        e = input('\nTurn on random sequencing errors? Y or N:').lower()
        if e == 'y':
            return True
        if e == 'n':
            return False

def askType():
    ans = input("Which sequencing method would you like to use?\n"
                "Input 1 for Pyrosequencing, 2 for Illumina, and 3 for Sanger:")
    return ans

def askCoverage(sequenceType):
    coverage = float(input("What is the coverage of your sequence protocol?\n"
                            "Read number will be rounded to nearest integer answer:"))
    if sequenceType == 1:
        if coverage < (700/len(template)):
            coverage = coverageErrorLoop(700/len(template))
            readNum = int((coverage*len(template))/700)   #700 = average read size for Pyrosequencing
            return readNum
        else:
            readNum = int((coverage*len(template))/700)
            return readNum
    elif sequenceType == 2:
        if coverage < (50/len(template)):
            coverage = coverageErrorLoop(50/len(template))
            readNum =  int((coverage*len(template))/175)   #175 = average read size for Illumina
            return readNum
        else:
            readNum =  int((coverage*len(template))/175)   #175 = average read size for Illumina
            return readNum
    else:
        if coverage < (400/len(template)):
            coverage = coverageErrorLoop(400/len(template))
            readNum = int((coverage*len(template))/650)    #650 = average read size for Sanger
            return readNum
        else:
            readNum = int((coverage*len(template))/650)    #650 = average read size for Sanger
            return readNum

def coverageErrorLoop(min):
    coverage = 0
    while coverage < min:
        coverage = float(input("Value entered is less than minimum coverage for this sequencing type."
                        "Enter another value:"))
    return coverage

print("Welcome to the Genome Sequencer 5000!")

sequenceType = askType()

errors = askErrors(sequenceType)

readNum = askCoverage(sequenceType)

print('Sequencing has begun!')

if errors == True:
    templateList = randsequence(sequenceType, readNum, template, preamble, antinucleotides, nucleotides)
else:
    templateList = sequence(sequenceType, readNum, template, preamble, antinucleotides, nucleotides)

starts = []
ends = []
sequence_start = None
for i, letter in enumerate(templateList):
    if letter is "X" and sequence_start is None:
        sequence_start = i
    elif letter is not "X" and sequence_start is not None:
        starts.append(sequence_start)
        ends.append(i-1)
        sequence_start = None
if sequence_start is not None: # string ended with an X
    starts.append(sequence_start)
    ends.append(len(sequence)-1)


print(starts)
print(ends)

