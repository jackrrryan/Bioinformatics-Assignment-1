'''
Created on Sep 19, 2014

@author: Jack Ryan B00321033
'''

from random import randint

#some variables used in all calculations
nucleotides = ['A','T','C','G']
antinucleotides = ['T','A','G','C']
preamble = template[0:72]   #This is all the text before the start of the sequence
template = template[73:]    #This is the actual sequence

def pyrosequencing(template, antinucleotides, nucleotides, errors, start, finish):
    
    antisense = ''

    for i in range (start, finish):
        templateN = template[i:i+1]
        randomerror = randint(0,1000)   #from reference. Error rate determined to be 99.9%
        if errors == True:
            randomerror = 1000
        if templateN == 'A':
            if randomerror < 1:     #1/1000 = 0.1%
                antisense = antisense + antinucleotides[randint(0,3)]
            else:
                antisense = antisense + antinucleotides[0]
        elif templateN == 'T':
            if randomerror < 1:
                antisense = antisense + nucleotides[randint(0,3)]
            else:
                antisense = antisense + antinucleotides[1]
        elif templateN == 'C':
            if randomerror < 1:
                antisense = antisense + nucleotides[randint(0,3)]
            else:
                antisense = antisense + antinucleotides[2]
        else:
            if randomerror < 1:
                antisense = antisense + nucleotides[randint(0,3)]
            else:
                antisense = antisense + antinucleotides[3]
    
    return antisense

def askErrors():
    e = ''
    errors = False
    #loop until user gives a proper 'y' or 'n' answer
    while e != 'y' or 'n':
        e = input('\nTurn on random sequencing errors? Y or N:').lower()
        if e == 'y':
            errors = True
            break
        if e == 'n':
            errors = False
            break

def askType():
    ans = ''
    while ans != 1 or 2 or 3:
        ans = input("Which sequencing method would you like to use?\n"
                    "Input 1 for Pyrosequencing, 2 for Illumina, and 3 for Sanger:")
        if ans == 1:
            return 1
        elif ans == 2:
            return 2
        elif ans == 3:
            return 3

def askCoverage(sequenceType):
    while !coverage.isInstance(coverage, double)
        coverage = input("What is the approximate coverage of your sequence protocol (in percentage)?:")


file = open('Assnt1_sampleinput.fna.txt', 'r')  #opening the FASTA file
template = file.read()  #reading the FASTA file to a string

print("Welcome to the Genome Sequencer 5000!")

errors = askErrors()

sequenceType = askType()

coverage = askCoverage(sequenceType)

print('Sequencing has begun!')