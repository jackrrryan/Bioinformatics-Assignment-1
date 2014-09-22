'''
Created on Sep 19, 2014

@author: Jack
'''

from random import randint

def readsequence(template, antinucleotides, nucleotides, errors, start, finish):
    
    antisense = ''

    for i in range (start, finish):
        templateN = template[i:i+1]
        randomerror = randint(0,1000)   #from reference. Error rate determined to be 99.5%
        if errors == True:
            randomerror = 1000
        if templateN == 'A':
            if randomerror < 5:     #5/1000 = 0.5%
                antisense = antisense + antinucleotides[randint(0,3)]
            else:
                antisense = antisense + antinucleotides[0]
        elif templateN == 'T':
            if randomerror < 5:
                antisense = antisense + nucleotides[randint(0,3)]
            else:
                antisense = antisense + antinucleotides[1]
        elif templateN == 'C':
            if randomerror < 5:
                antisense = antisense + nucleotides[randint(0,3)]
            else:
                antisense = antisense + antinucleotides[2]
        else:
            if randomerror < 5:
                antisense = antisense + nucleotides[randint(0,3)]
            else:
                antisense = antisense + antinucleotides[3]
    
    return antisense

file = open('Assnt1_sampleinput.fna.txt', 'r')

template = file.read()

preamble = template[0:72]
print(preamble)

template = template[73:]
print(template[:1000])

print("Welcome to the Genome Sequencer 5000!\nThis sequencer uses Pyrosequencing.")

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

pyroMaterials = ['primer','DNApolymerase','ATPsulfurylase','luciferase','apyrase','APS','luciferin']
nucleotides = ['A','T','C','G']
antinucleotides = ['T','A','G','C']

print('Sequencing has begun!')

fragments = []

for i in range(0, len(template)):
    start = randint(0,len(template))
    fragments.append(readsequence(template, antinucleotides, nucleotides,
                               errors, start, randint(start,len(template))))

print(fragments[0])
