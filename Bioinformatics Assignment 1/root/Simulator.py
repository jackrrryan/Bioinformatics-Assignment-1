'''
Created on Sep 19, 2014

@author: Jack Ryan B00321033
'''

from random import randint      #To be used in generating random read lengths for Illumina and Sanger sequencing,
                                #as well as generating errors in random sequencing

file = open('Assnt1_sampleinput.fna.txt', 'r')  #opening the FASTA file
template = file.read()  #reading the FASTA file to a string

#some variables used in all calculations
preamble = template[0:72]   #This is all the text before the start of the sequence
template = template[73:]    #This is the actual sequence

def make_templateList(template):        #Makes a list of bases out of the template strand
    template_with_gaps = []
    for i in range(0,len(template)):
        template_with_gaps.append(template[i])
    return template_with_gaps

def readSizePosition(sequenceType):     #Find the start and end of the current read by index
    if sequenceType == 1:      #1 stands for pyrosequencing
        readLength = 700  #700bp is length of pyrosequencing reads
        start = randint(0,len(template)-readLength) #read must end before last bp in sequence
        finish = start+readLength
    elif sequenceType == 2:     #2 stands for Illumina sequencing
        readLength = randint(50,300)    #Illumina can have reads between 50-300 bp, so generate random number
        start = randint(0,(len(template)-readLength))
        finish = start+readLength
    else:       #3 is the only other option, and that stands for Sanger sequencing
        readLength = randint(400,900)   #Sanger can have reads between 400-900 bp, so generate random number
        start = randint(0,(len(template)-readLength))
        finish = start+readLength
    return (start,finish)       #return a tuple of start and end values for this read

def sequence(sequenceType, readNum, template, templateList, preamble, antinucleotides):     #sequence the current read
    
    readsfile = open('reads.txt', 'w')      #Creating the reads.txt file, where all of the reads will be stored
    readsfile.write(preamble + '\n\nBegin sequencing:\n\n')    #Write the preamble to the reads file

    counter = 0     #counter used to mark in reads file the read number (ex. [Read12])

    for i in range(0,readNum):      #loop through the reads
        antisense = []              #Create a blank antisense list that stores the antisense nucleotides for each read
        sfpositions = readSizePosition(int(sequenceType))       #tuple for start and end of the read returned
        start = sfpositions[0]    # this is the start position of this specific read
        finish = sfpositions[1]   # this is the end position of this specific read

        for j in range(start,finish):   #loop through the length of the read
            templateN = template[j]     #value for the current nucleotide
            if templateN == 'A':
                antisense.append(antinucleotides[0])    #add the antisense nucleotide to the list
            elif templateN == 'T':
                antisense.append(antinucleotides[1])
            elif templateN == 'C':
                antisense.append(antinucleotides[2])
            else:
                antisense.append(antinucleotides[3])
        antisense = ''.join(antisense)      #turn the antisense list into a string
        readsfile.write(">Read{}[{}:{}]\n{}\n\n".format(counter,start,finish,antisense)) #formatting for the read header
        counter += 1
        template_with_gaps = assignGaps(start, finish, templateList)    #find and mark the location of this read in the template strand

    readsfile.flush()
    readsfile.close()
    return template_with_gaps       #return the template_with_gaps list that has reads marked

def randsequence(sequenceType, readNum, template, preamble, antinucleotides):   #sequence the current read with errors turned on (this method is not used in the program as is)

    #determine the proper error percentage for the sequencing mode.
    #This is done through an x/1000 calculation.
    if sequenceType == 1:
        errorpercent = 1    #1/1000, or 99.9% for Pyrosequencing
    elif sequenceType == 2:
        errorpercent = 20   #20/1000, or 98% for Illumina
    else:
        errorpercent = 1    #1/1000, or 99.9% for Sanger

    readsfile = open('reads.txt', 'w') #Creating the reads.txt file, where all of the reads will be stored
    readsfile.write(preamble + '\n\nBegin sequencing:\n\n') #Write the preamble to the reads file

    counter = 0 #counter used to mark in reads file the read number (ex. [Read12])

    for i in range(0,readNum):  #loop through the reads
        antisense = []          #Create a blank antisense list that stores the antisense nucleotides for each read
        sfpositions = readSizePosition(int(sequenceType))       #tuple for start and end of the read returned
        start = sfpositions[0]    # this is the start position of this specific read
        finish = sfpositions[1]   # this is the end position of this specific read

        for j in range(start,finish):       #loop through the length of the read
            templateN = template[j]         #value for the current nucleotide
            if templateN == 'A':
                if randint(0,1000) <= errorpercent:     #generate a random number from 0-1000, if it falls below the error percentage then turn on errors for this nucleotide
                    antisense = antisense + antinucleotides[randint(0,3)]       #generate a random # between 0-3 and search through antinucleotide list
                else:
                    antisense = antisense + antinucleotides[0]  #if error check fails, put correct nucleotide in
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
        antisense = ''.join(antisense)      #turn the antisense list into a string
        readsfile.write(">Read{}[{}:{}]\n{}\n\n".format(counter,start,finish,antisense)) #formatting for the read header
        counter += 1
        template_with_gaps = assignGaps(start, finish) #find and mark the location of this read in the template strand

    readsfile.flush()
    readsfile.close()
    return template_with_gaps

def assignGaps(start, finish, templateList):    #Mark the locations of the completed reads
    template_with_gaps = templateList       #copy the template to a new list
    for i in range(start, finish):      #run through the sequence
        template_with_gaps[i] = "X"     #add the "X" nucleotide where a read has taken place
    return template_with_gaps

def getNumOfReads(sequenceType, coverage):       #find the number of reads for this coverage value

    #n = CL/G

    if sequenceType == 1:       #pyrosequencing
        readNum = int((coverage*len(template))/700) #700 = average read size for Pyrosequencing
        return readNum
    elif sequenceType == 2:     #illumina sequencing
        readNum = int((coverage*len(template))/175)   #175 = average read size for Illumina
        return readNum
    else:                       #Sanger sequencing
        readNum = int((coverage*len(template))/650)    #650 = average read size for Sanger
        return readNum

def findnumofgaps(template_with_gaps):      #find the number of gaps in this sequence
    starts = []     #blank starts list to hold locations of starts of contigs
    ends = []       #blank ends list to hold locations of ends of contigs
    sequence_start = None       #start of a contig
    for i, letter in enumerate(template_with_gaps):
        if letter is "X" and sequence_start is None:
            sequence_start = i
        elif letter is not "X" and sequence_start is not None:
            starts.append(sequence_start)
            ends.append(i-1)
            sequence_start = None
    if sequence_start is not None: # string ended with an X
        starts.append(sequence_start)
        ends.append(len(template)-1)
    num_of_gaps = len(starts)+1
    return num_of_gaps

def runSimulator():

    antinucleotides = ['T','A','G','C']     #the four antinucleotides to the order A,T,C,G

    resultsfile = open('OtherResults.txt', 'w') #open t

    sequence_and_coverage_list = ['Pyrosequencing, 1', 'Pyrosequencing, 5', 'Pyrosequencing, 10',
                                  'Illumina, 1', 'Illumina, 5', 'Illumina, 10',
                                  'Sanger, 1', 'Sanger, 5', 'Sanger, 10']
    average_list = []

    num_of_gaps_list = []
    num_of_times_greater_than_2_list = []
    greater_gaps_counter = 0

    coverage_values = [0.5, 1, 5]

    templateList = make_templateList(template)   #this is a copy of the template string to be used in determining gaps in sequencing

    for i in range(1,4):
        for j in range(0,3):
            for k in range(0, 1000):
                readNum = getNumOfReads(i, coverage_values[j])
                template_with_gaps = sequence(i, readNum, template, templateList, preamble, antinucleotides)
                num_of_gaps_list.append(findnumofgaps(template_with_gaps))
                if findnumofgaps(template_with_gaps) > 2:
                    greater_gaps_counter += 1
            num_of_times_greater_than_2_list.append(greater_gaps_counter)
            greater_gaps_counter = 0
            average = float(sum(num_of_gaps_list)/len(num_of_gaps_list))
            print(average)
            average_list.append(average)
            num_of_gaps_list.clear()

    for i in range(0,9):
        print('{}: # times > 2: {} average #gaps: {}'.format(sequence_and_coverage_list[i], num_of_times_greater_than_2_list[i], average_list[i]))




print(runSimulator())


