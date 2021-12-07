from _npysearch import *
import os
from time import localtime, strftime

def test():
    print("This function works.")

    return None

def cigarString(query, alignment, target):

    # Modifying alignment so that it contains - for gaps instead of blank space
    alignment = "".join([alignment[i] if query[i] != "-" else query[i] for i in range(len(alignment))])
    alignment = "".join([alignment[i] if target[i] != "-" else target[i] for i in range(len(alignment))])

    flag = alignment[0]
    counter = 0
    matcher = {"|" : "=", " " : "X", "-" : "D"}
    cigar = ""
    for i, character in enumerate(alignment):
        if character == flag:
            counter +=1
        else:
            cigar = cigar + str(counter) + matcher[flag]

            flag = character
            counter = 1
    cigar = cigar + str(counter) + matcher[flag]
        
    return cigar

def readFasta(filepath):
    if not os.path.isfile(filepath):
        raise IOError("File does not exist.")

    sequences = {}
    with open(filepath, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[0] == ">":
                current_sequence = line.strip()[1:]
                sequences[current_sequence] = ""
            else:
                sequences[current_sequence] += line.strip()

    return sequences

def writeFasta(filepath, sequences):

    with open(filepath, 'w') as f:
        for sequence_name in sequences.keys():
            f.write("> " + sequence_name)
            f.write(sequences[sequence_name])

    return None

def writeCSV(filepath, csvPath):
    header = ["QueryId", "TargetId", "QueryMatchStart",
              "QueryMatchEnd", "TargetMatchStart", "TargetMatchEnd",
              "QueryMatchSeq", "TargetMatchSeq", "NumColumns", "NumMatches",
              "NumMismatches", "NumGaps", "Identity", "Alignment"]

    with open(csvPath, "w") as csvFile:
        csvFile.write(",".join(header) + "\n")
        with open(filepath, "r") as f:
            row = [None] * 14
            for i, line in enumerate(f):
                if i%13==4:
                    row[0] = line.strip().split()[2][1:].strip()
                    row[1] = line.strip().split()[2][1:].strip()
                if i%13==7:
                    row[2] = line.strip().split()[1].strip()
                    row[3] = line.strip().split()[-1].strip()
                    row[6] = line.strip().split()[3].strip()
                    start = line.index("+") + 2
                if i%13==8:
                    alignment = line[start:].strip("\n")
                if i%13==9:
                    row[4] = line.strip().split()[1].strip()
                    row[5] = line.strip().split()[-1].strip()
                    row[7] = line.strip().split()[3].strip()
                if i%13==11:
                    row[8] = line.split()[0].strip()
                    row[9] = line.split()[2].strip()
                    row[10] = str(int(row[8]) - int(row[9]))
                    row[11] = line.split()[5].strip()
                    row[12] = str(float(line.split()[4][1:-3])/100)
                    # row[13] = "-"
                    row[13] = cigarString(query = row[6], alignment = alignment, target = row[7])
                if i%13==12:
                    csvFile.write(",".join(row) + "\n")
                    row = [None] * 14

    return None

def readCSV(filepath):
    with open(filepath, "r") as f:
        header = f.readline().strip().split(",")
        
        data = [line.strip().split(",") for line in f.readlines()]
        data = list(map(list, zip(*data))) # Transposing list of lists

        intIndices = [2,3,4,5,8,9,10,11]
        data = [list(map(int, column)) if i in intIndices else column for i,column in enumerate(data)]
        floatIndices = [12]
        data = [list(map(float, column)) if i in floatIndices else column for i,column in enumerate(data)]

        dictionary = dict(zip(header, data))

    return dictionary



def blast(query,
          database,
          maxAccepts = 1,
          maxRejects = 16,
          minIdentity = 0.75,
          alphabet = "nucleotide",
          strand = "both",
          pathsUsed = False,
          outputToFile = False):
    
    startTime = strftime("%Y-%m-%d-%H:%M:%S", localtime())
    
    if pathsUsed:
        queryPath = query
        databasePath = database
    else:
        queryPath = "query_" + startTime + ".fasta"
        databasePath = "database_" + startTime + ".fasta"

        writeFasta(queryPath, query)
        writeFasta(databasePath, database)

    outputPath = "output_" + startTime + ".txt"

    if alphabet == "nucleotide":
        dna_blast(queryPath, databasePath, outputPath, maxAccepts, maxRejects, minIdentity, strand)
    elif alphabet == "protein":
        protein_blast(queryPath, databasePath, outputPath, maxAccepts, maxRejects, minIdentity)

    csvPath = "output_" + strftime("%Y-%m-%d-%H:%M:%S", localtime()) + ".csv"
    writeCSV(outputPath, csvPath)

    # Delete query, database, and output files
    if not pathsUsed:
        os.remove(queryPath)
        os.remove(databasePath)
    os.remove(outputPath)

    if outputToFile:
        return csvPath
    else:
        os.remove(csvPath)
        return readCSV(csvPath)






    