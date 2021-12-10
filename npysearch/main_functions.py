import os
from time import localtime, strftime

from _npysearch import *
from .auxillary_functions import *


def readFasta(filepath):
    """
    Simple FASTA file reader

    Input
    -----
    filepath  = str, path to the fasta file to be read
    
    Output
    ------
    sequences = dict, keys = sequence id, values = sequence. Both keys
                and values are strings.
    """

    # Ensure file exists
    if not os.path.isfile(filepath):
        raise IOError("File does not exist.")

    sequences = {}
    with open(filepath, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.strip()[0] == ">":
                current_sequence = line.strip()[1:]
                sequences[current_sequence] = ""
            elif line.strip() == "":
                continue
            else:
                sequences[current_sequence] += line.strip()

    return sequences

def writeFasta(filepath, sequences, wrapAfter = 0):
    """
    Simple FASTA file writer. Use wrapAfter to have sequences wrap 
    after a given number of characters.
    
    Input
    -----
    filepath  = str, path to the fasta file to be written
    sequences = dict, keys = sequence id, values = sequence. Both keys
                and values are strings.
    wrapAfter = int, whole number indicating number of characters in 
                each line of sequence in the file. 
                0 indicates no wrapping. (Default = 0)
    
    Output
    ------
    None
    """
    
    # Check that wrapAfter is a whole number
    assert isinstance(wrapAfter, int) and wrapAfter >= 0, "wrapAfter must be a whole number of type int."

    with open(filepath, 'w') as f:
        # No Wrapping
        if wrapAfter == 0:
            for sequence_name in sequences.keys():
                f.write("> " + sequence_name + "\n")
                f.write(sequences[sequence_name] + "\n")
        # Wrapping
        else:
            for sequence_name in sequences.keys():
                f.write("> " + sequence_name + "\n")
                indices = list(range(0, len(sequences[sequence_name], wrapAfter))) + [-1]
                for i in range(len(indices)-1):
                    f.write(sequences[sequence_name][indices[i]:indices[i+1]] + "\n")

    return None


def blast(query, database, maxAccepts = 1, maxRejects = 16, 
          minIdentity = 0.75, alphabet = "nucleotide", strand = "both",
          outputToFile = False):
    """
    Runs BLAST sequence comparison algorithm

    Input
    -----
    query        = dict or str. Either a dictionary of sequences with 
                   sequence ids as keys and sequences as values, or a
                   path str to the fasta file containing the sequences
    database     = dict or str. Either a dictionary of sequences with
                   sequence ids as keys and sequences as values, or a
                   path str to the fasta file containing the sequences
    maxAccepts   = int, number specifying the maximum accepted hits 
                   (Default = 1)
    maxRejects   = int, number specifying the maximum rejected hits 
                   (Default = 16)
    minIdentity  = float, number specifying the minimal accepted 
                   sequence similarity between the query and database
                   sequences (Default = 0.75)
    alphabet     = str, "nucleotide" or "protein" to specify the query
                   and database alphabet (Default = "nucleotide")
    strand       = str, specify the strand to search: "plus", "minus",
                   or "both". Only affects nucleotide searches. 
                   (Default = "both")
    outputToFile = boolean, set to True to get the results table as a
                   csv file in the working directory and False to 
                   return the results as a dictionary of lists
                   (Default = False)
    
    Output
    ------
    table        = dict of lists, results table in the form of a dict
                   of lists with column names as keys and columns as 
                   values. Contains 5 str, 8 int, and 1 float columns.
                   Can be converted easily to a pandas dataframe using
                   pandas.DataFrame.from_dict()
        OR

    csvPath      = str, path to the csv file containing the results, 
                   stored in the working directory
    """

    startTime = strftime("%Y-%m-%d-%H:%M:%S", localtime())
    
    # Checking what form the query was input in
    # str for path to fasta file and dict for sequences
    # If dict, write to file
    if type(query) == str:
        queryPath = query
        # Ensure file exists
        if not os.path.isfile(queryPath):
            raise IOError("Query file does not exist.")

    elif type(query) == dict:
        queryPath = "query_" + startTime + ".fasta"
        writeFasta(queryPath, query)
    
    else:
        raise TypeError("query must be of type string or dict")

    # Checking what form the database was input in
    # str for path to fasta file and dict for sequences
    # If dict, write to file
    if type(database) == str:
        databasePath = database
        # Ensure file exists
        if not os.path.isfile(databasePath):
            raise IOError("Database file does not exist.")

    elif type(database) == dict:
        databasePath = "database_" + startTime + ".fasta"
        writeFasta(databasePath, database)
    
    else:
        raise TypeError("database must be of type string or dict")

    outputPath = "output_" + startTime + ".txt"

    if alphabet == "nucleotide":
        dna_blast(queryPath, databasePath, outputPath, maxAccepts, maxRejects, minIdentity, strand)
    elif alphabet == "protein":
        protein_blast(queryPath, databasePath, outputPath, maxAccepts, maxRejects, minIdentity)

    csvPath = "output_" + strftime("%Y-%m-%d-%H:%M:%S", localtime()) + ".csv"
    writeCSV(outputPath, csvPath)

    # Delete query, database, and output files, if they were 
    # constructed in this function
    if type(query) == dict:
        os.remove(queryPath)
    if type(database) == dict:
        os.remove(databasePath)

    os.remove(outputPath)

    if outputToFile:
        return csvPath
    else:
        table = readCSV(csvPath)
        os.remove(csvPath)
        return table 
