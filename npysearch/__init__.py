"""
Please refer to the documentation provided in the README.md,
which can be found at npysearch's PyPI URL: https://pypi.org/project/npysearch/
"""

__all__ = ["blast", "read_fasta", "write_fasta", "cigar_string"]


import os
from time import localtime, strftime
from csv import reader
from _npysearch import *


def read_fasta(filepath):
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
    with open(filepath, "r") as f:
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


def write_fasta(filepath, sequences, wrapAfter=0):
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
    assert (
        isinstance(wrapAfter, int) and wrapAfter >= 0
    ), "wrapAfter must be a whole number of type int."

    with open(filepath, "w") as f:
        # No Wrapping
        if wrapAfter == 0:
            for sequence_name in sequences.keys():
                f.write(">" + sequence_name + "\n")
                f.write(sequences[sequence_name] + "\n")
        # Wrapping
        else:
            for sequence_name in sequences.keys():
                f.write("> " + sequence_name + "\n")
                indices = list(range(0, len(sequences[sequence_name], wrapAfter))) + [
                    -1
                ]
                for i in range(len(indices) - 1):
                    f.write(
                        sequences[sequence_name][indices[i] : indices[i + 1]] + "\n"
                    )

    return None


def blast(
    query,
    database,
    maxAccepts=1,
    maxRejects=16,
    minIdentity=0.75,
    alphabet="nucleotide",
    strand="both",
    outputToFile=False,
):
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
        write_fasta(queryPath, query)

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
        write_fasta(databasePath, database)

    else:
        raise TypeError("database must be of type string or dict")

    outputPath = "output_" + startTime + ".csv"

    if alphabet == "nucleotide":
        dna_blast(
            queryPath,
            databasePath,
            outputPath,
            maxAccepts,
            maxRejects,
            minIdentity,
            strand,
        )
    elif alphabet == "protein":
        protein_blast(
            queryPath, databasePath, outputPath, maxAccepts, maxRejects, minIdentity
        )

    # Delete query, database, and output files, if they were
    # constructed in this function
    if type(query) == dict:
        os.remove(queryPath)
    if type(database) == dict:
        os.remove(databasePath)

    if outputToFile:
        # Adding header to the output csv file from nsearch
        with open(outputPath, "r+") as f:
            old = f.read()
            f.seek(0)
            f.write("QueryId,TargetId,QueryMatchStart,QueryMatchEnd,TargetMatchStart,TargetMatchEnd,QueryMatchSeq,TargetMatchSeq,NumColumns,NumMatches,NumMismatches,NumGaps,Identity,Alignment\n"+old)

        return outputPath
    else:
        table = read_csv(outputPath)
        os.remove(outputPath)
        return table


def cigar_string(query, target):
    """
    Generates cigar string from query and target sequences, with '-'s
    for gaps in alignment

    Input
    -----
    query     = str, with '-' for gap in alignment
    target    = str, same size as query, with '-' for gap in alignment

    Output
    ------
    cigar     = str, cigar string with integer (count) followed by one
                of 'X' (mutation), 'D' (gap), or '=' (match)
    """

    # Ensuring that query and target lengths are equal
    assert len(query) == len(
        target
    ), "Query and Target strings need to be of the same size."

    # Checking of mismatches
    alignment = "".join(
        ["|" if query[i] == target[i] else " " for i in range(len(query))]
    )

    # Modifying alignment so that it contains - for gaps instead of
    # blank space
    alignment = "".join(
        [alignment[i] if query[i] != "-" else query[i] for i in range(len(alignment))]
    )
    alignment = "".join(
        [alignment[i] if target[i] != "-" else target[i] for i in range(len(alignment))]
    )

    flag = alignment[0]
    counter = 0
    matcher = {"|": "=", " ": "X", "-": "D"}
    cigar = ""
    for i, character in enumerate(alignment):
        if character == flag:
            counter += 1
        else:
            cigar = cigar + str(counter) + matcher[flag]
            flag = character
            counter = 1
    cigar = cigar + str(counter) + matcher[flag]

    return cigar


def read_csv(filepath):
    """
    Simple CSV reader function required for blast function to output
    results as a dictionary

    Input
    -----
    filepath   = str, path to the CSV file containing blast results

    Output
    ------
    dictionary = dict, keys = column names, values = columns.
                 Contains 5 str, 8 int, and 1 float columns.
    """

    with open(filepath, "r") as f:
        # Column names
        header = ['QueryId', 'TargetId', 'QueryMatchStart', 'QueryMatchEnd',
                  'TargetMatchStart', 'TargetMatchEnd', 'QueryMatchSeq',
                  'TargetMatchSeq', 'NumColumns', 'NumMatches',
                  'NumMismatches', 'NumGaps', 'Identity', 'Alignment']

        # Rest of the results as a list of lists
        data = [element for line in f.readlines() for element in list(reader([line], delimiter=",", quotechar='"'))]
        # Transposing the list of lists
        data = list(map(list, zip(*data)))

        # Converting type for the int and float columns
        intIndices = [2, 3, 4, 5, 8, 9, 10, 11]
        data = [
            list(map(int, column)) if i in intIndices else column
            for i, column in enumerate(data)
        ]
        floatIndices = [12]
        data = [
            list(map(float, column)) if i in floatIndices else column
            for i, column in enumerate(data)
        ]

        # Creating the dictionary
        dictionary = dict(zip(header, data))

    return dictionary
