

# This program will act as a GFF3 feature exporter. It will take 4 arguments
# related to GFF3 columns and use this information to return the related FASTA
# sequence.
# Author: James Matsumura

import sys,re

# Assign all arguments entered in the command line to a variable. The order in
# which the user enters these arguments is first with a path to a GFF3 file,
# then a type (i.e. gene), then an attribute (i.e. ID), then a value (ID seq).
source_gff = re.search(r"--source_gff=(.*)", sys.argv[1])
sg = source_gff.group(1)

type = re.search(r"--type=(.*)", sys.argv[2])
t = type.group(1)

attribute = re.search(r"--attribute=(.*)", sys.argv[3])
a = attribute.group(1)

value = re.search(r"--value=(.*)", sys.argv[4])
v = value.group(1)

# Combine the key=value pair into one variable.
pair = ("{0}={1}".format(a, v))

file = open(sg)

# Establish booleans to know whether or not a sequence needs to be isolated
# from the FASTA section dependent on how many matches are found and a boolean
# for when a matched FASTA sequence is found.
match_found = False
many_matches = False
found_seq = False

# Establish an empty string to build the sequence upon and a position to
# reference while isolating the sequence.
sequence = ""
pos = 1

# Begin to search through each line of the file searching for a match to the
# query. If a SINGLE match is found, isolate the FASTA sequence at the end of
# the file and output. If multiple or no matches are found, output this
# result to the user.
for line in file:

    # Remove all newline characters.
    line = line.rstrip()

    # Ignore comment lines.
    if line.startswith("#"):
        pass

    # If multiple matches are found, leave the for loop and tell the user.
    elif many_matches == True:
        break

    # Establish the start of the FASTA section. If only one match is found and
    # the sequence hasn't been isolated yet, isolate the sequence.
    elif line.startswith(">") and match_found == True and found_seq == False:
        # Ensure that the seqid matches exactly to the FASTA entry.
        s = re.search(r">(.*)", line)

        if seqid == s.group(1):
            found_seq = True


    # Build a string of the FASTA sequence found between two FASTA entries.
    elif found_seq == True:
        # If another '>' symbol is found the next FASTA entry has been reached
        # so all of the sequence should have been found. End the loop to print
        # the results.
        if line.startswith(">"):
            break

        # Iterate through lines by each character and only add to the string
        # when the characters are between the set start and end positions.
        for char in line:
            if start <= pos and pos <= end:
                sequence += char
                pos += 1

            else:
                pos += 1

    # If not a comment or not in the FASTA part, parse through the feature data
    # to try isolate a match.
    else:

        # Since GFF3 is tab-delimited, split on tabs.
        cols = line.split("\t")

        # If an entry is truly in GFF3 format, it will have 9 columns.
        if len(cols) == 9:

            # If the type and key=value pair is found in an entry then isolate
            # the data within it. Note that the code will only return the FASTA
            # sequence if only ONE match is found.

            # Handle the case where multiple matches are found.
            if cols[2] == t and pair in cols[8] and match_found == True:

                many_matches = True

            # Handle the case where the first match is found.
            elif cols[2] == t and pair in cols[8] and match_found == False:

                # Assign the start, end, and strand variables to those found in
                # this line and use this information to isolate the sequence in
                # the FASTA section. Also isolate the seqid for searching
                # through the FASTA part of the file.
                seqid = cols[0]
                start = int(cols[3])
                end = int(cols[4])
                strand = cols[6]
                match_found = True

# If multiple matches are found, print a warning.
if many_matches == True:
    print("WARNING: More than one feature matches the query.")

# If no matches are found, tell the user.
elif match_found == False:
    print("No features matched the query.")

# Provide the user with the matched feature if found.
elif match_found == True:

    # Create a header for the output that echoes the query input.
    print(">{0}:{1}:{2}".format(t, a, v))

    # Format the sequence printing so that it displays 60 characters per line
    # using the range function. The start is 0, the end is the length of the
    # sequence, and the step is 60 to enable each line to print groups of 60.
    for i in range(0, len(sequence), 60):

        # Print groups of 60 characters. Remember that the first number in the
        # [] will be inclusive while the second number will be exclusive. Thus,
        # the indices printed at each step will be 0-59, 60-119, etc. until the
        # sequence length is reached.
        print(sequence[i:i+60])