#!/usr/bin/python3
import re

filepath = "Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"
newFile = open("Sus_scrofa.Sscrofa11.1.dna.toplevel.Ymasked.fa","w+")

with open(filepath, 'r') as oldFile:
    checkpoint1 = False  # checkpoint is false to begin with, only true when hits the Y chr
    for line in oldFile:  # loop through lines in fasta file
        if line.startswith(">Y"):
            checkpoint1 = True  # checkpoint1 is true when you reach the Y chr
        elif line.startswith(">"):
            checkpoint1 = False  # turn off replacements for lines that start with >, but not >Y

        if checkpoint1:
            newLine = re.sub(r"[ACGT]", "N", line)  # modify line
        else:
            newLine = line  # don't modify line
        newFile.write(newLine)
newFile.close()
