#!/usr/bin/env python
"""
The purpose of this file is to parse through files containing promoter sequences and
gather motifs of sequences (re-occuring sequences), indicating abundance.
The promoter abundances are subsequently ranked in descending order.
"""
## Code annotated by K
import warnings
from CisRegModels import MYUTILS
import sys
import argparse

## Parses through a file containing many promoter sequnces and
## puts them in descending order of abundance.
## Building options for arguments to be passed.
## If contains 'required=True' the argument and its flag (i.e. '-i') must be passed.
## If contains 'required=False' the argument and its flag
## (i.e. '-o, -l, -v') are not necessarily required.
## Can execute 'python collapsePromoters.py --help' from command line to see options.
## Supporting information at https://docs.python.org/3/library/argparse.html
parser = argparse.ArgumentParser(description='Takes a file containing many promoter sequences and counts them, outputting them in decreasing order of abundance.')
parser.add_argument('-i',dest='inFP',	metavar='<inFile>',help='Input file of promoter sequences', required=True);
parser.add_argument('-o',dest='outFP', metavar='<outFile>',help='Where to output results [default=stdout]', required=False);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output errors/warnings [default=stderr]', required=False);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);

## initialize parser
args = parser.parse_args();

## Uses the smartGZOpen function in MYUTILS to read/parse through a file
## (indicated by 'r' argument).
inFile=MYUTILS.smartGZOpen(args.inFP,'r');

## initialize counter
promoterCounts = {};

## Creates log file of errors/warnings
## logFile = flexible framework to emit log messages
## (logging = tracking events when software runs)
## 'w' indicates file writing
if (args.logFP is not None):
	logFile=MYUTILS.smartGZOpen(args.logFP,'w');
	sys.stderr=logFile;

## Creates output directions (inculding warnings)
if (args.outFP is None):
	## system specific function - standard output
	outFile= sys.stdout;
else:
	if args.verbose>0: warnings.warn("Outputting to file "+args.outFP);
		## same as log file for errors/warnings
	outFile = MYUTILS.smartGZOpen(args.outFP,'w');

#raise Exception("Reached bad state=%d for '%s.%d' '%s' at line '%s'" %(state,mid,ver,tfid,line));
## Test if line can be read
## Adjust counter accordingly
for line in inFile:
	if line is None or line == "" or line[0]=="#":
		## If the above is valid, skip passed everything and restart the for loop
		continue
	## Returns a copy of the string with trailing characters removed
		## Based on argument passed
	line = line.rstrip();
	## Physical counter to add up abundance values as motifs occur per line
	if line not in promoterCounts:
		promoterCounts[line]=1;
	else:
		promoterCounts[line]+=1;
## Closes file currently associated with the output
inFile.close();

## Outputs counters (of abundances) in reversed, sorted order (descending)
for i in reversed(sorted(promoterCounts, key=promoterCounts.__getitem__)):
	outFile.write("%s\t%i\n"%(i, promoterCounts[i]));
## Close output
outFile.close();
## Close log if no output
if (args.logFP is not None):
	logFile.close();
