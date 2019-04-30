"""The purpose of this file is to parse through sequences and convert them into binary representation."""

#!/usr/bin/env python
import warnings
from  CisRegModels import MYUTILS
import sys
import argparse

## Parses through a file containing sequences and converts them into
## binary representation, excluding non-[ATCG] characters.
## Outputs one line per sequence.
## Building options for arguments to be passed.
## If contains 'required=True' the argument and its flag (i.e. '-i, -m') must be passed.
## If contains 'required=False' the argument and its flag
## (i.e. '-b, -o, -l, -v') are not necessarily required.
## Can execute 'python seqToOHC.py --help' from command line to see options.
## Supporting information at https://docs.python.org/3/library/argparse.html
parser = argparse.ArgumentParser(description='Converts a set of sequences into a one-hot-code (binary)
representation - excludes non [ATGC] chars.  Output in ACGT order, one line per sequence, base then
position.')
parser.add_argument('-i',dest='inFP',	metavar='<inFile>',help='Input file of sequences with a value in the
second column that will preceed the OHC output on each line, separated by a tab', required=True);
parser.add_argument('-m',dest='maxLen',	metavar='<maxSeqLen>',help='The maximum sequence length to consider
(truncated after this point)', required=True);
parser.add_argument('-b',dest='orientBack', action='count',help='Align sequences of different sizes to back
[default=front]?', required=False, default=0);
parser.add_argument('-o',dest='outFP', metavar='<outFile>',help='Where to output results [default=stdout]',
required=False);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output errors/warnings
[default=stderr]', required=False);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);

## initialize parser
args = parser.parse_args();

## Uses the smartGZOpen function in MYUTILS to read/parse through a file
## (indicated by 'r' argument).
inFile=MYUTILS.smartGZOpen(args.inFP,'r');
## Initialize max length integer
maxSeqLen = int(args.maxLen);

## Creates log file of errors/warnings
## logFile = flexible framework to emit log messages
## (logging = tracking events when software runs)
if (args.logFP is not None):
	logFile=MYUTILS.smartGZOpen(args.logFP,'w');
	sys.stderr=logFile;

## Creates output directions (inculding warnings)
if (args.outFP is None):
        ## system specific function - standard output, if output
	outFile= sys.stdout;
else:
	## Notes on how to exit the program
	if args.verbose>0: sys.stderr.write("Outputting to file "+args.outFP+"\n");
	## Same as log file for errors/warnings - write file
	outFile = MYUTILS.smartGZOpen(args.outFP,'w');

## Initialize sequence values that should be detected
## Ignores non-[ATCG] characters
## Specific order of bases when read
BASES = ['A','C','G','T'];


#raise Exception("Reached bad state=%d for '%s.%d' '%s' at line '%s'" %(state,mid,ver,tfid,line));
## Test if line can be read
## Parse through for specific characters
for line in inFile:
	if line is None or line == "" or line[0]=="#": continue
    	## Returns a copy of the string with trailing characters removed
    	## Based on argument passed (so strip tabs)
	data=line.rstrip().split("\t");
	curSeq = data[0].upper();
	curLabel = data[1];
	curSeqLen = len(curSeq);
    	## How to write string copy with specified characters removed
	outFile.write(curLabel+"\t");
	## Loop through all the values in bases for each line
	for b in BASES:
		if args.orientBack>0:
            		#output the front until we run out of chars, then print 0s, or truncate if it's longer than maxSeqLen
			if curSeqLen > maxSeqLen:
				curSeq = curSeq[:maxSeqLen];
			outFile.write("\t".join([ "1" if curSeq[i]==b else "0" for i in range(0,len(curSeq))]));
			if maxSeqLen > curSeqLen: # print 0s
				outFile.write("\t" + "\t".join(["0"]*(maxSeqLen - curSeqLen)));
        	#print any 0s at the front and also truncate from the front if too long.
		else:
			if maxSeqLen > curSeqLen: # print 0s
				outFile.write("\t".join(["0"]*(maxSeqLen - curSeqLen)) + "\t");
			elif curSeqLen > maxSeqLen:
				curSeq = curSeq[(len(curSeq) - maxSeqLen):];
			outFile.write("\t".join([ "1" if curSeq[i]==b else "0" for i in range(0,len(curSeq))]));
        	## if b does not equal character 'T' (which is the last character in the bases)
		##  print no values (just tab), ie move to the following line
		if b!="T":
			outFile.write("\t");
	## Print new line
	outFile.write("\n");
## Closes file currently associated with the output
inFile.close();
## Close output
outFile.close();

## Close log if no output
if (args.logFP is not None):
	logFile.close();

