#!/usr/bin/env python
## Takes an input file of sequences and translates them in the way defined by an input
## dictionary of the form 'from\tto'.
## Example inFPDict file:
## A	T
## G 	C
## ATGC 	TGCA
## Code annotated by Joe.
import warnings
from CisRegModels import MYUTILS
import sys
import argparse
## See my comments in mergeSeqsByBowtie.py for help on argparser.
parser = argparse.ArgumentParser(description='Translates the contents of a file from one thing to another using a dictionary.')
parser.add_argument('-is',dest='inFPSeqs',	metavar='<inFileSeqs>',help='Input file of sequences to translate.', required=True);
parser.add_argument('-id',dest='inFPDict',	metavar='<inFileDictionary>',help='Input file of translations to use in the form from\tto', required=True);
parser.add_argument('-o',dest='outFP', metavar='<outFile>',help='Where to output results [default=stdout]', required=False);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output errors/warnings [default=stderr]', required=False);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);

args = parser.parse_args();

inFileDict=MYUTILS.smartGZOpen(args.inFPDict,'r'); ## Open the file containing sequences to translate.
inFileSeqs=MYUTILS.smartGZOpen(args.inFPSeqs,'r'); ## Open the file containing the standard translations.

## Where to ouptut errors/warnings
if (args.logFP is not None): ## if the user has specified a log file.
	logFile=MYUTILS.smartGZOpen(args.logFP,'w'); ## write the log file.
	sys.stderr=logFile; ## stderr is where the interpreter's prompts (and their error messages) go.

## Where to output results.
if (args.outFP is None): ## if outFP hasn't been specified by the user, the results are just outputted to stdout.
	outFile= sys.stdout; ## stdout is used for the output of print statements and expressions
else: ## if outFP has been specified, the results are outputted there.
	if args.verbose > 0:
		sys.stderr.write("Outputting to file " + args.outFP + "\n");
	outFile = MYUTILS.smartGZOpen(args.outFP,'w');

translationDict = {};
#raise Exception("Reached bad state=%d for '%s.%d' '%s' at line '%s'" %(state,mid,ver,tfid,line));
for line in inFileDict: ## Iterating over the translations in the inFPDict input file.
	if line is None or line == "" or line[0]=="#":
		continue ## Skip the line, go to next iteration in for loop.
	data = line.rstrip().split("\t"); ## 'string.rstrip([chars])' removes whitespace from the right unless
									  ## characters are specified in the method argument.
									  ## Then '.split(separator)' returns a list of strings after breaking
									  ## the given string by the specified separator.
	translationDict[data[0]] = data[1]; ## The example inFPDict file given at top of script would give
										## a dictionary: translationDict = {'A' : 'T'} after the first iteration,
										## translationDict = {'A' : 'T', 'G' : 'C'} after the second iteration,
										## and translationDict = {'A' : 'T', 'G' : 'C', 'ATGC' : 'TGCA'}
inFileDict.close(); ## Close the file.

for line in inFileSeqs: ## Iterating over the sequences to be translated in the inFPSeqs input file.
	if line is None or line == "" or line[0]=="#":
		continue ## Skip the line
	seq=line.rstrip(); ## Remove whitespace from the right of the string.
	if seq in translationDict: ## if the sequence to be translated appears in the dictionary
		outFile.write(translationDict[seq]+"\n"); ## Write the translation of that sequence to an output file.
	else:
		outFile.write(seq + "\n"); ## if it's not in the dictionary, just write the original sequence to file.
inFileSeqs.close();

outFile.close();

if (args.logFP is not None): ## If a log file was specified, it needs to be closed. 
	logFile.close();
