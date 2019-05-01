#!/usr/bin/env python
"""
The purpose of this script is to ...
"""
## Code annotated by Guoyao
import warnings
from CisRegModels import MYUTILS
from CisRegModels.DNA import *
import sys
import argparse
## set up the parser command in python to run some command in shell
parser = argparse.ArgumentParser(description='Takes two fastq files for convergent reads and aligns them into one sequence in the strandedness of read1.')
parser.add_argument('-i1',dest='inFP1',	metavar='<inFile>',help='Input file of FastQ for read1', required=True);
parser.add_argument('-i2',dest='inFP2',	metavar='<inFile>',help='Input file of FastQ for read2', required=True);
parser.add_argument('-p',dest='overlap',	metavar='<overlap>',help='Approximately how much overlap should there be?', required=True);
parser.add_argument('-r',dest='overlapRange',	metavar='<overlapRange>',help='How many bases of flexibility should there be in the overlap?', required=True);
parser.add_argument('-o',dest='outFP', metavar='<outFile>',help='Where to output results [default=stdout]', required=False);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='0 [default=stderr]', required=False);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);

## set up the initial
args = parser.parse_args();

## open the files
inFile1=MYUTILS.smartGZOpen(args.inFP1,'r');
inFile2=MYUTILS.smartGZOpen(args.inFP2,'r');

## if the input is not None, use the sys command to direct the location of the inputs
if (args.logFP is not None):
	logFile=MYUTILS.smartGZOpen(args.logFP,'w');
	sys.stderr=logFile;

## set up the details of the output, including the locations and errors
if (args.outFP is None):
	### I don't know what does this mean, I assume it is when the output location is none, just print the name of the output
	outFile= sys.stdout;
else:
	if args.verbose>0: warnings.warn("Outputting to file "+args.outFP);
	outFile = MYUTILS.smartGZOpen(args.outFP,'w');

## the function to read next file, read files by line(only one line), remove line breaks
def getNextRead(inFile):
	name = inFile.readline().rstrip();
	seq = inFile.readline().rstrip();
	plus = inFile.readline().rstrip();
	quality = inFile.readline().rstrip();
	return (name,seq,quality);

## the function to test the alignment is match or not and return with the suitability
def testAlignment(seq1, seq2, offset):
	score = 0.0;
	for i in range(offset, len(seq1)):
		if seq1[i] == seq2[i-offset]:
			score+=1;
		#else:
		#	print("mismatch at %i"%(i));
	return score/(len(seq1)-offset)

## the function to return the alignment location of the sequences
def align(seq1, seq2, offset, overlapRange):
	for i in range(0,overlapRange+1):
		if testAlignment(seq1,seq2, offset + i) > 0.75:
			return (offset + i)
		elif testAlignment(seq1,seq2, offset - i) > 0.75:
			return (offset - i)

## the function to output the alignment sequence
def getConsensus(read1, read2, offset):
	consensus  = [];
	#print("from %i to %i"%(offset,len(read1[1])));
	for i in range(offset, len(read1[1])):
		if read1[2][i]> read2[2][i-offset]: #converts to ASCII, corresponds to quality; if equal, doesn't matter which is better since base is the same
			consensus.append(read1[1][i]);
		else:
			consensus.append(read2[1][i-offset]);
	#print("read 2 from %i"%((len(read1[1])-offset)));
	return read1[1][0:(offset)] + "".join(consensus) + read2[1][(len(read1[1])-offset):]

args.overlap = int(args.overlap);
args.overlapRange = int(args.overlapRange);

#tests
#read1 = ("@M01581:878:000000000-AP1UB:1:1101:13431:1072 1:N:0:1", "TGCATTTTTTTCACATCTACTGAACGGTGACGGTTGCGGCCCACGTGGCCCCCGGGAAAGCGGAGGGCGGCTGCT", "CCCCCGGGGGGGGGGGGGFEGFGGGFGFCCFEGGGGGGGCC+CBFG8FEFFEEFEG@F<FFE7=FEDFGDECFEF");
#read2 = ("@M01581:878:000000000-AP1UB:1:1101:13431:1072 2:N:0:1", "AACAGCCGTAACCTCGACCTCGTCCACGCAGGTCAAGCAGCCGCCCTCCGCTTTCCCGGGGGCCACGTGGGCCGC", "CCCCCGGGGGGGGFGDGEG@F;FFGGEFGGGGEFGFFGF@,CFEEFCFCC@CFFGFFC@7C+@+8CFG8C7+B>F");
#print("seq1 = "+read1[1]);
#print("seq2 = "+read2[1]);
#print("rc(seq2) = "+revcomp(read2[1]));
#print("qual2 = "+read2[2]);
#read2 = (read2[0], revcomp(read2[1]), read2[2][::-1]);
#print("rev(qual2) = "+read2[2]);
#print(str(len(read1[1]) - 40));
#print("Alignment = "+str(testAlignment(read1[1], read2[1],len(read1[1]) - args.overlap)));
#print("Alignment offset = "+str(align(read1[1], read2[1],len(read1[1]) - args.overlap, args.overlapRange )));
#print("Consensus = " + getConsensus(read1, read2, align(read1[1], read2[1],len(read1[1]) - args.overlap, args.overlapRange )));
#read1 = ("@M01581:878:000000000-AP1UB:1:1101:13431:1072 1:N:0:1", "TGCATTTTTTTCACATCTACTGAACGGTGACGGTTGCGGCCCACGTGGCCCCCGGGAAAGCGGAGGGCGGCTGCG", "CCCCCGGGGGGGGGGGGGFEGFGGGFGFCCFEGGGGGGGCC+CBFG8FEFFEEFEG@F<FFE7=FEDFGDECFEJ");
#print("Consensus = " + getConsensus(read1, read2, align(read1[1], read2[1],len(read1[1]) - args.overlap, args.overlapRange )));
#read1 = ("@M01581:878:000000000-AP1UB:1:1101:13431:1072 1:N:0:1", "TGCATTTTTTTCACATCTACTGAACGGTGACGGTTGCGGCCCACGTGGCCCCCGGGAAAGCGGAGGGCGGCTGCG", "CCCCCGGGGGGGGGGGGGFEGFGGGFGFCCFEGGGGGGGCC+CBFG8FEFFEEFEG@F<FFE7=FEDFGDECFE!");
#print("Consensus = " + getConsensus(read1, read2, align(read1[1], read2[1],len(read1[1]) - args.overlap, args.overlapRange )));
#read1 = ("@NB501583:183:HYK5YBGX2:1:11101:12773:1776 1:N:0:GGAGCTAC+CTAGTCGA", "CTAATCAAGTGCTAGCTAGACATCTATATATAACAAGCACAGAACCGTCTAATTGGTATTTTTCAGGACATTTTAAACATCCGTACAACGAGAACCCATACATTACTTTTTTTAATATTCTCGACCACATAGCCCCCAATGCATGTCCACACAGAAAGTGTTGCGCACAACCCCGGTAACCATCCCTATG", "AAAAAEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEE<AEEEAAEEEEEEEEEEEEEEEEEEE<AAEEEEAEEEEEE<EEEEAEEEEAEEEEEEEEAEEAEEE<EE<AAAAEE<AAE<<<AAAAEAA66<666AA");
#read2 = ("@NB501583:183:HYK5YBGX2:1:11101:12773:1776 2:N:0:GGAGCTAC+CTAGTCGA", "CTGCTCAGGATCCTAAACGTTATAATACACGAAATGAATGCCTCTCTTAACAGTATACTTCTTTAATGATGAGTTATGTCCAACCTAGGGGGACTACTGGTGCCCATAGGGA", "AAAAAEEEEEEEEEEEEEEEEEEAEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEA<");
#read2 = (read2[0], revcomp(read2[1]), read2[2][::-1]);
#print("Alignment = "+str(testAlignment(read1[1], read2[1],len(read1[1]) - args.overlap)));
#print("Alignment offset = "+str(align(read1[1], read2[1],len(read1[1]) - args.overlap, args.overlapRange )));
#print("Best alignment = "+str(testAlignment(read1[1], read2[1],align(read1[1], read2[1],len(read1[1]) - args.overlap, args.overlapRange ))));
#print("Consensus = " + getConsensus(read1, read2, align(read1[1], read2[1],len(read1[1]) - args.overlap, args.overlapRange )));

#exit();

## determinr the alignment is done or not
notDone = True;
i=0;
while notDone:
	i=i+1;
	d1 = getNextRead(inFile1);
	d2 = getNextRead(inFile2);
	if d1[0]=="" or d2[0]=="":
		if d1[0]=="" and d2[0]=="":
			notDone=False;
		else:
			raise Exception("One of the files ended prematurely at fastq entry %i (line ~%i)"%(i,i*4));
	else:
		## I dont know what is the revcomp command is, dearch in google is the reverse complement
		d2 = (d2[0], revcomp(d2[1]), d2[2][::-1]);
		offset = align(d1[1],d2[1], len(d1[1]) - args.overlap, args.overlapRange);
		if args.verbose>1:
			sys.stderr.write("Alignment offset = %i\n"%(offset));
		if offset is not None:
			consensus = getConsensus(d1, d2, offset);
			outFile.write("%s\t%s\t%i\n"%(d1[0],consensus,len(d1[1])-offset));
		else:
			sys.stderr.write("Skipping %s/%s because they don't align within given parameters\n"%(d1[0],d2[0]));

inFile1.close();
inFile2.close();
outFile.close();
if (args.logFP is not None):
	logFile.close();
## this is the same command as before but after the if the command changed, I don't know why
