"""The purpose of this is to organize all file reading/parsing and writing information.
This file is imported and called by all other files."""

## So files can be opened/read/saved as .gz files
import gzip

## Function to open and write files
def smartGZOpen(filename,mode):
    	"""This function writes/modifies and opens files.
	Input: A specific file and a 'mode' or state in
	which it should be opened.
	Output: Returns the file, opened."""
	## Modify file name to be accessible by gzip
    	## Incase there are too many/it is too long
	if len(filename)>3 and filename[-3:].lower()=='.gz':
		return gzip.open(filename,'%st'%(mode));
	else:
        	## Open file as it is
		return open(filename,mode);

## Function to read files
def smartGZForeach(filename):
    	"""This function parses through and reads files.
	Input: A specific file
	Output: Reads through each line of the file."""
	## Read each line in file
	inFile = smartGZOpen(filename,"r")
	for line in inFile:
        ## Return the line itself
		yield line
