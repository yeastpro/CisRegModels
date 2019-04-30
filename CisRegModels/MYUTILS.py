## File containing file reading/writing information
## Imported and called by all other files

## So files can be opened/read/saved as .gz files
import gzip

## Function to open and write files
def smartGZOpen(filename,mode):
    ## Modify file name to be accessible by gzip
    ## Incase there are too many/it is too long
	if len(filename)>3 and filename[-3:].lower()=='.gz':
		return gzip.open(filename,'%st'%(mode));
	else:
        ## Open file as it is
		return open(filename,mode);

## Function to read files
def smartGZForeach(filename):
    ## Read each line in file
	inFile = smartGZOpen(filename,"r")
	for line in inFile:
        ## Return the line itself
		yield line
