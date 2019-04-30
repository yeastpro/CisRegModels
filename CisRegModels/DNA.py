"""The purpose of this file is to initialize a translation table (of all variations of A,T,C,G)
for sequences to be compared to, and to calculate the levenshtein distances
between specified sequences"""

#from string import maketrans
#import string
## Build translation table of specified characters (like a dictionary)
rcTab = str.maketrans("ATGCatgc","TACGtacg");

## Build function to tie the specified input sequence to the translation table
## (initialized above)
def revcomp(seq):
	"""This function ties the specified input sequences to the
	translational table (initialized above).
	Input: a specified sequences
	Output: String where characters are mapped
	to a translational table"""
    ## Make the translation table a global variable
	global rcTab
    ## Returns string where each character in the input is mapped
    ## to a character in the translation table
    ## Takes everything in the dimension, but backwards ([::-1])
	return seq.translate(rcTab)[::-1];


## Build function to calculate levenshtein distance
## This is a string metric for measuring the difference between two sequences
## Returns the distance between specified sequences
def levenshtein(s1, s2):
	"""This function calculates teh levenshtein distance between two sequences.
	Input: Two specified sequences
	Output: Distance between the sequences"""
	if len(s1) < len(s2):
        ## Calls function with inputs swapped
		return levenshtein(s2, s1)
	# len(s1) >= len(s2)
	if len(s2) == 0:
		return len(s1)
	## Initialize 'counter' of sorts to be used in loop
    	## Only of s2 becasue compared in a nested for loop with s2
	previous_row = range(len(s2) + 1)
	## Initilaize for loop through s1
	for i, c1 in enumerate(s1):
		current_row = [i + 1]
        	## Initilaize nested for loop through s2
		for j, c2 in enumerate(s2):
            		# use j+1 instead of j since previous_row and current_row are one character longer
            		# than s2
            		## Defining different types of sequence modifications
			insertions = previous_row[j + 1] + 1 ## single/multiple base inserted into the sequence
			deletions = current_row[j] + 1 ## single/multiple base deleted from sequence
			substitutions = previous_row[j] + (c1 != c2) ## single or multiple bases substituted in for other bases
			current_row.append(min(insertions, deletions, substitutions)) ## append number of each per row
        	## Update 'counter'
		previous_row = current_row
	return previous_row[-1]

