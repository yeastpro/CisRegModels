"""
This file consists of various functions with the goal of reading and
writing items as part of a batch.
Used for converting sequences into binary (OHC) and into vectors for tensorflow
"""
## Code annotated by K
from threading import Thread;
from . import MYUTILS
import numpy as np;
import tensorflow as tf;
import sys
from tensorflow.python.framework import ops
from tensorflow.python.ops import math_ops

## This function determines the learning rate of the model
## Uses tensorflow (machine learning tool)
def cycle_learning_rate(learning_rate_low, learning_rate_high, global_step, period, name=None):
	"""This function determines the learning rate of the model,
	uses tensorflow.
	Input: High and low rate values, global variable for rate step, period
	Output: Rate of model (relative to all values)"""
	if global_step is None:
		raise ValueError("global_step is required for cycle_learning_rate.")
	with ops.name_scope(name, "CycleLR",
											[learning_rate_high, learning_rate_low, global_step,
											 period]) as name:
		## Converting np arrays to tensor
		learning_rate_high = ops.convert_to_tensor(learning_rate_high, name="learning_rate_high")
		learning_rate_low = ops.convert_to_tensor(learning_rate_low, name="learning_rate_low")
		## Initialize (gloabl) variables and types
		dtype = learning_rate_high.dtype
		global_step = math_ops.cast(global_step, dtype)
		period = math_ops.cast(period, dtype)
		step_size = tf.divide(tf.subtract(learning_rate_high,learning_rate_low),period) # increment by this until period is reached, then decrement
		cycle_stage = tf.abs(tf.subtract(tf.floormod(global_step,tf.multiply(math_ops.cast(2.0,dtype),period)),period))
		return tf.add(learning_rate_low, tf.multiply(tf.subtract(period,cycle_stage),step_size), name=name)

#abstract object for polymorphism
## Define class and function for getting specific batches
class BasicBatchGetter:
	def getNextBatch(self):
		"""This function in the class is used to get specific
		batches and define values.
		Input: Self variable
		Output: Class values"""
		self.curThread.join();
		curX = self.nextBatchX;
		curY = self.nextBatchY;
		curR = self.numRuns;
		self.curThread = Thread(target = self.prepareNextBatch);
		self.curThread.start()
		return(curX, curY, curR);

#asynchronous batch getting for one hot encoding of sequences
## Define class and batch for getting encoded sequences
## OHC format (binary)
## No specific return - value storage mechanism
## Class calles basicbatchgetter class (from above)
class BatchGetterOneHot(BasicBatchGetter):
	def __init__(self, inFP, batchSize, numRuns, seqLen):
		"""This function gets encoded sequences of batches.
		Input: Self, batch details, sequence length
		Output: Initialized self values"""
		## Define sequence parameters
		self.inFP = inFP;
		self.batchSize = batchSize;
		self.numRuns= numRuns;
		self.seqLen= seqLen;
		self.curFH = MYUTILS.smartGZOpen(self.inFP,'r')
		self.curThread = Thread(target = self.prepareNextBatch);
		self.curThread.start()

	## function for prepare another batch based on a previous one
	def prepareNextBatch(self):
		"""This function prepares a batch (based on previous).
		Input: Self variable
		Output: Stores self values and updated variables"""
		self.nextBatchX = np.zeros((self.batchSize,4,self.seqLen,1))
		self.nextBatchY = np.zeros((self.batchSize))
			## Initialize counter
		b=0
		while b < self.batchSize:
			line = self.curFH.readline()
			if line =="":
				if self.numRuns==1:
								## Update variables
					self.nextBatchX = self.nextBatchX[0:b,:,:,:]
					self.nextBatchY = self.nextBatchY[0:b]
					self.numRuns-=1;
					return;
				self.curFH.close();
						## Read file
				self.curFH = MYUTILS.smartGZOpen(self.inFP,'r')
				self.numRuns-=1;
				line = self.curFH.readline()
					## Open files
			if line is None or line[0]=="#": continue
			curData = np.fromstring(line, dtype=float, sep="\t")
			self.nextBatchY[b]=curData[0];
			self.nextBatchX[b,:,:,0] = curData[1:].reshape((4,self.seqLen))
					## Update counter
			b+=1


## Function using the basicbatchgetter class
## Optimizes original function
## Same as batchgetteronehot but not encoded
class BatchGetter(BasicBatchGetter):
	def __init__(self, inFP, batchSize, numRuns,numTFs, numKds):
		"""This function optimizes the original batch getter function.
		Same as the previous function but unencoded.
		Input: Self variable, batch values
		Output: Stores self values and updated variables"""
		## Defining parameters
		self.inFP = inFP;
		self.batchSize = batchSize;
		self.numRuns= numRuns;
		self.numTFs= numTFs;
		self.numKds= numKds;
		self.curFH = MYUTILS.smartGZOpen(self.inFP,'r')
		self.curThread = Thread(target = self.prepareNextBatch);
		self.curThread.start()

	def prepareNextBatch(self):
		"""This functionprepares a batch (based on previous, unencoded).
		Input: Self variable
		Output: Stores self values and updated variables"""
		self.nextBatchX = np.zeros((self.batchSize,self.numKds,self.numTFs)) +np.log(99999.0);
		self.nextBatchY = np.zeros((self.batchSize))
		b=0
		while b < self.batchSize:
			line = self.curFH.readline()
			if line =="":
				if self.numRuns==1:
					self.nextBatchX = self.nextBatchX[0:b,:,:]
					self.nextBatchY = self.nextBatchY[0:b]
					self.numRuns-=1;
					return;
				self.curFH.close();
				self.curFH = MYUTILS.smartGZOpen(self.inFP,'r')
				self.numRuns-=1;
				line = self.curFH.readline()
			if line is None or line[0]=="#": continue
			curData = line.split("\t");
			self.nextBatchY[b]=float(curData[0]);
			for t in range(1,len(curData)):
				curKds = [np.log(float(x)) for x in curData[t].split(";")]
				self.nextBatchX[b,0:min(self.numKds,len(curKds)),t-1] = curKds[0:min(self.numKds,len(curKds))];
			b+=1

## New batchgetter class
## Updated numKds values
## Otherwise same code as before
class BatchGetterFixedNumKds(BatchGetter):
	def prepareNextBatch(self):
		"""This function updates numKds values.
		Input: Self variable
		Output: Stores updated value"""
		self.nextBatchX = np.zeros((self.batchSize,self.numKds,self.numTFs)) +np.log(99999.0);
		self.nextBatchY = np.zeros((self.batchSize))
		b=0
		while b < self.batchSize:
			line = self.curFH.readline()
			if line =="":
				if self.numRuns==1:
					self.nextBatchX = self.nextBatchX[0:b,:,:]
					self.nextBatchY = self.nextBatchY[0:b]
					self.numRuns-=1;
					return;
				self.curFH.close();
				self.curFH = MYUTILS.smartGZOpen(self.inFP,'r')
				self.numRuns-=1;
				line = self.curFH.readline()
			if line is None or line[0]=="#": continue
			curData = np.fromstring(line, dtype=float, sep="\t")
			self.nextBatchY[b]=curData[0];
			self.nextBatchX[b,:,:] = np.transpose(curData[1:len(curData)].reshape((self.numTFs,self.numKds)))
			## Update counter
			b+=1

## Function with similar code as above
## Updated class, calling basicbatchgetter
## Used to conver sequences to vectors for tensorflow
class BatchGetterSeq2Vec(BasicBatchGetter):
	def __init__(self, inFP, batchSize, numRuns,seqLen, kmer2index, wordLen):
		"""This function converts sequences to vectors using tensor flow.
		Input: Self variable, batch information
		Output: Stores self variable information"""
		self.inFP = inFP;
		self.batchSize = batchSize;
		self.numRuns= numRuns;
		self.seqLen= seqLen;
		self.wordLen= wordLen;
		self.kmer2index= kmer2index;
		self.curFH = MYUTILS.smartGZOpen(self.inFP,'r')
		self.curThread = Thread(target = self.prepareNextBatch);
		self.curThread.start()

	def prepareNextBatch(self):
		"""This function uses the __init__ function values to
		actively convert sequences to vectors.
		Input: Self variable
		Output: Stores newly converted vector information"""
		self.nextBatchX = np.zeros((self.batchSize,self.seqLen-self.wordLen+1)).astype("int32");
		self.nextBatchY = np.zeros((self.batchSize))
		b=0
		while b < self.batchSize:
			line = self.curFH.readline()
			if line =="":
				if self.numRuns==1:
					self.nextBatchX = self.nextBatchX[0:b,:,:]
					self.nextBatchY = self.nextBatchY[0:b]
					self.numRuns-=1;
					return;
				self.curFH.close();
				self.curFH = MYUTILS.smartGZOpen(self.inFP,'r')
				self.numRuns-=1;
				line = self.curFH.readline()
			if line is None or line[0]=="#": continue
			curData = line.rstrip().split("\t");
			self.nextBatchY[b]=float(curData[0]);
			curSeq = curData[1];
			if len(curSeq) < self.seqLen:
				curSeq = "N"*(self.seqLen - len(curSeq)) + curSeq;  ### prepend Ns if the sequence is too short
			curSeq = curSeq[(len(curSeq)-self.seqLen):len(curSeq)]# trim distal bases if too long
			for si in range(0,self.seqLen-self.wordLen+1):
				self.nextBatchX[b,si] = self.kmer2index[curSeq[si:(si+self.wordLen)]] #fill X with the indeces of the various k-mers
			b+=1
