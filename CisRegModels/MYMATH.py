"""
The purpose of this function is to store all important and repeatedly used
mathematical functions in a single file that is accessible by all other files.
"""
## Code annotated by K
import math
import numpy
from . import MYUTILS

## Function to initialize a matrix from specified file
## Returns column names, row names and the matrix itself
def inputMatrix(inFile, inType = numpy.float, colnames=True, rownames=True):
    	"""This function initializes a matrix from specified files
	and returns the matrix and additional values.
	Input: File opening function, file type columns and rows
	Output: Matrix, with columns and rows"""
	## Initialize list
	colNames = []
	## Initialize counters
	nRows = 0;
	nCols = 0;
	maxString=1;
	first=colnames;
    	## Parse through file and store only relevant characters
	for line in inFile:
		if line is None or line=="": continue;
        	## Data after trailing characters are stripped
		data = line.rstrip().split("\t");
        	## Specifying column information
		if first:
			first=False;
			colNames = data;
		else:
			nCols = len(data);
			if inType==numpy.str:
				for i in range(1,len(data)):
					maxString = max(maxString,len(data[i]))
            		## Update counter
			nRows+=1;
	## Specifying row information (separate of for loop)
	if rownames:
		nCols = nCols-1;
		rowLabs = [""] * nRows
	else:
		rowLabs = [];
	if inType==numpy.str:
		dataMatrix = numpy.empty((nRows,nCols), dtype = numpy.dtype('|S%i'%(maxString)));
	else:
		dataMatrix = numpy.empty((nRows,nCols), dtype = inType)
    	## Update counter
	nRows = 0;
	first=colnames;
	inFile.seek(0,0);
    	## Same as for loop above, but for specifying row information
	for line in inFile:
		if line is None or line=="": continue;
		data = line.rstrip().split("\t");
		if first:
			first=False;
		else:
			if rownames:
				rowLabs[nRows] = data[0];
				dataMatrix[nRows,:] = numpy.array([data[1:len(data)]]).astype(inType)
			else:
				dataMatrix[nRows,:] = numpy.array([data[0:len(data)]]).astype(inType)
			nRows+=1;
	return (rowLabs, colNames, dataMatrix);

## Function to write and save newly formed matrix
def saveMatrix(outFileName, rowLabs, colLabs, dataMatrix):
	"""This function writes and saves a matrix.
	Input: File writing function, generated matrix, row and column values
	Output: Writes and saves matrix"""
	outFile = MYUTILS.smartGZOpen(outFileName, "w");
	outFile.write("\t".join(colLabs)+"\n");
	for i in range(0,len(rowLabs)):
		outFile.write(rowLabs[i]);
		for j in range(0,dataMatrix.shape[1]):
			outFile.write("\t%g"%dataMatrix[i,j]);
		outFile.write("\n");
	outFile.close();

## Function to scale the matrix/data to a specified length/size
def scaleLength(data,desiredLength):
	"""This function scales a matrix.
	Input: Matrix and desired length
	Output: Scaled matrix"""
	scaleMiddle(data,0,length(data),desiredLength);

## Function to scale the data to a mean/a desired length
## Returns a scaled matrix
def scaleMiddle(data,fStart,fEnd,desiredLength,func="MEAN"): #from start to end-1
	"""This function scales a matrix to a mean or desired length.
	Input: Matrix, desired length and range, calles python's mean function
	Output: Scaled matrix"""
	#now scale middle bin accordingly to desiredLength.
	thisMid = reshapeData(data, fStart, fEnd, desiredLength, func);
	#print(data[:fStart]);
	#print(thisMid)
	scaledData = list(data[:fStart]) + thisMid + list(data[fEnd+1:]);
	return scaledData;

## Make sure all components are floats
def isfloat(x):
	"""This function converts all components to a float.
	Input: Data (matrix, list etc.)
	Output: Data (all components converted to floats)"""
	return issubclass(type(x), numpy.float) or isinstance(x,float)

## Reshapes data to make it easier to use
## Combines data using average
def reshapeData(curProfile,fStart,fEnd,desiredLength,func):
	"""This function reshapes data to make it more usable.
	Input: Profile to compare to/data to be scaled, desired length and range, matrix function
	Output: Combines input data using an average, outputs scaled data"""
	#combines data by average
	#print(fStart)
	#print(fEnd)
	#print(desiredLength)
	curLen = fEnd-fStart+1;
	#print(curLen)
    	## Initialize values
	thisMid = [0.0]*desiredLength; ## Float
	theCounts = [0]*desiredLength; ## Int
	mapRatio=(0.0+desiredLength)/curLen;
    	#print(mapRatio)
    	## If input length is greater than desired length
	if curLen>=desiredLength:#compress
		for j in range(0,curLen):
			if ((fStart+j)>=len(curProfile) or math.isnan(curProfile[fStart+j]) or not numpy.isreal(curProfile[fStart+j])):
				continue  #skip if data has a NaN
			dest=(mapRatio*j);
			ceilPct =0;
			floor=int(math.floor(dest));
			ceil=int(math.floor(dest+mapRatio));
			if ceil>floor:
				ceilPct=(dest+mapRatio- math.floor(dest+mapRatio))/mapRatio;
			floorPct=1-ceilPct;
			#print([curLen, desiredLength, dest, mapRatio, floorPct, ceilPct, floor, ceil].join("\t")+"\n");
			#p [fStart, j, floorPct, floor];
			#p curProfile;
        		## Update counters
			thisMid[floor]+=floorPct*curProfile[fStart+j];
			theCounts[floor]+=floorPct;
			if ceil<desiredLength:
				thisMid[ceil]+=ceilPct*curProfile[fStart+j];
				theCounts[ceil]+=ceilPct;
    	## If input length is less than desired lengtbh
	else: #expand
		for j in range(0,curLen):
			if ((fStart+j)>=len(curProfile) or  math.isnan(curProfile[fStart+j]) or not numpy.isreal(curProfile[fStart+j])):
				continue;
			dest = (mapRatio*(j));
			nextDest = dest+mapRatio;
			first=int(math.floor(dest))
			last=int(math.floor(nextDest))
			totalPct=0.0;
			for l in range(first,last+1):
				if l==desiredLength:
					continue;
				elif l==first:
					curPct = (1-(dest-math.floor(dest)))/mapRatio;
				elif l==last:
					curPct = (nextDest-math.floor(nextDest))/mapRatio;
				else:
					curPct=1.0/mapRatio;
				#print([curLen, desiredLength, dest, nextDest, mapRatio, first, last, l, curPct].join("\t")+"\n");
				theCounts[l]+=curPct;
				thisMid[l]+=curPct*curProfile[fStart+j];
				totalPct+=curPct;
        	    		## If not in range, throw exception - error
			if abs(totalPct-1.0)>=0.001:
				print("ERROR: Total Percent is not 1: "+totalPct.to_s()+"\n");
    	# Replace those without data withno data with NaN and average the rest
    	## Use values from above for loops to return thisMid value
    	## thisMid gives reshaped data value
	for i in range(0, len(thisMid)):
		if theCounts[i]==0:
			thisMid[i]=float('NaN');
		elif func=="MEAN":
			thisMid[i] = thisMid[i]/theCounts[i];
	return thisMid;

## Function returns angle (in 2D) between sequences
## Returns d{theta}
def Angle2D(x1, y1, x2, y2):
    """This function takes the end points on vectors and computes
    the angle between them.
    Input: Endpoints of vectors
    Output: d{theta} angle between the vectors"""
    ## Solve for theta values
    theta1 = np.arctan2(y1,x1);
    theta2 = np.arctan2(y2,x2);
    dtheta = theta2 - theta1;
    ## if dtheta is greater than pi, subtract
    while (dtheta > np.pi):
        dtheta -= 2*np.pi;
    ## if dtheta is less than pi, add
    while (dtheta < -np.pi):
        dtheta += 2*np.pi;
    return(dtheta);

## This function determines whether or not the specified point values
## are inside the specified polygon
def isInPolygon(polygonX, polygonY, pointsX, pointsY):
    """This function determine whether specific points are within
    a specified (polygonal) area
    Input: Area boundaries, specified points
    Output: Confirmation of all the specified points that exist within the boundaries"""
    ## Initialize list
    allInside = [];
    ## Initialize for loop to test pointsX
    for i in range(0,len(pointsX)):
        angle=0;
        ## Initialize nester for loop to test pointsY
        for j in range(0,len(polygonX)):
            angle += Angle2D(polygonX[j]-pointsX[i], polygonY[j]-pointsY[i],polygonX[(j+1)%len(polygonX)]-pointsX[i],
polygonY[(j+1)%len(polygonX)]-pointsY[i]);
        ## Append values to initialized list if the angle is greater than or equal to pi
        ## This means that the point lies inside the polygon
        if (np.abs(angle) >= np.pi):
            allInside.append(i)
    return allInside;
#
#def reshapeDataAvg(curProfile,fStart,fEnd,desiredLength):
#	#combines data by average
#	curLen = fEnd-fStart
#	thisMid = [0.0]*desiredLength;
#	theCounts = [0.0]*desiredLength;
#	mapRatio=(0.0+desiredLength)/curLen;
#	if curLen>=desiredLength:#compress
#		for j in range(0,curLen):
#			dest=(mapRatio*(j));
#			ceilPct =0;
#			floor=dest.floor();
#			ceil=(dest+mapRatio).floor();
#			if ceil>floor:
#				ceilPct=(dest+mapRatio- (dest+mapRatio).floor())/mapRatio;
#			floorPct=1-ceilPct;
#			#print([curLen, desiredLength, dest, mapRatio, floorPct, ceilPct, floor, ceil].join("\t")+"\n");
#			#p [cPIndex, j, floorPct, floor];
#			#p curProfile;
#			if curProfile[cPIndex+j]==NaN  or  curProfile[cPIndex+j] is None:
#				continue  #skip if data has a NaN
#
#			thisMid[floor]+=floorPct*curProfile[cPIndex+j];
#			theCounts[floor]+=floorPct;
#			if ceil<desiredLength:
#				thisMid[ceil]+=ceilPct*curProfile[cPIndex+j];
#				theCounts[ceil]+=ceilPct;
#	else: #expand
#		for j in range(0,curLen):
#			dest = (mapRatio*(j));
#			nextDest = dest+mapRatio;
#			first=dest.floor();
#			last=(nextDest).floor()
#			totalPct=0.0;
#			first.upto(last){|l|
#				if l==desiredLength:
#					continue;
#				elif l==first:
#					curPct = (1-(dest-dest.floor()))/mapRatio;
#				elif l==last:
#					curPct = (nextDest-nextDest.floor())/mapRatio;
#				else:
#					curPct=1.0/mapRatio;
#				#print([curLen, desiredLength, dest, nextDest, mapRatio, first, last, l, curPct].join("\t")+"\n");
#				if (curProfile[cPIndex+j] is None  or  curProfile[cPIndex+j]==NaN):
#					continue;
#				theCounts[l]+=curPct;
#				thisMid[l]+=curPct*curProfile[cPIndex+j];
#				totalPct+=curPct;
#			if totalPct-1.0>=0.001:
#				print("ERROR: Total Percent is not 1: "+totalPct.to_s()+"\n");
#	for i in range(0, thisMid.length()): # Replace those without data withno data with NaN and average the rest
#		if theCounts[i]==0:
#			thisMid[i]=NaN;
#		else:
#			thisMid[i] = thisMid[i]/theCounts[i];
#	return thisMid;
#
#def reshapeDataSum(curProfile,curLen, cPIndex):
#	#combines data by sum -unsure if this is actually working correctly because it's a weird method - I think it mght also average
#	curLen = fEnd-fStart
#	thisMid = [0.0]*desiredLength;
#	mapRatio=(0.0+desiredLength)/curLen;
#	if curLen>=desiredLength: #compress
#		for j in range(0,curLen):
#			dest=(mapRatio*(j));
#			ceilPct =0;
#			floor=dest.floor();
#			ceil=(dest+mapRatio).floor();
#			if ceil>floor:
#				ceilPct=(dest+mapRatio- (dest+mapRatio).floor())/mapRatio;
#			floorPct=1-ceilPct;
#			#print([curLen, desiredLength, dest, mapRatio, floorPct, ceilPct, floor, ceil].join("\t")+"\n");
#			#p [cPIndex, j, floorPct, floor];
#			#p curProfile;
#			if curProfile[cPIndex+j]==NaN: #skip if data has a NaN
#				continue
#			#p([curProfile.length(),j,cPIndex, curLen]);
#			thisMid[floor]+=floorPct*curProfile[cPIndex+j];
#			if ceil<desiredLength:
#				thisMid[ceil]+=ceilPct*curProfile[cPIndex+j];
#	else: #expand
#		for j in range(0,curLen):
#			dest = (mapRatio*(j));
#			nextDest = dest+mapRatio;
#			first=dest.floor();
#			last=(nextDest).floor()
#			totalPct=0.0;
#			first.upto(last){|l|
#				if l==desiredLength:
#					continue;
#				elif l==first:
#					curPct = (1-(dest-dest.floor()))/mapRatio;
#				elif l==last:
#					curPct = (nextDest-nextDest.floor())/mapRatio;
#				else:
#					curPct=1.0/mapRatio;
#				#print([curLen, desiredLength, dest, nextDest, mapRatio, first, last, l, curPct].join("\t")+"\n");
#				if curProfile[cPIndex+j]!=NaN:
#					thisMid[l]+=curPct*curProfile[cPIndex+j];
#				totalPct+=curPct;
#			if totalPct-1.0>=0.001:
#				print("ERROR: Total Percent is not 1: "+totalPct.to_s()+"\n");
#	return thisMid;
