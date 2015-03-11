#! /usr/bin/python
import fileinput, string, sys

# This script takes a gnuplot data file and extracts data of only a specified index of each block.
# The result is written to the standard output and can be used as an input to gnuplot.

try:
	if(len(sys.argv) < 3) :
		print "Usage: python extract.py <filename> index=<node index> for extracting indices"
		print "or     python extract.py <filename> x=<position> for extracting positions along y-axis"
		print "or     python extract.py <filename> y=<position> for extracting positions along x-axis"
		print "or     python extract.py <filename> x=<position> y=<position> for extracting single positions"

	filename = sys.argv[1]
	index = 0
	useXPosition = False
	useYPosition = False
	
	for i in range(2,len(sys.argv)) :
		#print "i=", i
		str = sys.argv[i]
		#print str

		arg2 = str.split("=")
		if arg2[0] == "index" :
			index = arg2[1]
			#print "Searching for index ", index
		elif arg2[0] == "x" :
			useXPosition = True
			xPosition = float(arg2[1])
			#print "Searching for x position ", xPosition
		elif arg2[0] == "y" :
			useYPosition = True
			yPosition = float(arg2[1])
			#print "Searching for y position ", yPosition
		else :
			print "Usage: python extract.py <filename> index=<node index> for extracting indices"
			print "or     python extract.py <filename> x=<position> for extracting positions along y-axis"
			print "or     python extract.py <filename> y=<position> for extracting positions along x-axis"
			print "or     python extract.py <filename> x=<position> y=<position> for extracting single positions"
			exit

	#print "useXPosition: ", useXPosition
	#print "useYPosition: ", useYPosition

	cols = ""
	time = 0.0

	nTimeSteps = 0
	indexCounter = 0
	currentXPosition = 0.0
	currentYPosition = 0.0
	lastXPosition = 0.0
	lastYPosition = 0.0
	
	foundLastY = False

	f = open(filename, 'r')
	for line in f :
		if '#' in line :
			indexCounter = 0
			nTimeSteps += 1
			splitted = line.split()
			if(len(splitted) > 1) :
				time = [2]
			if(not(useXPosition and useYPosition)) :
				if(nTimeSteps != 1) :
					print "\n"
				print line,
		elif ((not useXPosition) and (not useYPosition)) and (line != "\n") :
			if indexCounter == int(index) :
				cols = [time]
				cols.extend(line.split())
				for col in cols :
					print col,
				print ""
			indexCounter += 1
		elif (useXPosition or useYPosition) and (line != "\n") :
			# Initialize flags for successful search
			if (useXPosition) :
				foundXPosition = False
			else :
				foundXPosition = True
			if (useYPosition) :
				foundYPosition = False
			else:
				foundYPosition = True
			
			lastXPosition = currentXPosition
			currentXPosition = float(line.split()[0])

			temp = lastYPosition
			lastYPosition = currentYPosition
			currentYPosition = float(line.split()[1])
			# same y values in this and previous line => restore y value of the preceding block
			if(currentYPosition == lastYPosition) :
				lastYPosition = temp
			if useXPosition :
				# Reset x position when entering a new column
				if(lastXPosition > currentXPosition) :
					lastXPosition = 0.0;
				#print lastXPosition, " | ", xPosition, " | ", currentXPosition
				if (currentXPosition == xPosition) or (lastXPosition < xPosition and currentXPosition > xPosition):
					foundXPosition = True
			if useYPosition :
				if (currentYPosition == yPosition) or (lastYPosition < yPosition and currentYPosition > yPosition):
					foundYPosition = True
					#print "lastYPosition = ", lastYPosition	
					#print "currentYPosition = ", currentYPosition
					if (currentXPosition == 0 and foundLastY and not(useXPosition and useYPosition)) :
						print ""
					foundLastY = True
				else :
					foundLastY = False
			# Print lines containing at desired positions
			if(foundXPosition and foundYPosition) :
				#cols = [time]
				#cols.extend(line.split())
				cols = line.split()
				k = -1
				# Print time instead of position
				if (useXPosition and useYPosition) :
					print time,
				for col in cols :
					k = k+1
					if (useXPosition and k == 0) :
						continue
					if (useYPosition and k == 1) :
						continue		
					print col,
				print ""			

			indexCounter += 1
		else :
			if line == "\n" :
				#print "empty line"
				indexCounter -= 1
			else :
				print "Error!"
				exit(1)
				
	# Close input file
	f.close()
except IOError as e:
	# Disbaled because of annoying error messages when calling this script from gnuplot
	#print "I/O error({0}): {1}".format(e.errno, e.strerror)
	f.close()
