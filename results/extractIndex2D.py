#! /usr/bin/python
import fileinput, string, sys

# This script takes a gnuplot data file and extracts data of only a specified index of each block.
# The result is written to the standard output and can be used as an input to gnuplot.

if(len(sys.argv) < 3) :
	print "Usage: ./extractIndex <filename> index=<node index> for extracting indices"
	print "or     ./extractIndex <filename> pos=<x>,<y> for extracting positions"

filename = sys.argv[1]
str = sys.argv[2]

index = 0
position = [0.0, 0.0]
usePosition = True

arg2 = str.split("=")
if arg2[0] == "index" :
	usePosition = False
	index = arg2[1]
	print "Searching for index ", index
elif arg2[0] == "pos" :
	usePosition = True
	posString = arg2[1].split(",")
	print posString
	position[0] = float(posString[0])
	position[1] = float(posString[1])
	print "Searching for position ", position
else :
	print "Usage: ./extractIndex <filename> index=<node index> for extracting indices"
	print "or     ./extractIndex <filename> pos=<node index> for extracting positions"
	exit


cols = ""
time = 0.0

nTimeSteps = 0
indexCounter = 0
currentPosition = [0.0, 0.0]
lastPosition = [0.0, 0.0]

referencePosition = 50.0

f = open(filename, 'r')
for line in f :
	if '#' in line :
		indexCounter = 0
		nTimeSteps += 1
		time = line.split()[2]
		#print "[", nTimeSteps, "]", time
	elif (not usePosition) and (line != "\n") :
		if indexCounter == int(index) :
			cols = [time]
			cols.extend(line.split())
			for col in cols :
				print col,
			print ""
		indexCounter += 1
	elif (usePosition) and (line != "\n") :
		lastPosition = currentPosition
		currentPosition = [float(line.split()[0]), float(line.split()[1])]
		print "currentPosition = ", currentPosition, "lastPosition = ", lastPosition
		if (currentPosition == position) or (lastPosition < position and currentPosition > position):
			cols = [time]
			cols.extend(line.split())
			for col in cols :
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

f.close()
print "e\n"
