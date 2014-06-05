## Script to annotate a bed file with the distance to the nearest TSS file
## Adds the distance to the file as a new column, the rightmost field
## Gregory Moyerbrailean, Wayne State University 2013
##
## Required files
##   TSS locations: a two-column (pref. tab-separated) file specifying 
##    				chromosome and TSS position.
## How to use
## 	$ python getNearestTss.py input.bed output.bed
## 	input.bed must be a bam file without a header. As long as it is tab-
##	separated and has chr, start, and end, it doesn't matter how many
##  columns it has - getNearestTss will append a column to the far right
##  with the distance to the nearest TSS for each record

import sys, gzip

def getTssLocs(fName='hg19.tss.loci.gz', delim='\t'):
	'''Return a dict of TSS positions (k - chr, v - position list)'''
	fd = gzip.open(fName,'rb')
	tssFile = fd.read().strip()
	tssFile = tssFile.split('\n')
	if ( '' in tssFile ): tssFile.remove('')
	tssDict = dict()
	for tss in tssFile:
		c, p = tss.split(delim)
		if c in tssDict:
			tssDict[c].append(int(p)) # make non-redund?
		else:
			tssDict[c] = [int(p)]
	# Sort the positions in case the file wasn't sorted
	for k in tssDict.keys():
		tssDict[k].sort()
	return tssDict

## For each loci, use binary search to find closest TSS
def getDistToTss(pos,tss):
	'''Return the distance of a loci to the nearest tss'''
	if ( pos < tss[0] ): # Before first TSS
		return ( tss[0] - pos )
	elif ( pos > tss[-1] ): # After last TSS
		return ( pos - tss[-1])
	else:
		mn=0
		while True:
			mx=len(tss)
			if len(tss) == 2: # Find closest
				return min(abs(pos-tss[0]), abs(pos-tss[1]))
			elif len(tss) == 1:
				return abs(pos-tss[0])
			m = (mn + mx) / 2
			if tss[m] < pos:
				tss = tss[m+1:]
			elif tss[m] > pos:
				tss = tss[:m]
			else:
				return 0 # loci is a TSS!

## Get file names from argv
fName = sys.argv[1]
oName = sys.argv[2]

# Get TSS locations
tssName='hg19.tss.loci.gz'
if len(sys.argv) > 3:
	tssName = sys.argv[3]
tDict = getTssLocs(tssName)

## Open input. Can be gzipped or not
if ( fName[-3:] == '.gz' ):
	fd = gzip.open(fName, 'rb')
else:	
	fd = open(fName)
loci = fd.read().split('\n')
fd.close()

delim = '\t' # Delimiter TODO: make accept custom delimiters
res = ''	 # Empty string to hold the output

## Check if the first line is a header (note: this is not bed-standard).
## If so, grab it so as to add it to the output
header=''
lineOne = loci[0]
lineOneSplit = lineOne.split(delim)
if not (lineOneSplit[1].isdigit() and lineOneSplit[2].isdigit()):
	header = lineOne + delim + 'tss.distance' + '\n'
	loci = loci[1:]
res += header


## Annotate with the nearest TSS
if '' in loci: loci.remove('')
for line in loci:
	l = line.split(delim)
	assert len(l) > 2 # At least 3 cols (bed format)
	c = l[0] # chr
	p = ((int(l[1]) + int(l[2])) / 2 ) # middle of position, e.g, motif
	
	## Some file may have non-standard chrs, e.g, "chr11_gl000202_random"
	## If there aren't TSS data for them, skip
	try:
		dist = getDistToTss(p,tDict[c])
	except KeyError:
		continue
	
	nLine = line + delim + str(dist)
	res += nLine + '\n'

## Output new file with annotations appeneded as new column
if ( oName[-3:] == '.gz' ):
	fd = gzip.open(oName, 'wb')
else:	
	fd = open(oName, 'w')
fd.write(res)
fd.close()
