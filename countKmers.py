## Script to count all k-mers in a fasta file
## Gregory Moyerbrailean, Wayne State University 2014
##
## How to use
## 	$ kmers.py [-h] [-o OUT] [-v] file kmer
## 	positional arguments:
## 	  file               a fasta file to read in
## 	  kmer               k-mer length
## 	optional arguments:
## 	  -h, --help         show this help message and exit
## 	  -o OUT, --out OUT  output file (default: kmers.txt)
## 	  -v, --verbose      verbose mode (default: off)

import argparse

## Get options from command line, print help if necessary
parser = argparse.ArgumentParser(prefix_chars='-',
								 description='Count k-mers in a fasta file')
parser.add_argument("file", metavar = "file", type = str, 
					help = "a fasta file to read in")
parser.add_argument("kmer", metavar = "kmer", type = int,
					help = "k-mer length")
parser.add_argument("-o", "--out", default = "kmers.txt", type = str,
					help="output file (default: kmers.txt)")
parser.add_argument("-v", "--verbose", help="verbose mode (default: off) ",
                    action="store_true")
args = parser.parse_args()
seqFile = args.file
outFile = args.out
k = args.kmer
v = args.verbose

## Get input file
if v:
	print "Reading file..."
fasta = open(seqFile).read()

# Process each sequence individually
kMers = dict()
recs=fasta.split('>')
if '' in recs: recs.remove('')

l = len(recs)
n = 1
for r in recs:
	if v:
		print "Processing %d of %d" % (n, l)
	r = r.split('\n')
	if '' in r: r.remove('')
	seqName = r[0]
	seqStr = ''.join(r[1:])

	i = 0
	while i+k < len(seqStr):
		ss = seqStr[i:i+k]
		ss = ss.upper() # normalize case
		if ss in kMers:
			kMers[ss] += 1
		else:
			kMers[ss] = 1
		i += 1
	n += 1

## Output the results
if v:
	print "Writing results..."
fd = open(outFile, 'w')
for k,v in kMers.items():
	fd.write('\t'.join([k,str(v)]) + '\r\n')
fd.close()
