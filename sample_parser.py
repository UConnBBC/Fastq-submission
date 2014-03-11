#	sample_parser.py
#
#	description: Parses FASTQ read-files into individual files by sample name, 
#		using a merged FASTA file to define sample:header relationships.
#
#! /usr/bin/python

import time
import gzip
import sys
import re
import os

class FASTQrecord:
	hash = -1
	source = ''
	
	def __init__(self, h, s):
		self.hash = h #holds a hash of the header of a record from a FASTQ file
		self.source = s #the full data associated with a record (4 lines), header included, for output
			
args = sys.argv
narg = len(sys.argv)
if narg != 3:
	exit('Invalid number of arguments. Usage is:\n\t'+args[0]+' [sample_definition_file] [read_file(s)]\n\te.g. '+args[0]+' seqs.fna R1.fastq,R2.fastq\n\tNOTE: multiple read-files should be comma-delimited with no spaces. It is reccomended that the read files have simple names as these are used to generate the names of the corresponding output, split by sample.')

sampledefs = args[1]
fastqfiles =  args[2].split(',')#['R1.fastq','R2.fastq']

t0 = time.time()
with open(os.getcwd()+'/'+sampledefs) as f:
	lines = f.readlines()
f.close()

print 'Loaded '+sampledefs+'. Now pairing sample IDs and headers...'
i = 0
di = {'null':[],}
for l in lines:
	if l[0] == '>': 
		m = re.match( r'(.*?)_(\d+) (.*?) ',l) #match sample id and header
		m1 = m.group(1) #sample ID
		m1 = m1[1:len(m1)] #trim off '>'
		#m2 = m.group(2) #this is just the position of sequence in def file, not useful
		m3 = m.group(3) #header to be found in read files
		
		try:
			di[m1].append(hash(m3))
		except:	
			di[m1] = [hash(m3),]
		
		i+=1 #counts the number of sample:header pairings made
		
del di['null']

print 'Sample information from '+sampledefs+':'
foo = 0
for k in di.keys():
	uids = di.get(k)
	j = len(uids)
	foo = foo + j
	print '\tSample key \''+k+'\' has '+str(j)+' sequences.'
print '\tTOTAL # of sequences across all samples is '+str(foo)
print 'Paired sample:header information for '+str(len(di.keys()))+' sample IDs.'

#Now that we have the linkage between headers and samples... find the headers in each
#	read file (fastq) and then write them to the appropriate file based on sample ID.
#	Each sample ID will thus have an individual fastq file for each read file provided via
#	comnmand-line argument

for fs in fastqfiles: 
	#Goes through read file, one at a time, to look for sequences we have a sample ID:header pair
	t3 = time.time()
	din = di
	dout = {'' : [],} 
	#If header on record in read file has a matching sampleID, put in list for that sampleID key
	n = 0
	l = 0
	missing = [] 
	#If a record in a read file does not have a sample ID, it goes into this list
	flag = -1
	suffix = fs.split('.')[0]
	p = []
	#This list 'p' is the parsed list of FASTQrecords for the current read file
	with open(os.getcwd()+'/'+fs, "rU") as inner:
		while True:
			try: 
				#This logic relies on each FASTQ record being 4 lines, delineated by \n
				#	The first line is assumed to contain header information
				l1 = inner.next()
				m = re.match(r'^@(.*?) ',l1)
				p.append(FASTQrecord(hash(m.group(1)),l1+inner.next()+inner.next()+inner.next()))
				n+=1
   			except StopIteration: 
   			 	break 
   	#Loops through every record in the read file to try and match to the sampleID:header definitions
   	print 'Attempting to define '+str(n)+' records in '+fs+'...'
	for record in p:
		flag = -1 #means record has not yet been found
		#r1 = time.time() #Uncomment to check length of time required for each record def
		h = record.hash
		for k in din:
			if h in din[k]:
				flag = 1 #flag = 1 means record was found
				l+=1 #counts
				try:
					dout[k].append(record) #Add the record to the proper sample output list
				except:
					dout[k] = [record,]	#If key in output dict does not exist, create it
				break #stop searching for a sampleID for this record, we found it
			
		if flag != 1: #If cannot find current record in dict of lists of hashes
			missing.append(record) 
			#add to list of records we failed to have a def for
		#else:
			#print (time.time() - r1)*1000 
			## how long (ms) does it take to find a single record?
			
	t4 = time.time()
	m, s = divmod((t4-t3), 60)	
	h, m = divmod(m, 60)
	print 'Defined '+str(l)+'/'+str(n)+' ('+str(round(100.*l/n))+'%) records from '+fs+" in %d:%02d:%02d" % (h, m, s)	

	count = 0 
	#We have looped through all records in the current file so write the dicts to file
	#	each key gets a file, plus one file for undefined records, for each read file
	print 'Generating output...'
	for k in dout:
		if k != '':
			cur = len(dout[k])
			count+=cur
			#output format is SAMPLEID_READFILEID.fastq. e.g. A3Nc_R1.fastq 
			with gzip.open(os.getcwd()+'/'+k+'_'+suffix+'.fastq.gz', 'wb') as out:
				for o in dout[k]:
					out.write(o.source)
			out.close()
			print '\t wrote '+str(cur)+' sequences to '+k+'_'+suffix+'.fastq'
	print '\tTOTAL defined records written to file: '+str(count)
	with gzip.open(suffix+'_ERR.fastq.gz', 'wb') as out:
		for o in missing:
			out.write(o.source)
	out.close()
	print '\tMISSING definitions for '+str(len(missing))+' records in '+fs+'. Wrote to '+suffix+'_ERR.fastq.'
runt = time.time()-t0
m, s = divmod(runt, 60)
h, m = divmod(m, 60)
print "Done! Total run time: %d:%02d:%02d" % (h, m, s)
exit(0)
