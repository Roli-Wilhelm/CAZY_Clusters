#!/usr/bin/python
import sys, os, re, getopt, glob, numpy as np
import timeit
import fileinput
import collections
import os.path
from Bio import SeqIO
from subprocess import call
from itertools import izip_longest
from itertools import groupby
from operator import itemgetter

start = timeit.default_timer()

Usage = """

PURPOSE:	
	Identify and recover clusters of CAZymes by running a .fna file through CAZy, dbCANN (HMM) and Cardenas models for AA families and dyp (HMM), then filter to include only clusters and recover sequences

REQUIRED ARGUMENTS:
		#BIOLOGICAL PARAMETERS
		-c	Number of CAZymes present to qualify as a cluster  					[DEFAULT: 3]

		-A	Add X nt to the front of your CAZY cluster sequence  					[DEFAULT: 0]

		-f	Filter out all partial genes as defined by prodigal (KEEP_PARTIAL | REMOVE_PARTIAL)   	[DEFAULT: KEEP_PARTIAL]

                -e      Maximum E-value  									[DEFAULT : 1e-05]

		-d	Force CAZymes to occur sequentially on contig to be considered a cluster		[DEFAULT: N]

		-N	Search for potential novel CAZymes by returning only clusters that: 			[DEFAULT: N]
			possess an unannotated ORF within a sequence of CAZymes  (Y | N)


		NOTE:	In cases where multiple, non-overlapping CAZyme hits are found within a single predicted ORF, they are all reported. However, the order in which they appear may not be correct. The ordering of CAZymes
			not found in the same ORF (aka ordinal) are reported in the correct order.

		#COMPUTATIONAL PARAMETERS
                -I	A file in .tsv format with two columns :  SampleID     Fasta File (DNA)

		-O	Specify OUTPUT DIRECTORY

		-p	Number of processors to use


REQUIREMENTS:

	SOFTWARE:
		-DIAMOND (built with v. 0.7.9)

		-HMMER (built with v. 3.0)

		-PRODIGAL (built with  v. 2.6.2)

		-bedtools (built with v 2.17.0)

		
Usage:  ./CAZY_CLUSTER_SIMPLIFIED.py -I CAZY_OUTPUT_LIST.tsv -O OUTPUT -c 3 -f REMOVE_PARTIAL -d Y
			
"""

# Read command line args
myopts, args = getopt.getopt(sys.argv[1:],"I:p:O:f:e:N:A:c:d:")

LIST=''
PROCESSORS=''
OUTPUT=''
FILTER=''
E_VALUE=''
FORE=''
CLUSTER=''
DISTANCE=''
NOVEL=''

for o, a in myopts:
    if o == '-I':
        LIST= a
    if o == '-p':
        PROCESSORS= a
    if o == '-O':
        OUTPUT= a
    if o == '-e':
        E_VALUE= a
    if o == '-f':
        FILTER = a
    if o == '-c':
        CLUSTER= a
    if o == '-A':
        FORE= a
    if o == '-d':
        DISTANCE  = a
    if o == '-N':
        NOVEL  = a


# Set Defaults
if len(sys.argv)<2:
        print Usage
        sys.exit()

if not CLUSTER:
	CLUSTER = 3

if not E_VALUE:
	E_VALUE = float(1e-05)

if not FORE:
	FORE = "0"

if not DISTANCE:
	DISTANCE = "N"

if not NOVEL:
	NOVEL = "N"
	DISTANCE == "Y" # Distance has to be set to yes in order to identify cluster sequences

if not FILTER:
	FILTER = "KEEP_PARTIAL"

# Build Functions
def SHARED_SET(a, b):
    return not set(a).isdisjoint(b)

def GUNZIP(x):
        print "\n-- Unzipping "+x+" --\n"
        call(["pigz","-d","-p",PROCESSORS,x])

        #Remove ".gz" file extension from name
        x = re.sub(".gz", "", x)
        return x

def GZ(x):
        print "\n-- Zipping "+x+" --\n"
        call(["pigz","-p",PROCESSORS, x])

def grouper(n, iterable, fillvalue=None):
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def clean_list(x):
	x = re.sub("'","",str(x))
	x = re.sub("\[","",x)
	x = re.sub("\]","",x)
	return x


## CREATE OUTPUT DIRECTORY
if os.path.exists(OUTPUT):
	pass
else:
	os.mkdir(OUTPUT)


## IMPORT FILES into DICTIONARY
FILE_DICT = {}

with open(LIST) as f:
       	next(f)

        for line in f:
               	line = line.strip("\r\n")
       	        line = line.split("\t")

		if re.search(".gz",line[1]):
			line[1] = GUNZIP(line[1])
			ZIP_ME = "YES"

		FILE_DICT[line[0]] = line[1]


'''#################################
## PREDICT PROTEINS WITH PRODIGAL
#################################
print "\n-----------------------------------\nRUNNING PRODIGAL\n-----------------------------------\n"

for key, value in FILE_DICT.iteritems():
	FILE = value
	NAME = key

	os.system(' '.join([
		"prodigal",
		"-a",
		OUTPUT + "/" + NAME + ".faa",
		"-o",
		OUTPUT + "/" + NAME + ".prodigal",
		"-i",
		FILE,
		"-p meta",
		">",
		"/dev/null 2>&1"
	]))

	###################################
	## RUN THE VARIOUS ANALYSIS METHODS
	###################################

	print "Processing: "+re.sub("\[|\]","",str(NAME))+"\n"

	# BLAST
	print "\n--------------------\nYOU ARE NOW BLASTING:  "+NAME+"\n--------------------\n"

	if os.path.exists("TEMP"):
		pass
	else:
		os.mkdir("TEMP")

	os.system(' '.join([
		"diamond blastp",
		"-d",
		"/home/roli/CAZy/CAZy_2015_08_19",
		"-q",
		OUTPUT + "/" + NAME + ".faa",
		"-a",
		OUTPUT + "/" + NAME + ".daa",
		"-p",
		PROCESSORS,
		"-t TEMP"
	]))

	os.system(' '.join([
		"diamond view",
		"-a",
		OUTPUT + "/" + NAME + ".daa",
		"-o",
		OUTPUT + "/" + NAME + ".di.blast.foo",
	]))

	try:
		input = open(OUTPUT + "/" + NAME + ".di.blast.foo","r")
		output = open(OUTPUT + "/" + NAME + ".di.blast","w")

		for line in input:
			line = line.strip()
			output.write(line+"\t"+NAME+"\n")

		input.close()
		output.close()
	
	except:
		print "No BLAST Hits"

	# AA HMM SCAN
	print "\n-----------------------------------------------\nYOU ARE NOW SCANNING FOR AA GENES WITH HMMSCAN:  "+NAME+"\n----------------------------------------------\n"

	HMM_LIST = ["ary.oh.oxi.hmm","dyp2.hmm","lacc.hmm","sm.lacc.hmm","van.oh.oxi.hmm","ver.perox.hmm"]
	HMM_DIRECTORY = "/home/roli/CAZy/AA_hmm/"

	os.system(' '.join([
		"cp",
		'/usr/local/hmmscan-parser.sh',
		"."
	]))

	## Compiled Output file
	output = open(OUTPUT + "/" + NAME + ".AA.tbl","w")

	for HMM in HMM_LIST:
		# Run Scan		
		os.system(' '.join([
			"hmmscan",
			HMM_DIRECTORY+HMM,
			OUTPUT + "/" + NAME + ".faa",
			">",
			OUTPUT + "/" + NAME + "." + re.sub(".hmm","",HMM) + ".AA"
		]))

		# Find HMM model length and store in file used by dbCAN 'hmmscan-parser.sh
		# TAKEN FROM CARDENAS' HMM.search.and.parse.and.extract.py
		hmm_fileout = open(OUTPUT+"/"+HMM+".length.txt",'w')
		hmm_filein = open(HMM_DIRECTORY+HMM,'r')

		for line in hmm_filein:
			if line.startswith('NAME'):
				line = line.strip('\n')
				line = line.split(' ')
				name = line[2]
				hmm_fileout.write('%s\t' %name)

			else:
				if line.startswith('LENG'):
					line = line.strip('\n')
					line = line.split(' ')
					len = line[2]
					hmm_fileout.write('%s\n' %len)
				else:
					continue

		hmm_fileout.close()
		hmm_filein.close()
		os.system(' '.join(['cp',OUTPUT+"/"+HMM+".length.txt",'./all.hmm.ps.len']))  #THIS FILE IS NEEDED BY dbCAN HMM_parser

		# Parse Results

		os.system(' '.join([
			'sh',
			'/usr/local/hmmscan-parser.sh',
			OUTPUT + "/" + NAME + "." + re.sub(".hmm","",HMM) + ".AA",
			">",
			OUTPUT + "/" + NAME + "." + re.sub(".hmm","",HMM) + ".AA.tbl.foo",
		]))

		# Append the HMM Name + Compile
		try:
			input = open(OUTPUT + "/" + NAME + "." + re.sub(".hmm","",HMM) + ".AA.tbl.foo","r")
	
			for line in input:
				line = line.strip()
				output.write(line+"\t"+re.sub(".hmm","",HMM)+"\t"+NAME+"\n")
	
		except:
			print "No Cardenas HMM's found"
	input.close()
	output.close()

	# dbCAN HMM SCAN
	print "\n---------------------------------------\nYOU ARE NOW SCANNING dbCAN WITH HMMSCAN:  "+NAME+"\n---------------------------------------\n"

	os.system("cp /home/roli/CAZy/dbCAN/all.hmm.ps.len .")
	os.system("cp /home/roli/CAZy/dbCAN/hmmscan-parser.sh .")

	# Run Scan
	os.system(' '.join([
		"hmmscan",
		"/home/roli/CAZy/dbCAN/dbCAN-fam-HMMs.txt",
		OUTPUT + "/" + NAME + ".faa",
		">",
		OUTPUT + "/" + NAME + ".dbCAN"			
	]))

	# Parse Results
	os.system(' '.join([
		"sh",
		"hmmscan-parser.sh",
		OUTPUT + "/" + NAME + ".dbCAN",
		">",
		OUTPUT + "/" + NAME + ".dbCAN.tbl.foo"
	]))

	# Append the HMM Name + Compile
	try:
		input = open(OUTPUT + "/" + NAME + ".dbCAN.tbl.foo","r")
		output = open(OUTPUT + "/" + NAME + ".dbCAN.tbl","w")

		for line in input:
			line = line.strip()
			output.write(line+"\t"+NAME+"\n")
	
		input.close()
		output.close()

	except:
		print "No dbCAN HMM's found"

	os.system("rm "+OUTPUT+"/*.AA")
	os.system("rm "+OUTPUT+"/*.daa")
	os.system("rm "+OUTPUT+"/*.foo")
	os.system("rm "+OUTPUT+"/*.dbCAN")

##########################################
## PARSE AND COMPILE BLAST and HMM OUTPUTS 
##########################################

# HEADER INFO
header = "QUERY, HIT, CAZY_FAMILY, %_IDENTITY_OR_COVERAGE, LENGTH, Q_START, Q_END, ID_METHOD, E_VALUE, SOURCE\n"
header = re.sub(",","\t",header)

output1 = open(OUTPUT+"/ALL_CAZYmes.compiled.out","w")
output1.write(header)

# GO THROUGH EACH METHOD
for key, value in FILE_DICT.iteritems():
	FILE = value
	NAME = key

	print "\n-----------------\nNOW PARSING FILE:  "+NAME+".CAZYmes.compiled.out\n-----------------\n"

	output2 = open(OUTPUT+"/"+NAME+".CAZYmes.compiled.out","w")
	output2.write(header)

	# BLAST
	# LOAD CAZy Dictionary	
	CAZY_ID = {}
	with open(r"/home/roli/CAZy/CAZy_2015_08_19.fixed.pkl", "rb") as input:
		CAZY_ID = pickle.load(input)

	# ID with CAZy Dictionary
	try:
		with open(OUTPUT + "/" + NAME + ".di.blast") as f:
			for line in f:
				line = line.strip()
				line = line.split("\t")

				QUERY = line[0]
				HIT = line[1]
				CAZY_FAMILY = CAZY_ID[HIT][1]
				IDENTITY = line[2]
				LENGTH = line[3]
				Q_START = line[6]
				Q_END = line[7]
				ID_METHOD = "BLAST"
				E_VALUE = line[10]

				output1.write(QUERY+"\t"+HIT+"\t"+CAZY_FAMILY+"\t"+IDENTITY+"\t"+LENGTH+"\t"+Q_START+"\t"+Q_END+"\t"+ID_METHOD+"\t"+E_VALUE+"\t"+NAME+"\n")
				output2.write(QUERY+"\t"+HIT+"\t"+CAZY_FAMILY+"\t"+IDENTITY+"\t"+LENGTH+"\t"+Q_START+"\t"+Q_END+"\t"+ID_METHOD+"\t"+E_VALUE+"\t"+NAME+"\n")
	except:
		pass

	# AA HMM SCAN
	try:
		with open(OUTPUT + "/" + NAME + ".AA.tbl") as f:
			for line in f:
				line = line.strip()
				line = line.split("\t")
				
				QUERY = line[0]
				HIT = line[1]
				E_VALUE = line[2]
				Q_START = line[3]
				Q_END = line[4]
				LENGTH = str(abs(int(Q_END)-int(Q_START))+1)
				ID_METHOD = "CARDENAS_AA"
				IDENTITY = str(round(float(line[7]),2))
				CAZY_FAMILY = line[8]
				
				output1.write(QUERY+"\t"+HIT+"\t"+CAZY_FAMILY+"\t"+IDENTITY+"\t"+LENGTH+"\t"+Q_START+"\t"+Q_END+"\t"+ID_METHOD+"\t"+E_VALUE+"\t"+NAME+"\n")
				output2.write(QUERY+"\t"+HIT+"\t"+CAZY_FAMILY+"\t"+IDENTITY+"\t"+LENGTH+"\t"+Q_START+"\t"+Q_END+"\t"+ID_METHOD+"\t"+E_VALUE+"\t"+NAME+"\n")
	except:
		pass
			
	# dbCAN HMM SCAN
	try:
		with open(OUTPUT + "/" + NAME + ".dbCAN.tbl") as f:
			for line in f:
				line = line.strip()
				line = line.split("\t")
			
				QUERY = line[0]
				HIT = line[1]

				# To Deal with the Super Weirdness of parsing (the lines look normal in text editing, but do not split properly and cause a value error)
				try:
					E_VALUE = str(line[2])
					Q_START = line[3]
					Q_END = line[4]
					LENGTH = str(abs(int(Q_END)-int(Q_START))+1)
					ID_METHOD = "dbCAN"
					IDENTITY = str(round(float(line[7]),2))
					CAZY_FAMILY = re.sub(".hmm","",line[1])
			
					output1.write(QUERY+"\t"+HIT+"\t"+CAZY_FAMILY+"\t"+IDENTITY+"\t"+LENGTH+"\t"+Q_START+"\t"+Q_END+"\t"+ID_METHOD+"\t"+E_VALUE+"\t"+NAME+"\n")
					output2.write(QUERY+"\t"+HIT+"\t"+CAZY_FAMILY+"\t"+IDENTITY+"\t"+LENGTH+"\t"+Q_START+"\t"+Q_END+"\t"+ID_METHOD+"\t"+E_VALUE+"\t"+NAME+"\n")
			
				except:
					print line

	except:
		pass

output2.close()
output1.close()'''

###########################################
## FILTER BLAST RESULTS AND ASSIGN CLUSTERS
###########################################
for key, value in FILE_DICT.iteritems():
	FILE = value
	NAME = key

	print "FILE NAME IS: "+FILE
	print "NAME IS: "+NAME

	GENE_INFORMATION ={}

	print "\n-------------------\nPROCESSING:  "+NAME+"\n-------------------\n"

	##################################################
	## FILTER BASED ON PRODIGAL PREDICTED COMPLETENESS
	##################################################

	# Open Prodigal FAA output to parse information on completeness into dictionary
	with open(OUTPUT+"/"+NAME+".faa","r") as f:
		for line in f:
			if line.startswith(">"):
				line = line.strip()
				line = line.split("#")

				# STORE INFO
				GENE = line[0]
				GENE = GENE.strip()
				GENE = re.sub(">","",GENE)

				GENE_START = line[1]
				GENE_START = GENE_START.strip()

				GENE_END = line[2]
				GENE_END = GENE_END.strip()

				SENSE = line[3]
				SENSE = SENSE.strip()

				# FIND COMPLETENESS INFO
				ID = line[4]
				ID = ID.strip()
				PARTIAL = ID.split(";")[1]
				PARTIAL = PARTIAL.split("=")[1]

				# CREATE FILTER DICTIONARIES
				if FILTER == "KEEP_PARTIAL":
					if PARTIAL != "11":	#DEFAULT TO NOT KEEP GENES WITH PARTIAL ON BOTH ENDS
						GENE_INFORMATION[GENE] = [GENE_START,GENE_END,SENSE]

				elif FILTER == "REMOVE_PARTIAL":
					if PARTIAL == "00":
						# By storing individual contig hits, we filter all partial genes
						GENE_INFORMATION[GENE] = [GENE_START,GENE_END,SENSE]


	############################################
	## FILTER BLAST HITS and PREPARE CONTIG DICT
	############################################
	OPEN_ME = OUTPUT+"/"+NAME+".CAZYmes.compiled.out"
	output = open(OUTPUT+"/"+NAME+".CAZYmes.filtered.out","w")

	header = "QUERY, HIT, CAZY_FAMILY, %_IDENTITY_OR_COVERAGE, LENGTH, Q_START, Q_END, ID_METHOD, E_VALUE, SOURCE\n"
	header = re.sub(",","\t",header)
	output.write(header)

	# INITIALIZE DICTIONARIES	
	CONTIG_DICT = {}
	GENE_DICT = {}		# Used to verify that CAZYmes from within the same ORF are different
	CLUSTER_DICT = {}

	# POPULATE CONTIG_DICT and FILTER BY PARTIAL
	with open(OPEN_ME) as f:
		next(f)

		for line in f:
			line = line.strip()
			original_line = line	# SAVE TO WRITE to FILTERED FILE
			line = line.split("\t")

			# PARSE NAME OF CONTIG (i.e. separate CONTIG from ORDINAL ID)
			if np.size(line[0].split("_")) > 2:
				CONTIG = line[0].split("_")
				ORDINAL = CONTIG[len(CONTIG)-1]
				CONTIG = '_'.join(CONTIG[0:(len(CONTIG)-1)])

			else:
				CONTIG = line[0].split("_")[0]
				ORDINAL = line[0].split("_")[1]			

			GENE = line[0]
			HIT_START = int(line[5])
			HIT_END = int(line[6])

			## CAZYme HITS FROM HMMs NEED TO BE STRIPPED (or else they're counted as unique CAZymes)
			CAZYme = line[2] 
			CAZYme = CAZYme.strip()

			# FILTER E_VALUE
			if E_VALUE > float(line[8]):
				if GENE_INFORMATION.has_key(GENE): ## FILTER OUT THOSE WHICH WERE PARTIAL (IF CHOSEN)
					if GENE_DICT.has_key(GENE):

						## SELECT TOP BLAST HIT(S) || NOTE: THIS SCRIPT PERMITS MULTIPLE CAZYMES PER PRODIGAL ORF **ONLY** IF THEY ARE i) DIFFERENT CAZYMES and ii) DO NOT OVERLAP BY MORE THAN 30 aa
						if CAZYme not in GENE_DICT[GENE][0]:
							RECORD_START = GENE_DICT[GENE][1]
							RECORD_END = GENE_DICT[GENE][2]
							
							# Calculate Total Length of Record Gene
							#RECORD_LENGTH = float(len(range(RECORD_START, RECORD_END)))

							# Calculate Length of Overlap
							OVERLAP = list(set(range(HIT_START, HIT_END)) & set(range(RECORD_START, RECORD_END)))
							OVERLAP_LENGTH = float(len(OVERLAP))

							#if OVERLAP_LENGTH/RECORD_LENGTH < 0.05:  
							if OVERLAP_LENGTH < 30:

								# Add to Main Dictionary
								CONTIG_DICT[CONTIG] = CONTIG_DICT[CONTIG] + [GENE,CAZYme,HIT_START,HIT_END,ORDINAL]

								# Write to 'filtered' file
								output.write(original_line+"\n")

								# Modify GENE_DICT information
								GENE_DICT[GENE][0].append(CAZYme)
								TOTAL_RANGE = sorted(list(set(range(HIT_START, HIT_END) + range(RECORD_START, RECORD_END))))
								GENE_DICT[GENE][1] = TOTAL_RANGE[0]		# Set New RECORD_START
								GENE_DICT[GENE][2] = TOTAL_RANGE[len(TOTAL_RANGE)-1]	# Set New RECORD_END


					else:
						if not CONTIG_DICT.has_key(CONTIG):
							CONTIG_DICT[CONTIG] = [GENE,CAZYme,HIT_START,HIT_END,ORDINAL]
							output.write(original_line+"\n")

						else:
							CONTIG_DICT[CONTIG] = CONTIG_DICT[CONTIG] + [GENE,CAZYme,HIT_START,HIT_END,ORDINAL]
							output.write(original_line+"\n")
		
						GENE_DICT[GENE] = [[CAZYme],HIT_START,HIT_END]	
	output.close()	


	##################
	## ASSIGN CLUSTERS 
	##################
	output = open(OUTPUT+"/"+NAME+".CAZYme.clusters.summary", "w")

	### ADD FLAG
	header = "SAMPLE_NAME, CAZY_COUNT, CAZY_LIST, CONTIG_ID, ORDINAL_IDs, GENE_COORDINATES, CLUSTER_COORDINATES, STRAND_SENSE\n"

	header = re.sub(",","\t",header)
	output.write(header)

	## RUN THROUGH DATA CONTIG BY CONTIG
	for key, value in CONTIG_DICT.iteritems():
		if np.size(value)/5 >= int(CLUSTER):	# FILTER BASED ON CLUSTER SIZE (multiply by 5 since each GENE hit contains five elements... tried to make this less awkward, but couldn't find an answer)

			# PARSE INFORMATION AND PERFORM FINAL FILTERS 
			CONTIG = key
			CAZYmes = []
			ORDINALS = []
			SENSE = []
			
			# BASE FILTERING OF CONTIGS and CLUSTERS BASED ON ORDINALS (i.e. the gene order assigned by PRODIGAL)
			for GENE in grouper(5,value):
				ORDINALS.append(int(GENE[4]))

			### SELECT FOR SEQUENTIAL CLUSTERS AND NOVEL CLUSTERS - IF SPECIFIED
			### 
			if DISTANCE == "Y":
				KEEP_CONTIG = "N"
				ORDINALS_TO_KEEP = []
				SEQUENCE_TEST = list(set(ORDINALS))
				SEQUENCE_TEST = sorted(SEQUENCE_TEST)

				for key, group in groupby(enumerate(SEQUENCE_TEST), lambda (index, item): index - item):
					group = map(itemgetter(1), group)

					## Only consider sequential clusters greater than or equal to cluster cut-off
					if len(group) >= CLUSTER:
						KEEP_CONTIG = "Y"  #IF NO GROUP IS > THAN CLUSTER SIZE: DISCARD	
	
						### INCLUDE CLUSTERS WITH FILTER FOR NOVEL 
						if NOVEL == "Y":
							GROW = "Y"

							while GROW == "Y":
		
								if (max(group) + 1) in SEQUENCE_TEST:
									group =  group + [(max(group) + 1)]

								elif (max(group) + 2) in SEQUENCE_TEST:
									group =  group + [(max(group) + 1),(max(group) + 2)]

								elif (min(group) - 2) in SEQUENCE_TEST:
									group =  group + [(min(group) - 1),(min(group) - 2)]
									group = [x for x in group if x > 0]

								else:
									ORDINALS_TO_KEEP.extend(group)
									GROW = "N"
						else:
							ORDINALS_TO_KEEP.extend(group)
	
				if KEEP_CONTIG == "Y":
					## REMOVE ORDINALS NOT IN CONTIGUOUS CLUSTERS
					items = set(ORDINALS_TO_KEEP)
					ORDINALS = [i for i in ORDINALS if i in items]
				else: 
					ORDINALS = []

			else:
				KEEP_CONTIG = "Y"  ## DEFAULT IS TO ACCEPT CONTIG

			GENE_LIST = []
			count = 0
			for GENE in grouper(5,value):
				if int(GENE[4]) in ORDINALS:
					# GET INDVIDUAL GENE SENSE from GENE_INFORMATION
					SENSE.append(GENE_INFORMATION[GENE[0]][2])

					#CAPTURE LIST OF CAZYmes
					CAZY = GENE[1]
					CAZY = CAZY.strip()				
					CAZYmes.append(CAZY)
	
					#SUMMARIZE START AND STOP COORDINATES
					GENE_START = GENE_INFORMATION[GENE[0]][0]
					GENE_END = GENE_INFORMATION[GENE[0]][1]

					COORDINATE = GENE_START+"-"+GENE_END
	
					if count == 0:
						MIN = int(GENE_START)
						MAX = int(GENE_END)
						count = count + 1
		
					else:		
						if int(GENE_START) < int(MIN):
							MIN = GENE_START
	
						if int(GENE_END) > int(MAX):
							MAX = GENE_END

					if COORDINATE not in GENE_LIST:
						GENE_LIST.append(COORDINATE)

			COUNT = len(ORDINALS)

			## Clean List for Writing to File
			CAZYmes = clean_list(CAZYmes)
			ORDINALS = clean_list(ORDINALS)
			GENE_LIST = clean_list(GENE_LIST)
			SENSE = clean_list(SENSE)
	

			if KEEP_CONTIG == "Y":
				output.write(NAME+"\t"+str(COUNT)+"\t"+str(CAZYmes)+"\t"+CONTIG+"\t"+str(ORDINALS)+"\t"+str(GENE_LIST)+"\t"+str(MIN)+"-"+str(MAX)+"\t"+str(SENSE)+"\n")


########################################################################
### RETRIEVE FULL DNA SEQUENCE FOR EACH CONTIG CONTAINING A CAZY_CLUSTER
########################################################################

# Cycle Through All Files and Create Dictionary of PULs, Then Create and Excecute BedTools to Recover the Raw Sequences
for key, value in FILE_DICT.iteritems():
	FASTA_FILE = value
	NAME = key
	CLUSTER_FILE = OUTPUT+"/"+NAME+".CAZYme.clusters.summary"

	output = open(OUTPUT+"/"+NAME+".full.bed", "w")
	FULL_INFO_DICT = {}

        with open(CLUSTER_FILE) as f:
		next(f)

		# Cycle Through All Files and Create Dictionary of PULs, Then Create and Excecute BedTools to Recover the Raw Sequences
		for line in f:
			line = line.strip()
			line = line.split("\t")

			NUMBER_OF_CAZYMES = line[1]
			COORDS = line[6]
			COORDS = COORDS.split("-")		
			START = COORDS[0]
			STOP = COORDS[1]

			# SAVE INFO TO WRITE BETTER FASTA LABELS in BEDTOOLS OUTPUT
			CAZYME = line[2]
			CAZYMES = CAZYME.split(",")
			CONTIG = line[3]

			FULL_INFO_DICT[CONTIG] = CAZYMES			

			# PREPARE BEDTOOLS "BED" FILE
			if NUMBER_OF_CAZYMES >= CLUSTER:
				START = START.strip()
				STOP = STOP.strip()
		
				TEST = int(START) 

				if TEST-int(FORE) > 0:
					START = str(TEST-int(FORE))
				else:
					START = "1"

				output.write(CONTIG+"\t"+START+"\t"+STOP+"\n")
		
	output.close()

	os.system(' '.join([
		"fastaFromBed",
		"-fi",
		FASTA_FILE,
		"-bed",
		"./"+OUTPUT+"/"+NAME+".full.bed",
		"-fo",
		"./"+OUTPUT+"/"+NAME+".recovered."+FORE+".fa"
	]))

	# Write Better Names for Retrieved Records
	output = open(OUTPUT+"/"+NAME+".recovered.foo","w")

	for record in SeqIO.parse(open(OUTPUT+"/"+NAME+".recovered."+FORE+".fa",'rU'),'fasta'):
		RECORD = record.id
		FASTA = str(record.seq)

		#FOO = re.split(":|-",RECORD)
		FOO = re.split(":",RECORD)
		CONTIG = FOO[0]
			
		CAZYMES = '_'.join(FULL_INFO_DICT[CONTIG])
		CAZYMES = re.sub(" ","",CAZYMES)

		if int(FORE) > 0:
			NEW_RECORD = RECORD+"_"+CAZYMES+"_incl_"+FORE+"_nt_upstream"
		else:
			NEW_RECORD = RECORD+" "+FULL_INFO_DICT[CONTIG]

		output.write(">"+NEW_RECORD+"\n"+FASTA+"\n")

	output.close()

	if int(FORE) > 0:
		os.system(' '.join([
			"mv",
			OUTPUT+"/"+NAME+".recovered.foo",
			OUTPUT+"/"+NAME+".recovered."+FORE+".fa"
		]))
	else:
		os.system(' '.join([
			"mv",
			OUTPUT+"/"+NAME+".recovered.foo",
			OUTPUT+"/"+NAME+".recovered.fa"
		]))

	try:
		ZIP_ME
		GZ(FASTA_FILE)
	except:
		print "NOTHING TO ZIP"

	
## REMOVE TEMPORARY FILES
os.system("rm -r TEMP")
os.system("rm -r all.hmm.ps.len")
os.system("rm -r hmmscan-parser.sh")
os.system("rm "+OUTPUT+"/*length.txt")
os.system("rm hmmscan-parser.sh")
os.system("rm all.hmm.ps.len")
os.system("rm "+OUTPUT+"/"+"*.bed")
os.system("rm  find . -size 0 -delete") #Remove dud bedtools recovery attempts

	
end = timeit.default_timer()
print end - start	
