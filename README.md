Fastq-submission
================

Parses FASTQ reads into individual files, organized by sample name. Merged FASTA file is provided as input to define sample:header relationships.

-- INPUT FILE FORMATS --
sample_parser.py takes a minimum of two files as input:

	1) a .FNA format file with Illumina sequence identifiers, delineated by ':' 
		
			>[sample id]_[sequence index] [instrument name]:[run id]:[flowcell id]:[flowcell lane]:[tile #]:[xpos]:[ypos] 
			
			e.g.
				>A3Nc_0 M00704:71:000000000-A5D7M:1:1101:17028:1428 1:N:0:0
	
	2 a FASTQ format file with Illumina sequence identifiers
			
			e.g.
				@M00704:71:000000000-A5D7M:1:1101:16105:1358 1:N:0:0
			

-- USAGE -- 
sample_parser.py takes two arguments:

	>$ python sample_parser.py [sample_definition_file] [read_file(s)]
	
		1) the file name of .FNA files defining relationships between sequence identifiers and sample IDs
		2) a comma-delimited (no spaces) list of one or more FASTQ files to parse
	
	e.g.
	>$ python sample_parser.py seqs.fna R1.fastq,R2.fastq
	
	
-- OUTPUT --
sample_parser.py generates one .fastq file per sample ID per .fastq file. The data within these files follows 
the format for a multiline FASTQ file. 

Output file nomenclature is of the form:
	
		[SAMPLE ID]_[READ FILE].fastq.gz
		
		Where 'SAMPLE ID' is the ID that is related to the header of all the records in that file and 'READ FILE' is the name of the .fastq file that sequences were obtained from.

Any records found in the read files that lack a sample ID definition are sent to an [READ FILE]_ERR.fastq file. 

All files are written to disk using the GZIP utility and, thus, are appended with the extension '.gz'

