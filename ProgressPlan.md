# CNV array pipeline

## Aims

- get the snakemake pipeline publishable until summer break
	- fix all current issues
	- decide on parts of whishlist to incorporate:
		- moasacism? (if possible without to much effort)
		- SNP distance between samples
	- show that it works // test expectations:
		- comparison with WGS calls?
		- comparison with expert knowledge
		- training adapting parameter settings from specific training set, evaluate on other samples
		- cellline relationships: later lines should have more validated calls, not lose any(?)
		- (potentially ?) include tests with external data as well
		- include curated data ressource for tool validation
		
- start parallel developement of interactive shiny tool
	- interactivity for annotation
	- keep it simple for now

- long term outlook / project continuation
	2) Meta analysis: collect data from multiple cores & look for patterns of CNVs (&SNPs?) with metadata features
	3) Publicly usable tool: build interface and open up our pipeline to outside users

## Timeline

~ April: 	coding 
~ May:		test expectations
~ June:		writing paper
~ July:		buffer


## Known Issues / Bugs:

- tables for individual calls miss reportable & reference annotation
	-> annotation needs to be done also for the 'pre-overlap' calls which are plotted OR they need inherit annotation somehow (better, but harder)
- There may be some calls without probes? / position mismatches?
- column selection in report tables & Chr should be factor

- multi-overlaps lead to probelms
	-> need better merging &/ tool-overlap


## Current state / TODO

- preprocessing & report need gene track overlay
	- use a bed file? where to get it from?
- SNP clustering step (only run if enabled in config?)
	- add SNP distance to reference into QC overview table
	
? unclear prio -> maybe needed for good output ?
- merge by probe distance
- improve tool overlap (need at least tools for testing?)
- ref (& tool?) overlap with min 80% (reciprocal?) overlap
- define hotspots / regions of interest
- mosaicism (R-GADA / others)
- individual tool stats
- hotspot regions / regions of interest
- circO plots
- installation procedure: as much as possible should be in bioconda

- decide on a proper set of static data files (fasta, gtf/gff, bpm, csv, egt, pfb, GCmodel)
	-> GCmodel should be calculated from fasta somehow?
	-> pfb calculation from egt->vcf seems to work well
	-> bpm | csv might be optional (need one but maybe not both) 

Lower Priority goals:
- PennCNV trio calling
- make certain tools optional
- DT import to HTML
- outout of CNV vcf // IGV readble files
- ensure rerun (of report? of everything?) on parameter change in config ?!
- other filter options?
- additional / more modern CNV callers
- machine learning approach for CNV annotation
- addition of bed files (i.e. USCS tracks) on top of gene overlap (refseq enchancers, ...)

