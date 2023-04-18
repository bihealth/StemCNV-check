# CNV array pipeline

## Aims

- get the snakemake pipeline publishable until summer break
	- fix all current issues
	- decide on parts of whishlist to incorporate:
		- moasacism? (if possible without to much effort)
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

- None

## Current state / TODO

- define hotspots / regions of interest

? unclear prio -> maybe needed for good output ?
- merge by probe distance
- ref (& tool?) overlap with min 80% (reciprocal?) overlap
- individual tool stats
- circO plots
- installation procedure: as much as possible should be in bioconda

- decide on a proper set of static data files (fasta, gtf/gff, bpm, csv, egt, pfb, GCmodel)
	-> GCmodel should be calculated from fasta somehow?
	-> pfb calculation from egt->vcf seems to work well
	-> bpm | csv might be optional (need one but maybe not both) 

Lower Priority goals:
- PennCNV trio calling
- output of CNV vcf // IGV readable files
- ensure rerun (of report? of everything?) on parameter change in config ?!
- other filter options?
- additional / more modern CNV callers
- machine learning approach for CNV annotation
- addition of bed files (i.e. USCS tracks) on top of gene overlap (refseq enchancers, ...)

