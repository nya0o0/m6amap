## Goal of the project: Is it clear to you from the proposal.md how the goal can be accomplished using Python and the specified packages?

The goal is clear to me from the proposal.md. The project's aim is to integrate RNA modification prediction from m6Anet with biological annotation from M6ADD and other resourses like GO and KEGG into a single pipeline. Python is a suitable languate for this task, as it enables to write functions and pack them in a package cleanly. The use of command-line tools, APIs, and interactive plotting libraries looks to work well.

The goal is also unique as described that it is trying to connect predicted m6A sites to biological meaning, which has been not integrated into a single tool.

## The Data: Is it clear to you from the proposal.md what the data for this project is, or will look like?

From the proposal.md, the data imput requirements are:
- BAM file (aligned reads)
- Transcriptome reference (FA)
- FASTA read file
- GTF gene annotation file

External resources include:
- GENCODE annotations
- M6ADD for annotation
- Gene Oncology and KEGG pathway databases (APIs)

The output data will look like a CSV table with full annotations and an interactive dashbord.

The output file's columns will be:
- transcript_id
- gene_id
- transcript_position
- n_reads	
- probability_modified	
- kmer	
- mod_ratio	
- sample_id	
- group_id	
- genome_pos	
- chromosome	
- dist_up_exon_junc	
- dist_down_exon_junc	
- region	
- gene_name	
- transcript_type	
- Associated Disease (from M6ADD)	
- Pathways

## The code: (Look at the Python code files in detail first and try to comprehend a bit of what is written so far)

The current codebase includes a working setup structure (setup.py, src/, installable with pip install -e .). Code does not seem to be written yet. From the proposal.md, however, it's clear what the major components should be.

## Does the current code include a proper skeleton (pseudocode) for starting this project?
What can this code do so far?
Given the project description, what are some individual functions that could be written to accomplish parts of this goal?
Code contributions/ideas: See description below for what to write here.

It seems like the implementation logic is still in development.

potential functions
- run_m6anet(BAM): Run m6Anet on the given input and output a raw CSV
- parse_m6anet(filepath): Load and process m6Anet CSV
- annotate_w_m6add(df): Cross-reference M6ADD to annotate modification sites with disease and pathway info.
- query_go(gene_id): Get GO terms for a gene
- query_kegg(gene_name): Get pathway annotations for a gene
- generate_network_graph(df): Use plotly to build a disease-gene-pathway network graph as discribed
- launch_dashbord(df): Let users interact with the data visually




