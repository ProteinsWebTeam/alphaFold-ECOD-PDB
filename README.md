# alphaFold-ECOD-PDB
Workflow to compare AlphaFold models with PDB structures and ECOD domains

## Requirements

The following tools should be installed in the working directory:
- Foldseek `$ wget --no-check-certificate https://mmseqs.com/foldseek/foldseek-linux-sse41.tar.gz; tar xvzf foldseek-linux-sse41.tar.gz; mv foldseek/bin alphaFold-ECOD-PDB/; rm -r foldseek; rm foldseek-linux-sse41.tar.gz;`
- TMalign `$ cd alphaFold-ECOD-PDB; wget https://zhanggroup.org/TM-align/TMalign.gz; gunzip TMalign.gz; chmod u+x TMalign; rm TMalign.gz`

## Run the comparison

Fill the config file (`config_comp.ini`). 
- *ecod_pdb_file* should be `<path>/ecod.latest.F70.pdb.tar.gz`
- *ecod_file* and matches_file should have `.txt` extension.
- *output_file* should have `.m8` extension
- *af2pfam_file* should have `.tsv` extension

Execute `$ ./pipeline.sh config_comp.ini`

The pipeline is divided in 4 steps:

### Find the AlphaFold models 

This is executed by the python script `find_pfam_duf.py`.
Only the AlphaFold models matching Pfam without PDB structures are used for the comparison.

### Download files from ECOD

2 files need to be downloaded:
- The compressed file containing the PDB files (http://prodata.swmed.edu/ecod/distributions/ecod.latest.F70.pdb.tar.gz). The PDB files are organised in multiple subdirectories, the script put them in one unique directory (`ecod_pdb`)
- The ECOD domain file (http://prodata.swmed.edu/ecod/distributions/ecod.latest.domains.txt). It contains the PDB/ECOD domains mapping.

### Run Foldseek

Run Foldseek using the AlphaFold pdb files and ECOD pdb files downloaded.

### Extract relevant matches

This is executed by the python script `find_relevant_matches.py`.
From the output of Foldseek, we only consider matches with an e-value < e-05.
Additionnally, TMalign is run and only matches with TMscores > 0.6 for normalized by length of Chain_2 are kept. 
Annotations to Pfam domains and ECOD domains are also provided when available.
