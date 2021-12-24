# miRAW Wrapper

##  Introduction
This Python script provides a series of pre- and post- processing functions to simplify usage of the **miRAW** target prediction tool. Pre-processing involves setting up shell scripts for bulk target prediction tasks. Post-processing involves filtering target prediction results according to different parameters

## Pre-processing & Setting up a miRAW run
The basic requirements for setting up a miRAW run are:
 
* a list of miRNAs
* a list of 3'UTR targets
* a candidate site selection model

In addition, paths to the miRAW jar file and trained model need to be specified.

The complete set of parameters that can be specified are as follows:

### Experiment Name
`"-e", "--exptName"`
This is the name given to the script, to make it easier to identify.

### Deep Learning model
`-d/--dlModel`
See the miRAW README.md on how to build your own model. But miRAW comes with one trained model for human miRNAs.

### Out Folder
`"-o", "--outFolder"`
the full path to where the results will be written. This allows the user to generate scripts to execute **miRAW** on a remote machine.

### Candidate Site Selection Model (CSSM)
`"-c", "--cssm"`
the type of candidate site model used in the prediction. This can be set to this can be **Regular**, **Pita**, **targetScan** or **Personalized**

### UnifiedFile Path
`-u/--unifiedFile`
miRAW takes a single input file that contains both the query **miRNAs** and **3'UTR** sequences. You can specify the path to an existing file or you can generate one by specifying both a miRNA file `-m/--miRFile` and a 3'UTR file `-3/--3UTR`. Both these files should be in **fasta** format

### miRNA Fasta File
`"-m", "--mirnaFile"`
the full path to a file containing the **miRNA** sequences in FASTA format.

### 3'UTR Fasta File
`"-3", "--3UTR"`
the full path to a file containing the **3'UTR** sequences in FASTA format.

### How to write the results
`"-t", "--split"`
For large numbers of miRNAs and 3'UTR sequences, the generated results can be difficult to parse out. Therefore, they can be broken down and written out in different ways.

- `-t/--split=0`
: use existing unified File (specified with `-u/--unifiedFile`)
- `-t/--split=1`
: write by miRNA (single file for all predictions) 
- `-t/--split=2`
: write by 3'UTR (single file for all predictions)
- `-t/--split=3`
: split by miRNA (separate file for each miRNA)
- `-t/--split=4`
: split by 3'UTR (separate file for each 3'UTR)



### jar file location
`"-j", "--jarLoc"`
the full path to where the **miRAW.jar** executable is stored on the local machine. The jar file is packaged with the script. This uses more disk space, but it simplifies the script building process.


miRAW uses a sliding window to analyze a 3'UTR, so the default values will use a sliding window of 40nt with a step size of 5nt. Generally, there isn't much advantage in changing these values.

### maximum site length
`"-W", "--maximum_site_length"`
target site size (in nucleotides). default 40nt 

### seed alignment offset
`"-S", "--seed_alignment_offset"`
step_size (in nucleotides). default 5nt.

### energy filtering
`"-E", "--energy_filtering"`
for conservative models (i.e. pita and targetScan) performance improves in the absence of energy filtering. for miRAW models using extended seed regions, performance improves when energy filtering is used. This flag allows you turn to turn on / off energy filtering. Allowed values are (`True`/`False`).



### Example parameters
```
-e 186712beta 
-d /Users/simonray/DropboxUiO/dropData/miraw/bestModel.bin 
-j /Users/simonray/DropboxUiO/dropData/miraw/miRAW.jar
-o /Users/simonray/data/186712beta
--cssm Personalised;10:1:7
-m /Users/simonray/DropboxUiO/dropData/mirbase/22.1/hsa_22pt1uniq.fa
-3 /Users/simonray/DropboxUiO/dropData/ensembl/ensembl_hsa_3utrs.uniq.fa
--split=4
--energy_filtering=True
```


### Useful commands

The following will flatten a FASTA file so that each sequence is only on one line
```
 awk '/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < ensembl_hsa_3utrs.uniq.fa > ensembl_hsa_3utrs.uniq_flat.fa
```


```
grep -w -A 1 -f test.txt ensembl_hsa_3utrs.uniq_flat.fa |grep -v -- "^--$"
```


