# Locus-Splitter
## Split target-enrichment sequence data into its partitions for alignment and downstream analyses

### Basic use case
Sometimes a target-enrichment protocol will involve baits tiled across a conserved region so that the sequence reads can extend into some variable regions. This often involves one or more exons but not always (some UCEs or CNEEs). For the cases where the data can be neatly described: a conserved core (c), variable left flank (l), and variable right flank (r), I wanted to separate these regions as well as possible. This is in part to maybe help alignment, but there might be analyses where we do not want the cores or *vice versa*.

<div align="center">

<p>llllllllllllllllccccccccccccccccrrrrrrrrrrrrrrrr<br>
llllllllllllllllccccccccccccccccrrrrrrrrrrrrrrrr<br>
llllllllllllllllccccccccccccccccrrrrrrrrrrrrrrrr<br>
llllllllllllllllccccccccccccccccrrrrrrrrrrrrrrrr</p>

</div>

Starting with fasta files, the below script will split the l, c, and r regions into separate fasta files. This is accomplished by running the `splitLoci.pl` script twice.
```
perl splitLoci.pl --controlFile splitLoci.ctl --template split-template.sh --runmode 0
perl splitLoci.pl --controlFile splitLoci.ctl --template split-template.sh --runmode 1
```
When `runmode = 0`, some blast output is generated for each individual against a reference database that should be the c regions from related individuals. This means if the c region is an exon, the reference should be the whole exons, not just bits of the sequence used to synthesize the baits. Running the script again as `runmode = 1` uses the blast coordinates to split the sequences.

You probably want to do something with the sequences now. I align the l, c, and r regions separately. I also translate the c to get amino acid alignments and then back-translate to keep the codons in frame. For closely-related species, this is not really necessary. It might be helpful for some larger studies where the exon-intron boundaries can get a little messy though.
```
perl splitLoci.pl --controlFile splitLoci.ctl --template split-template.sh --runmode 2
```

### The control file

User options are set in the control file - splitLoci.ctl
```
ROOT\_DIR = YOUR\_PATH/SplitLoci/
FASTA\_INPUT = YOUR\_PATH/SplitLoci/inputFastaFiles
BLAST\_OUTPUT = YOUR\_PATH/SplitLoci/blastOutput
FASTA\_OUTPUT = YOUR\_PATH/SplitLoci/splitFastaFiles
SUMMARYSTATS\_OUT = YOUR\_PATH/SplitLoci/splitStats
INDIVIDUAL\_LIST = YOUR\_PATH/SplitLoci/speciesList.txt
REFERENCE\_FASTA = YOUR\_PATH/SplitLoci/referenceFastaFiles
ALIGNMENT\_IN = YOUR\_PATH/SplitLoci/alignedRegions
ALIGNMENT\_OUT = YOUR\_PATH/SplitLoci/alignedLoci
SCHEDULER = sbatch

BLAST = tblastx 
DBTYPE = nucl
EVALUE = 1e-2
DELIMITER = _
```

A lot of the options are just directories so the script knows where to find and put files. About the input directories:
* FASTA\_INPUT = all of your fasta files should be here. They should have the \*.fa or \*.fasta ending
* ALIGNMENT\_IN = where the *aligned* fasta files can be found. This will be something you do on your own after `runmode = 2`. The ALIGNMENT\_IN directory should have three sub-directories that contain the sequences implied:
	* leftFlank
	* core
	* rightFlank
* INDIVIDUAL\_LIST = a simple list of individuals or a subset included in the fasta files. All that matters is the first thing per line is the taxon name. It could look something like this:
```
speciesA	ID0001	-49.85	20.39	my favorite species
speciesB	ID0032 	-48.7	19.51	DNA from voucher
.
.
.
speciesN	ID9408	-47.8	21.55	outgroup
```

Some other options:
* Scheduler = relevant to the template file described below. If on a cluster, this would be the command for scheduling a job (sbatch in the case of slurm). You can still use this on a local desktop! Just change "sbatch" to "bash" and you will execute basic bash scripts.
* BLAST = the blast program you want to use. I use tblastx to go from nucleotide sequences to a nucleotide reference where **the reference must be in-frame starting at position 1**. blastn would be another option for nucleotide-to-nucleotide. This could just as easily be blastx if you only have the amino acids for the reference sequences.
* DBTYPE = necessary for configuring the blast database of reference sequence. *nucl* for nucleotide or *prot* for amino acid.
* EVALUE = The evalue threshold for if a hit is output or not. I keep this low because some core regions that are very short will have high evalues. Maybe having one at all is a bad idea, but maybe it helps swat away some garbage too?
* DELIMITER = a charater used to split extraneous information from a taxon label so it matches stuff in the INDIVIDUAL\_LIST. My motivation for this is that I append an identifier like __1  and __2 to taxon labels to keep track of individual alleles. Setting `DELIMITER = __` saves me from having to enter every individual allele into the list and bins the blast jobs to cut down on I/O.  

### The template file

Environment variables and cluster directives (if applicable) can be set in the template file - split-template.sh
```
#!/bin/bash
#SBATCH --job-name=__RUNID__
#SBATCH --output=__LOGFILE__
#SBATCH --mail-user=YOUR_EMAIL
#SBATCH --mail-type=FAIL
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=2000M
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=YOUR_PARTITION_OR_DELETE
module load blast/2.10.0
```

Here is have the bits necessary to make slurm launch a single blast job per individual. This might be helpful if you have many individuals or a lot of data. But it is possible to use the template to make things run locally too!

```
#!/bin/bash
export PATH=PATH_TO_BLAST_PROGRAMS:$PATH
```

And by setting `scheduler = bash`, every blast job will run serially on my desktop or laptop with some basic UNIX environment. In this case, blast is compile somewhere else, but I can make it available by exporting the path at runtime. I like executing things through a template this way so nobody is forced to do system-wide installs.