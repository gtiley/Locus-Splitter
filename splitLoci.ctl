#input and output paths
ROOT_DIR = YOUR_PATH/SplitLoci/
FASTA_INPUT = YOUR_PATH/SplitLoci/inputFastaFiles
BLAST_OUTPUT = YOUR_PATH/SplitLoci/blastOutput
FASTA_OUTPUT = YOUR_PATH/SplitLoci/splitFastaFiles
SUMMARYSTATS_OUT = YOUR_PATH/SplitLoci/splitStats
INDIVIDUAL_LIST = YOUR_PATH/SplitLoci/speciesList.txt
REFERENCE_FASTA = YOUR_PATH/SplitLoci/referenceFastaFiles
ALIGNMENT_IN = YOUR_PATH/SplitLoci/alignedRegions
ALIGNMENT_OUT = YOUR_PATH/SplitLoci/alignedLoci
SCHEDULER = sbatch

#paths to software or just the binaries if already in your path
BLAST = tblastx 
DBTYPE = nucl
EVALUE = 1e-2
DELIMITER = _
