#!/usr/bin/perl -w
use strict;

#----------------------------------------------------------------------------------------#
#George P. Tiley
#17 April 2022
#contact: g.tiley@kew.org
#Split loci from phased fasta output into their left intron, exon, and right intron
#reuses the control file and template for PATE since the blast can be distributed
#----------------------------------------------------------------------------------------#

# Accepted AND necessary commands
my @checkArgs = ("controlFile","template","runmode");
my %passedArgs = ();
if (scalar(@ARGV) == 0)
{
die "/*--------INPUT PARAMETERS--------*/\n
--controlFile STRING <control file>
--template STRING <template file for distributing on a cluster>
--runmode INT <0 = run blast || 1 = split fasta files after blast is complete || 2 = concatenate locus parts back to original locus structure from a seperate directory (e.g. after alignment)>

\n/*--------EXAMPLE COMMAND--------*/\n
    perl splitLoci.pl --controlFile splitLoci.ctl --template template.sh --runmode 0\n
	
\n/*--------NOTES--------*/\n
This script assumes you want to split the loci into a left, middle, and right chunk. This is approriate for some types of target enrichment data, but not others that might be tiled across multiple exons.\n";
}

elsif (scalar(@ARGV) > 0)
{
    for my $i (0..(scalar(@ARGV) - 1))
    {
		if ($ARGV[$i] eq "--controlFile")
		{
	    	$passedArgs{controlFile} = $ARGV[$i+1];
		}
		if ($ARGV[$i] eq "--runmode")
		{
	    	$passedArgs{runmode} = $ARGV[$i+1];
	    	if (($passedArgs{runmode} != 0) && ($passedArgs{runmode} != 1) && ($passedArgs{runmode} != 2))
	    	{
	    		die ("--runmode must be either 0 or 1 or 2\n")
	    	}
		}
		if ($ARGV[$i] eq "--template")
		{
        	$passedArgs{template} = $ARGV[$i+1];
		}
	}
	foreach my $arg (@checkArgs)
	{
		if (! exists $passedArgs{$arg})
		{
	    	die "/*--------MISSING PARAMETER--------*/\nMissing command line argument: $arg\n\n";
		}
	}
}

my %controlArgs = ();
open FH1,'<',"$passedArgs{controlFile}";
while (<FH1>)
{
	if (/ROOT_DIR\s+\=\s+(\S+)/)
    {
		$controlArgs{ROOT_DIR} = $1;
    }
	if (/FASTA_INPUT\s+\=\s+(\S+)/)
    {
		$controlArgs{FASTA_INPUT} = $1;
    }
    if (/BLAST_OUTPUT\s+\=\s+(\S+)/)
    {
		$controlArgs{BLAST_OUTPUT} = $1;
		if ($passedArgs{runmode} == 0)
		{
			system "mkdir $controlArgs{BLAST_OUTPUT}";
		}
    }
    if (/FASTA_OUTPUT\s+\=\s+(\S+)/)
    {
		$controlArgs{FASTA_OUTPUT} = $1;
		if ($passedArgs{runmode} == 1)
		{
			system "mkdir $controlArgs{FASTA_OUTPUT}";
			system "mkdir $controlArgs{FASTA_OUTPUT}/leftFlank";
			system "mkdir $controlArgs{FASTA_OUTPUT}/rightFlank";
			system "mkdir $controlArgs{FASTA_OUTPUT}/core";
		}
    }
    if (/SUMMARYSTATS_OUT\s+\=\s+(\S+)/)
    {
		$controlArgs{SUMMARYSTATS_OUT} = $1;
		if ($passedArgs{runmode} == 1)
		{
			system "mkdir $controlArgs{SUMMARYSTATS_OUT}";
		}
    }
    if (/INDIVIDUAL_LIST\s+\=\s+(\S+)/)
    {
        $controlArgs{INDIVIDUAL_LIST} = $1;
    }
    if (/REFERENCE_FASTA\s+\=\s+(\S+)/)
    {
        $controlArgs{REFERENCE_FASTA} = $1;
    }
    if (/ALIGNMENT_IN\s+\=\s+(\S+)/)
    {
        $controlArgs{ALIGNMENT_IN} = $1;
    }
    if (/ALIGNMENT_OUT\s+\=\s+(\S+)/)
    {
        $controlArgs{ALIGNMENT_OUT} = $1;
        if ($passedArgs{runmode} == 2)
		{
			system "mkdir $controlArgs{ALIGNMENT_OUT}";
			system "mkdir $controlArgs{ALIGNMENT_OUT}/loci";
			system "mkdir $controlArgs{ALIGNMENT_OUT}/concatenated";
		}
    }
    if (/SCHEDULER\s+\=\s+(\S+)/)
    {
        $controlArgs{SCHEDULER} = $1;
    }
    if (/BLAST\s+\=\s+(\S+)/)
    {
        $controlArgs{BLAST} = $1;
    }
    if (/DBTYPE\s+\=\s+(\S+)/)
    {
        $controlArgs{DBTYPE} = $1;
    }
    if (/EVALUE\s+\=\s+(\S+)/)
    {
        $controlArgs{EVALUE} = $1;
    }
    if (/DELIMITER\s+\=\s+(\S+)/)
    {
        $controlArgs{DELIMITER} = $1;
    }
}
close FH1;

if ($passedArgs{runmode} == 0)
{
####
#Create the necessary blast databases from the reference sequences
	system "mkdir $controlArgs{BLAST_OUTPUT}/db";
	open OUT1,'>',"$controlArgs{BLAST_OUTPUT}/db/makeblastdb.sh";
	open FH1,'<',"$passedArgs{template}";
	while (<FH1>)
	{
		my $line = $_;
		chomp $line;
		$line =~ s/__RUNID__/db/;
		$line =~ s/__LOGFILE__/db.log/;
		print OUT1 "$line\n";
	}
	close FH1;
	print OUT1 "cd $controlArgs{BLAST_OUTPUT}/db\n\n";
	
	my @refSeqFastas = glob("$controlArgs{REFERENCE_FASTA}/*.fa*");
	foreach my $rsf (@refSeqFastas)
	{
		my @temp = ();
		@temp = split(/\//,$rsf);
		if ($temp[scalar(@temp)-1] =~ m/(\S+)\.fa.*/)
		{
			my $prefix = $1;
			print OUT1 "makeblastdb -in $rsf -out $prefix -parse_seqids -dbtype $controlArgs{DBTYPE}\n";
		}
		else
		{
			die "Could not split fasta file name when building databases:\n $rsf\n\n. Check that there are no unexpected symbols, which will be important for downstream things too.\n";
		}
	}
	close OUT1;
	#Alternate strategy to set dependency but then limits to slurm
	#my $dbid = `$controlArgs{SCHEDULER} $controlArgs{BLAST_OUTPUT}/db/makeblastdb.sh`;
	#chop $dbid
    system "$controlArgs{SCHEDULER} $controlArgs{BLAST_OUTPUT}/db/makeblastdb.sh";
	#simple solution, just wait 2 minutes for db to process. Probably some scenarios where the db gets stuck or takes a long time though.
	sleep(120);
	
########
#Go through each individual, make a single fasta for each, blast and return the best hit and alignment
########

	my %indList = ();
	my %indArrays = ();
	my %indChecker = ();
	open FH1,'<',"$controlArgs{INDIVIDUAL_LIST}";
	while (<FH1>)
	{
		if (/^(\S+)/)
		{
			my $ind = $1;
			$indList{$ind} = 1;
			system "mkdir $controlArgs{BLAST_OUTPUT}/$ind";
		}
	}
	close FH1;

	my %seqs = ();
	my %loci = ();
	my @inputFastas = glob("$controlArgs{FASTA_INPUT}/*.fa*");
	foreach my $ipfa (@inputFastas)
	{
		my @temp = ();
		@temp = split(/\//,$ipfa);
		if ($temp[scalar(@temp)-1] =~ m/(\S+)\.fa.*/)
		{
			my $prefix = $1;
			if (! exists $loci{$prefix})
			{
				$loci{$prefix} = 1;
			}
			
			my $tax = "";
			my $seq = "";
			open FH1,'<',"$ipfa";
			while (<FH1>)
			{
				if (/^>(\S+)/)
				{
					$tax = $1;
					$seqs{$tax}{$prefix} = "";
					if (exists $indList{$tax})
					{
						if (! exists $indChecker{$tax})
						{
							push @{$indArrays{$tax}}, $tax;
							$indChecker{$tax} = 1;
						}
					}
					elsif (! exists $indList{$tax})
					{
						my @checkID = ();
						@checkID = split(/__/,$tax);
#						print "$prefix\t$tax\t$checkID[0]\n";
						if (exists $indList{$checkID[0]})
						{
							if (! exists $indChecker{$tax})
							{
								push @{$indArrays{$checkID[0]}}, $tax;
								$indChecker{$tax} = 1;
							}
						}
						elsif (! exists $indList{$checkID[0]})
						{
							die "Problem - Could not assign $tax to an individul in the list!\n";
						}
					}
				}
				elsif (/(\S+)/)
				{
					$seq = $1;
					$seqs{$tax}{$prefix} = $seqs{$tax}{$prefix} . $seq;
				}
			}
		}
		else
		{
			die "Could not split fasta file name when reading input sequences:\n $ipfa\n\n. Check that there are no unexpected symbols, which will be important for downstream things too.\n";
		}
	}
	
	foreach my $ind (sort keys %indList)
	{
########
#blast jobs are executed through 1 script per individual to cut down on I/O
########	
		open OUT1,'>',"$controlArgs{BLAST_OUTPUT}/$ind/blast.sh";
		open FH1,'<',"$passedArgs{template}";
		while (<FH1>)
		{
			my $line = $_;
			chomp $line;
			$line =~ s/__RUNID__/$ind.blast/;
			$line =~ s/__LOGFILE__/$ind.blast.log/;
			print OUT1 "$line\n";
		}
		close FH1;
		print OUT1 "cd $controlArgs{BLAST_OUTPUT}/$ind\n\n";
		foreach my $locus (sort keys %loci)
		{
#			print "$ind\t$locus\t$indArrays{$ind}[0]\n";
			my $checkforseqs = 0;
			foreach my $tax (@{$indArrays{$ind}})
			{
				if (exists $seqs{$tax}{$locus})
				{
					$checkforseqs = 1;
				}
			}
			if ($checkforseqs == 1)
			{
#				print " ==> It exists $ind\t$locus\t$indArrays{$ind}[0]\n";
				open OUT2,'>',"$controlArgs{BLAST_OUTPUT}/$ind/$locus.fasta";
				foreach my $tax (@{$indArrays{$ind}})
				{
					if (exists $seqs{$tax}{$locus})
					{
						print OUT2 ">$tax\n$seqs{$tax}{$locus}\n";
					}
				}
				close OUT2;
########
#Special handling for hybpiper output
#It like to smush _supercontig behind the locus name
#So remove for the database since it is likely without the appendix
				my $dblocusname = $locus;
				$dblocusname =~ s/\_supercontig//;				
########				
				print OUT1 "$controlArgs{BLAST} -query $controlArgs{BLAST_OUTPUT}/$ind/$locus.fasta -db $controlArgs{BLAST_OUTPUT}/db/$dblocusname -max_target_seqs 5 -max_hsps 100 -evalue $controlArgs{EVALUE} -outfmt \"6 qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe\" -out $controlArgs{BLAST_OUTPUT}/$ind/$locus.blast.out\n";
			}
		}
		close OUT1;
		system "$controlArgs{SCHEDULER} $controlArgs{BLAST_OUTPUT}/$ind/blast.sh";
	}
}
	
	
####

elsif ($passedArgs{runmode} == 1)
{
########
#First get all sequences then split by locus based on blast coordinates
########

	my %indList = ();
	my %indArrays = ();
	my %indChecker = ();
	my %splitStats = ();
	open FH1,'<',"$controlArgs{INDIVIDUAL_LIST}";
	while (<FH1>)
	{
		if (/^(\S+)/)
		{
			my $ind = $1;
			$indList{$ind} = 1;
		}
	}
	close FH1;
	
	my %seqs = ();
	my %loci = ();
	my @inputFastas = glob("$controlArgs{FASTA_INPUT}/*.fa*");
	foreach my $ipfa (@inputFastas)
	{
		my @temp = ();
		@temp = split(/\//,$ipfa);
		if ($temp[scalar(@temp)-1] =~ m/(\S+)\.fa.*/)
		{
			my $prefix = $1;
			if (! exists $loci{$prefix})
			{
				$loci{$prefix} = 1;
			}
			
			my $tax = "";
			my $seq = "";
			open FH1,'<',"$ipfa";
			while (<FH1>)
			{
				if (/^>(\S+)/)
				{
					$tax = $1;
					$seqs{$tax}{$prefix} = "";
					if (exists $indList{$tax})
					{
						if (! exists $indChecker{$tax})
						{
							push @{$indArrays{$tax}}, $tax;
							$indChecker{$tax} = 1;
						}
					}
					elsif (! exists $indList{$tax})
					{
						my @checkID = ();
						@checkID = split(/__/,$tax);
						if (exists $indList{$checkID[0]})
						{
							if (! exists $indChecker{$tax})
							{
								push @{$indArrays{$checkID[0]}}, $tax;
								$indChecker{$tax} = 1;
							}
						}
						elsif (! exists $indList{$checkID[0]})
						{
							die "Problem - Could not assign $tax to an individul in the list!\n";
						}
					}
				}
				elsif (/(\S+)/)
				{
					$seq = $1;
					$seqs{$tax}{$prefix} = $seqs{$tax}{$prefix} . $seq;
				}
			}
		}
		else
		{
			die "Could not split fasta file name when reading input sequences:\n $ipfa\n\n. Check that there are no unexpected symbols, which will be important for downstream things too.\n";
		}
	}

########
# Store the flanking coordinates for each seq to memory
########
	
	my %coordinates = ();
	foreach my $ind (sort keys %indList)
	{
		foreach my $locus (sort keys %loci)
		{
			my $checkforseqs = 0;
			foreach my $tax (@{$indArrays{$ind}})
			{
				if (exists $seqs{$tax}{$locus})
				{
					$checkforseqs = 1;
				}
			}
			if ($checkforseqs == 1)
			{
				open FH1,'<',"$controlArgs{BLAST_OUTPUT}/$ind/$locus.blast.out";
				my %maxbitscores = ();
				my %isinframe = ();
				my %foundhit = ();
				
				while (<FH1>)
				{
					my $line = $_;
					chomp $line;
					##Query_Sequence	Reference_Sequence	PercentMatch	MatchLength	Mismatches	Gaps	QueryStart	QueryEnd	ReferenceStart	ReferenceEnd	Evalue	BitScore
					##MSV2176__REF	BPKH_Neurachne_annularis__REF	94.545	55	3	0	266	102	167	3	1.64e-33	130
					if ($line =~ m/(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)/)
					{
						my $querySeq = $1;
						my $qstart = $2;
						my $qend = $3;
						my $rstart = $4;
						my $rend = $5;
						my $bitscore = $6;
						my $rframe = $7;
						
						if ($controlArgs{BLAST} eq "tblastx")
						{
							if (! exists $maxbitscores{$querySeq})
							{
								$maxbitscores{$querySeq} = $bitscore;
								$foundhit{$querySeq} = 0;
								$isinframe{$querySeq} = 0;
								if ($rframe == 1)
								{
									$isinframe{$querySeq} = 1;
								}
							}
							elsif (exists $maxbitscores{$querySeq} && $isinframe{$querySeq} == 0)
							{
								if (($maxbitscores{$querySeq} - $bitscore) < (0.5 * $maxbitscores{$querySeq}))
								{
									if ($rframe == 1)
									{
										$isinframe{$querySeq} = 1;
									}
								}
							}
						}
						
						elsif ($controlArgs{BLAST} ne "tblastx")
						{
							if (! exists $maxbitscores{$querySeq})
							{
								$maxbitscores{$querySeq} = $bitscore;
								$isinframe{$querySeq} = 1;
								$foundhit{$querySeq} = 0;
							}
						}
						
						if ($isinframe{$querySeq} == 1 && $foundhit{$querySeq} == 0)
						{
						
							my $dquery = $qend - $qstart;
							my $dref = $rend - $rstart;
						
							if (($dquery < 0 && $dref > 0) || ($dquery > 0 && $dref < 0))
							{
								$seqs{$querySeq}{$locus} = Complement($seqs{$querySeq}{$locus});
								if (($dquery < 0 && $dref > 0))
								{
									$coordinates{$querySeq}{$locus}{start} = $qend;
									$coordinates{$querySeq}{$locus}{end} = $qstart;
								}
								else
								{
									$coordinates{$querySeq}{$locus}{start} = $qstart;
									$coordinates{$querySeq}{$locus}{end} = $qend;
								}
							}
							elsif (($dquery < 0 && $dref < 0) || ($dquery > 0 && $dref > 0))
							{
								if ($dquery < 0 && $dref < 0)
								{
#									print "//\nCheck that $querySeq $locus looks correct\n//\n\n";
									$coordinates{$querySeq}{$locus}{start} = $qend;
									$coordinates{$querySeq}{$locus}{end} = $qstart;
								}
								else
								{
									$coordinates{$querySeq}{$locus}{start} = $qstart;
									$coordinates{$querySeq}{$locus}{end} = $qend;
								}
							}
							elsif ($dquery == 0 || $dref == 0)
							{
								print "//--------WARNING!--------//\n$querySeq for locus $locus returned 0 length blast hit against best reference\n\nSplitting of other sequences will procede but requires some investigation to check if incident is isolated or systemic\n//------------------------//\n";
							}
							$foundhit{$querySeq} = 1;
						}
					}
				}
			}
		}
	}
	
	foreach my $locus (sort keys %loci)
	{
		open OUT1,'>',"$controlArgs{FASTA_OUTPUT}/leftFlank/$locus.fasta";
		foreach my $ind (sort keys %indList)
		{
			foreach my $tax (@{$indArrays{$ind}})
			{
				if (exists $seqs{$tax}{$locus} && exists $coordinates{$tax}{$locus}{start})
				{
						my $chunk = substr($seqs{$tax}{$locus},0,($coordinates{$tax}{$locus}{start} - 1));
						print OUT1 ">$tax\n$chunk\n";
						
						push @{$splitStats{$tax}{$locus}}, length($chunk);
				}
				else
				{
					if (exists $seqs{$tax}{$locus} && ! exists $coordinates{$tax}{$locus}{start})
					{
						print "Problem, $tax has  a sequence but no blast output for $locus, check logs?\nThis is likely caused by the blast hit not clearing the very permissive e-value threshold.\nEither there is not a good reference sequnce, or maybe the data for this locus is dubious and should be removed anyway?\n";
					}
					elsif (! exists $seqs{$tax}{$locus} && exists $coordinates{$tax}{$locus}{start})
					{
						print "Problem, $tax has blast output for $locus but no actual sequence, very strange!\nThere are perhaps some problematic characters in the taxon names like a \"__\", which is used by PATÉ for identifying alleles from a single individual.\n";
					}
				}
			}
		}
		close OUT1;
		open OUT1,'>',"$controlArgs{FASTA_OUTPUT}/core/$locus.fasta";
		foreach my $ind (sort keys %indList)
		{
			foreach my $tax (@{$indArrays{$ind}})
			{
				if (exists $seqs{$tax}{$locus} && exists $coordinates{$tax}{$locus}{start})
				{
						my $dquery = ($coordinates{$tax}{$locus}{end} - $coordinates{$tax}{$locus}{start});
						my $chunk = substr($seqs{$tax}{$locus},($coordinates{$tax}{$locus}{start} - 1),($dquery+1));
						print OUT1 ">$tax\n$chunk\n";
						
						push @{$splitStats{$tax}{$locus}}, length($chunk);
				}
			}
		}
		close OUT1;
		open OUT1,'>',"$controlArgs{FASTA_OUTPUT}/rightFlank/$locus.fasta";
		foreach my $ind (sort keys %indList)
		{
			foreach my $tax (@{$indArrays{$ind}})
			{
				if (exists $seqs{$tax}{$locus} && exists $coordinates{$tax}{$locus}{start})
				{
						my $dquery = (length($seqs{$tax}{$locus}) - $coordinates{$tax}{$locus}{end});
						my $chunk = substr($seqs{$tax}{$locus},$coordinates{$tax}{$locus}{end},$dquery);
						print OUT1 ">$tax\n$chunk\n";
						
						push @{$splitStats{$tax}{$locus}}, length($chunk);
				}
			}
		}
		close OUT1;
	}
	
	open OUT1,'>',"$controlArgs{SUMMARYSTATS_OUT}/averageSplitStats.txt";
	print OUT1 "ID\tnLociSplit\tnFailedSplit\tnMissing\tleftLength\tcoreLength\trightLength\n";
	foreach my $ind (sort keys %indList)
	{
		foreach my $tax (@{$indArrays{$ind}})
		{
			my %avgStats = ();
			open OUT2,'>',"$controlArgs{SUMMARYSTATS_OUT}/$tax.splitstats";
			print OUT2 "Locus\tleftLength\tcoreLength\trightLength\n";
			$avgStats{nls} = 0;
			$avgStats{fs} = 0;
			$avgStats{nm} = 0;
			$avgStats{ll} = 0;
			$avgStats{cl} = 0;
			$avgStats{rl} = 0;
			foreach my $locus (sort keys %loci)
			{
				if (exists $seqs{$tax}{$locus} && exists $coordinates{$tax}{$locus}{start})
				{
					my $ll = ($coordinates{$tax}{$locus}{start} - 1);
					my $cl = ($coordinates{$tax}{$locus}{end} - $coordinates{$tax}{$locus}{start}) + 1;
					my $rl = (length($seqs{$tax}{$locus}) - $coordinates{$tax}{$locus}{end});

					$avgStats{nls} = $avgStats{nls} + 1;
					$avgStats{ll} = $avgStats{ll} + $ll;
					$avgStats{cl} = $avgStats{cl} + $cl;
					$avgStats{rl} = $avgStats{rl} + $rl;
					print OUT2 "$locus\t$ll\t$cl\r$rl\n";
				}
				elsif (exists $seqs{$tax}{$locus} && ! exists $coordinates{$tax}{$locus}{start})
				{
					$avgStats{fs} = $avgStats{fs} + 1;
					print OUT2 "$locus\tFS\tFS\rFS\n";
				}
				elsif (! exists $seqs{$tax}{$locus} && ! exists $coordinates{$tax}{$locus}{start})
				{
					$avgStats{nm} = $avgStats{nm} + 1;
					print OUT2 "$locus\tNA\tNA\rNA\n";
				}
			}
			close OUT2;
			if ($avgStats{nls} > 0)
			{
				$avgStats{ll} = $avgStats{ll}/$avgStats{nls};
				$avgStats{cl} = $avgStats{cl}/$avgStats{nls};
				$avgStats{rl} = $avgStats{rl}/$avgStats{nls};
			}
			print OUT1 "$tax\t$avgStats{nls}\t$avgStats{fs}\t$avgStats{nm}\t$avgStats{ll}\t$avgStats{cl}\t$avgStats{rl}\n";
		}
	}
	close OUT1;
}

########
#Put the flanks and core back together in order and write out RAxML-style and Nexus-formated (e.g. IQTREE) partition files per gene. A single concatenated file and partition file is output as well. I do not agree with this if combined with the phasing output that are intended for some downstream MSC analyses, but it makes the script more general for people that want this outside of PATÉ.
########

elsif ($passedArgs{runmode} == 2)
{
	my %taxa = ();
	my %loci = ();
	my %seqs = ();
	my %seqLengths = ();
	my %partitions = ();
	
	my @leftFlanks = glob("$controlArgs{ALIGNMENT_IN}/leftFlank/*.fa*");
	foreach my $ipfa (@leftFlanks)
	{
		my @temp = ();
		@temp = split(/\//,$ipfa);
		if ($temp[scalar(@temp)-1] =~ m/(\S+)\.fa.*/)
		{
			my $locus = $1;
			if ($controlArgs{DELIMITER} ne "0")
			{
				my @prefixSplitter = ();
				@prefixSplitter = split(/$controlArgs{DELIMITER}/,$locus);
				$locus = $prefixSplitter[0];
#				print "$locus\n"
			}
			if (! exists $loci{$locus})
			{
				$loci{$locus} = 1;
			}
			
			my $tax = "";
			my $seq = "";
			open FH1,'<',"$ipfa";
			while (<FH1>)
			{
				if (/^>(\S+)/)
				{
					$tax = $1;
					if (! exists $taxa{$tax})
					{
						$taxa{$tax} = 1;
					}
					if (! exists $seqs{$tax}{$locus})
					{
						$seqs{$tax}{$locus} = "";
						$seqLengths{$tax}{$locus} = 0;
						$partitions{$tax}{$locus}[0] = $seqLengths{$tax}{$locus};
					}
					elsif (exists $seqs{$tax}{$locus})
					{
						print "//--------WARNING!--------//\nLocus $locus has repeating taxon label $tax\n\nThis will cause problems here and downstream.\n//------------------------//\n";
					}
				}
				elsif (/(\S+)/)
				{
					$seq = $1;
					$seqs{$tax}{$locus} = $seqs{$tax}{$locus} . $seq;
					$seqLengths{$tax}{$locus} = $seqLengths{$tax}{$locus} + length($seq);
					$partitions{$tax}{$locus}[1] = $seqLengths{$tax}{$locus};
#					print "$tax\t$locus\t$partitions{$tax}{$locus}[0]\t$partitions{$tax}{$locus}[1]\n";
				}
			}
			close FH1;
		}
	}
	my @cores = glob("$controlArgs{ALIGNMENT_IN}/core/*.fa*");
	foreach my $ipfa (@cores)
	{
		my @temp = ();
		@temp = split(/\//,$ipfa);
		if ($temp[scalar(@temp)-1] =~ m/(\S+)\.fa.*/)
		{
			my $locus = $1;
			if ($controlArgs{DELIMITER} ne "0")
			{
				my @prefixSplitter = ();
				@prefixSplitter = split(/$controlArgs{DELIMITER}/,$locus);
				$locus = $prefixSplitter[0];
			}
			if (! exists $loci{$locus})
			{
				$loci{$locus} = 1;
				print "//--------WARNING!--------//\nLocus $locus has an aligned core but no left flanking sequence\n\nCheck data if flanking regions should have been assembled, split, and aligned\n//------------------------//\n";
			}
			
			my $tax = "";
			my $seq = "";
			open FH1,'<',"$ipfa";
			while (<FH1>)
			{
				if (/^>(\S+)/)
				{
					$tax = $1;
					if (! exists $taxa{$tax})
					{
						$taxa{$tax} = 1;
					}
					if (exists $seqs{$tax}{$locus})
					{
						$partitions{$tax}{$locus}[2] = $partitions{$tax}{$locus}[1] + 1;
					}
					elsif (! exists $seqs{$tax}{$locus})
					{
						print "//--------WARNING!--------//\n$tax at locus $locus has an aligned core but no left flanking sequence\n\nCheck data if flanking regions should have been assembled, split, and aligned\nThe left flank will be replaced by Ns in the concatenated alignment//------------------------//\n";
						$seqs{$tax}{$locus} = "";
						my $probablyCorrectPartitionStart = 0;
						my $probablyCorrectPartitionEnd = 0;
						foreach my $ind (sort keys %taxa)
						{
							if (exists $partitions{$ind}{$locus}[0])
							{
								$probablyCorrectPartitionStart = $partitions{$ind}{$locus}[0];
								$probablyCorrectPartitionEnd = $partitions{$ind}{$locus}[1];
							}
						}
						for my $i ($probablyCorrectPartitionStart..$probablyCorrectPartitionEnd)
						{
							$seqs{$tax}{$locus} = $seqs{$tax}{$locus} . "N";
						}
						$partitions{$tax}{$locus}[0] = $probablyCorrectPartitionStart;
						$partitions{$tax}{$locus}[1] = $probablyCorrectPartitionEnd;
						$partitions{$tax}{$locus}[2] = $partitions{$tax}{$locus}[1] + 1;
						$seqLengths{$tax}{$locus} = $partitions{$tax}{$locus}[1];
					}
				}
				elsif (/(\S+)/)
				{
					$seq = $1;
					$seqs{$tax}{$locus} = $seqs{$tax}{$locus} . $seq;
					$seqLengths{$tax}{$locus} = $seqLengths{$tax}{$locus} + length($seq);
					$partitions{$tax}{$locus}[3] = $seqLengths{$tax}{$locus};
				}
			}
			close FH1;
		}
	}
	my @rightFlanks = glob("$controlArgs{ALIGNMENT_IN}/rightFlank/*.fa*");
	foreach my $ipfa (@rightFlanks)
	{
		my @temp = ();
		@temp = split(/\//,$ipfa);
		if ($temp[scalar(@temp)-1] =~ m/(\S+)\.fa.*/)
		{
			my $locus = $1;
			if ($controlArgs{DELIMITER} ne "0")
			{
				my @prefixSplitter = ();
				@prefixSplitter = split(/$controlArgs{DELIMITER}/,$locus);
				$locus = $prefixSplitter[0];
			}
			if (! exists $loci{$locus})
			{
				$loci{$locus} = 1;
				print "//--------WARNING!--------//\nLocus $locus has an aligned right flank but no core or left flank\n\nThis should not happen and there was a problem, potentially the references are not appropriate or the wrong blast program was used\n//------------------------//\n";
			}
			
			my $tax = "";
			my $seq = "";
			open FH1,'<',"$ipfa";
			while (<FH1>)
			{
				if (/^>(\S+)/)
				{
					$tax = $1;
					if (! exists $taxa{$tax})
					{
						$taxa{$tax} = 1;
					}
					if (exists $seqs{$tax}{$locus})
					{
						$partitions{$tax}{$locus}[4] = $partitions{$tax}{$locus}[3] + 1;
					}
					elsif (! exists $seqs{$tax}{$locus})
					{
						print "//--------WARNING!--------//\n$tax at locus $locus has an aligned right flank but no left flanking or core sequence\n\nCheck data if flanking regions should have been assembled, split, and aligned\nThe left flank will be replaced by Ns in the concatenated alignment//------------------------//\n";
						$seqs{$tax}{$locus} = "";
						$seqLengths{$tax}{$locus} = 0;
						my $probablyCorrectPartitionStart1 = 0;
						my $probablyCorrectPartitionEnd1 = 0;
						my $probablyCorrectPartitionStart2 = 0;
						my $probablyCorrectPartitionEnd2 = 0;
						foreach my $ind (sort keys %taxa)
						{
							if (exists $partitions{$ind}{$locus}[0])
							{
								$probablyCorrectPartitionStart1 = $partitions{$ind}{$locus}[0];
								$probablyCorrectPartitionEnd1 = $partitions{$ind}{$locus}[1];
								$probablyCorrectPartitionStart2 = $partitions{$ind}{$locus}[2];
								$probablyCorrectPartitionEnd2 = $partitions{$ind}{$locus}[3];
							}
						}
						for my $i ($probablyCorrectPartitionStart1..$probablyCorrectPartitionEnd2)
						{
							$seqs{$tax}{$locus} = $seqs{$tax}{$locus} . "N";
						}
						$partitions{$tax}{$locus}[0] = $probablyCorrectPartitionStart1;
						$partitions{$tax}{$locus}[1] = $probablyCorrectPartitionEnd1;
						$partitions{$tax}{$locus}[2] = $probablyCorrectPartitionStart2;
						$partitions{$tax}{$locus}[3] = $probablyCorrectPartitionEnd2;
						$partitions{$tax}{$locus}[4] = $partitions{$tax}{$locus}[3] + 1;
						$seqLengths{$tax}{$locus} = $partitions{$tax}{$locus}[3];
					}
				}
				elsif (/(\S+)/)
				{
					$seq = $1;
					$seqs{$tax}{$locus} = $seqs{$tax}{$locus} . $seq;
					$seqLengths{$tax}{$locus} = $seqLengths{$tax}{$locus} + length($seq);
					$partitions{$tax}{$locus}[5] = $seqLengths{$tax}{$locus};
#					print "$tax\t$locus\t$partitions{$tax}{$locus}[0]\t$partitions{$tax}{$locus}[1]\t$partitions{$tax}{$locus}[2]\t$partitions{$tax}{$locus}[3]\t$partitions{$tax}{$locus}[4]\t$partitions{$tax}{$locus}[5]\n";
				}
			}
			close FH1;
		}
	}

########
#Write out new fasta files and partition files
########
	my %concatenatedSeqs = ();
	my %concatenatedFlanks = ();
	my %concatenatedPartitions = ();
	my $totalLength = 0;
	
	foreach my $locus (sort keys %loci)
	{
		open OUT1,'>',"$controlArgs{ALIGNMENT_OUT}/loci/$locus.fasta";
		open OUT2,'>',"$controlArgs{ALIGNMENT_OUT}/loci/$locus.partitions.raxml";
		open OUT3,'>',"$controlArgs{ALIGNMENT_OUT}/loci/$locus.partitions.nexus";
		open OUT4,'>',"$controlArgs{ALIGNMENT_OUT}/loci/$locus.flanksOnly.fasta";
		
		my @partitionChecker = (-1,-1,-1,-1,-1,-1);
		my $partitionsPrinted = 0;
		my $locusLength = 0;
		my $flankOnlyLength = 0;
		
		foreach my $tax (sort keys %taxa)
		{
			if (exists $seqs{$tax}{$locus})
			{
				print OUT1 ">$tax\n$seqs{$tax}{$locus}\n";
				my $flankSeq1 = substr($seqs{$tax}{$locus},$partitions{$tax}{$locus}[0],($partitions{$tax}{$locus}[1] - $partitions{$tax}{$locus}[0] + 1));
				my $flankSeq2 = substr($seqs{$tax}{$locus},$partitions{$tax}{$locus}[4],($partitions{$tax}{$locus}[5] - $partitions{$tax}{$locus}[4] + 1));
				my $superFlank = $flankSeq1 . $flankSeq2;
				print OUT4 ">$tax\n$superFlank\n";
				if (length($seqs{$tax}{$locus}) > $locusLength)
				{
					$locusLength = length($seqs{$tax}{$locus});
					$flankOnlyLength = length($superFlank);
				}		
				
				if (! exists $concatenatedSeqs{$tax})
				{
					$concatenatedSeqs{$tax} = $seqs{$tax}{$locus};
					
				}
				elsif (exists $concatenatedSeqs{$tax})
				{
					$concatenatedSeqs{$tax} = $concatenatedSeqs{$tax} . $seqs{$tax}{$locus};
					$concatenatedFlanks{$tax} = $concatenatedFlanks{$tax} . $superFlank;
				}
				
				if (! exists $concatenatedFlanks{$tax})
				{
					$concatenatedFlanks{$tax} = $superFlank;
				}
				elsif (exists $concatenatedSeqs{$tax})
				{
					$concatenatedFlanks{$tax} = $concatenatedFlanks{$tax} . $superFlank;
				}
				if ($partitionChecker[0] == -1)
				{
#					print "New $tax for $locus with partitions:\n";
					for my $i (0..5)
					{
						$partitionChecker[$i] = $partitions{$tax}{$locus}[$i];
#						print "\t$i = $partitions{$tax}{$locus}[$i]\n";
					}
				}
				elsif ($partitionChecker[0] != -1)
				{
					for my $i (0..5)
					{
						if ($partitionChecker[$i] != $partitions{$tax}{$locus}[$i])
						{
							print "//--------SERIOUS PROBLEM!--------//\n\nPartition coordinates between sequences starting at individual $tax for locus $locus do not match others in alignment.\n\nThis implies that sequences are not aligned and requires inspection!\n//--------------------------------//\n";
						}
					}
				}
			}
		}
		foreach my $tax (sort keys %taxa)
		{
			if (! exists $seqs{$tax}{$locus})
			{
				if (! exists $concatenatedSeqs{$tax})
				{
					$concatenatedSeqs{$tax} = "";
					$concatenatedFlanks{$tax} = "";
				}
				for my $i (0..($locusLength - 1))
				{
					$concatenatedSeqs{$tax} = $concatenatedSeqs{$tax} . "N";
				}
			}
		}
		
		if ($partitionsPrinted == 0)
		{
			#do the stuff
			my $shiftedp1 = $partitionChecker[0] + 1;
			my $shiftedp2 = $partitionChecker[1];
			my $shiftedp3 = $partitionChecker[2];
			my $shiftedp4 = $partitionChecker[3];
			my $shiftedp5 = $partitionChecker[4];
			my $shiftedp6 = $partitionChecker[5];
			
#			print "$locus\t$shiftedp1\t$shiftedp2\t$shiftedp3\t$shiftedp4\t$shiftedp5\t$shiftedp6\n";
					
			print OUT2 "DNA, $locus.leftFlank = $shiftedp1-$shiftedp2\nDNA, $locus.core = $shiftedp3-$shiftedp4\nDNA, $locus.rightFlank = $shiftedp5-$shiftedp6\n";
				
			print OUT3 "#nexus\nbegin sets;\ncharset $locus.leftFlank = $shiftedp1-$shiftedp2;\ncharset $locus.core = $shiftedp3-$shiftedp4;\ncharset $locus.rightFlank = $shiftedp5-$shiftedp6;\nend;";
					
			
					
			$totalLength = $totalLength + $locusLength;
					
			$concatenatedPartitions{$locus}[0] = ($totalLength - $locusLength + $shiftedp1);
			$concatenatedPartitions{$locus}[1] = ($totalLength - $locusLength + $shiftedp2);
			$concatenatedPartitions{$locus}[2] = ($totalLength - $locusLength + $shiftedp3);
			$concatenatedPartitions{$locus}[3] = ($totalLength - $locusLength + $shiftedp4);
			$concatenatedPartitions{$locus}[4] = ($totalLength - $locusLength + $shiftedp5);
			$concatenatedPartitions{$locus}[5] = ($totalLength - $locusLength + $shiftedp6);
			
#			print "TotalLength = $totalLength\n";
#			print "LocusLength = $locusLength\n";
#			print "$concatenatedPartitions{$locus}[0]\t$concatenatedPartitions{$locus}[1]\t$concatenatedPartitions{$locus}[2]\t$concatenatedPartitions{$locus}[3]\t$concatenatedPartitions{$locus}[4]\t$concatenatedPartitions{$locus}[5]\n";
			
			$partitionsPrinted = 1;
		}
		close OUT1;
		close OUT2;
		close OUT3;
	}
	
	open OUT1,'>',"$controlArgs{ALIGNMENT_OUT}/concatenated/concatenated.fasta";
	open OUT2,'>',"$controlArgs{ALIGNMENT_OUT}/concatenated/concatenated.flanksOnly.fasta";
	foreach my $tax (sort keys %taxa)
	{
		print OUT1 ">$tax\n$concatenatedSeqs{$tax}\n";
		print OUT2 ">$tax\n$concatenatedFlanks{$tax}\n";
	}
	close OUT1;
	close OUT2;
	
	open OUT1,'>',"$controlArgs{ALIGNMENT_OUT}/concatenated/concatenated.partitions.raxml";
	open OUT2,'>',"$controlArgs{ALIGNMENT_OUT}/concatenated/concatenated.partitions.nexus";
	print OUT2 "#nexus\nbegin sets;\n";
	foreach my $locus (sort keys %loci)
	{
		print OUT1 "DNA, $locus.leftFlank = $concatenatedPartitions{$locus}[0]-$concatenatedPartitions{$locus}[1]\nDNA, $locus.core = $concatenatedPartitions{$locus}[2]-$concatenatedPartitions{$locus}[3]\nDNA, $locus.rightFlank = $concatenatedPartitions{$locus}[4]-$concatenatedPartitions{$locus}[5]\n";			
		print OUT2 "charset $locus.leftFlank = $concatenatedPartitions{$locus}[0]-$concatenatedPartitions{$locus}[1];\ncharset $locus.core = $concatenatedPartitions{$locus}[2]-$concatenatedPartitions{$locus}[3];\ncharset $locus.rightFlank = $concatenatedPartitions{$locus}[4]-$concatenatedPartitions{$locus}[5];\n";
	}
	print OUT2 "end;\n";
	close OUT2;
	close OUT1;
}

##########################################################################################
#Subroutine for complementing strings of nucleotides - does not handle ambiguity codes
##########################################################################################
sub Complement
{
	my @CharArray = split//,$_[0];
	my $CompGene = "";
	foreach my $Char (@CharArray)
	{
		my $CompChar = "";
		if ($Char eq "A")
		{
			$CompChar = "T";
		}
		if ($Char eq "T")
		{
			$CompChar = "A";
		}
		if ($Char eq "C")
		{
			$CompChar = "G";
		}
		if ($Char eq "G")
		{
			$CompChar = "C";
		}
		if ($Char ne "A" and $Char ne "T" and $Char ne "C" and $Char ne "G" )
		{
			$CompChar = "N";
		}
		$CompGene = $CompGene . $CompChar;
	}
	$CompGene = reverse($CompGene);
	return $CompGene;
}
exit;