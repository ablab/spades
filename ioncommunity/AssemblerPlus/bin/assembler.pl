#!/usr/bin/perl -w
# Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved
use strict;
 
#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 2
#of the License, or (at your option) any later version.
 
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
 
#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.



#force scope
{
    my ($bamFile, $reference, $fractionOfReads, $assemblerPath, $outputDirectory, $urlRoot, $sampleName, $miraversion, $runMira, $runSpades, $spadesversion, $spadesOptions, $quastOnly, $qref) = @ARGV;

    $fractionOfReads=1 if $fractionOfReads eq ''||$fractionOfReads==0||$fractionOfReads>1;
    $runMira=0 if $runMira eq ''; 
    $runSpades=0 if $runSpades eq '';

    $miraversion="4.0" if $miraversion eq ''; # should be already set if called by launch.sh
    $spadesversion="3.0.0" if $spadesversion eq '';

    $reference=$ENV{'TSP_FILEPATH_GENOME_FASTA'} if $reference eq '';
    $qref=$ENV{'TSP_FILEPATH_GENOME_FASTA'} if $qref eq '';
    my $type = defined $ENV{"PLUGINCONFIG__TYPE"} ? $ENV{"PLUGINCONFIG__TYPE"} : "accurate";

    my $sffExtractPath = $assemblerPath."/sff_extract";
    my $downsampleBAMPath = $assemblerPath."/DownsampleSam.jar";

    my $spadesPath = $assemblerPath."/SPAdes-$spadesversion-Linux/bin/spades.py";

    #########################################################################################################
    ############################## edit the line below when upgrading QUAST #################################
    my $quastPath = $assemblerPath."/quast-2.3/quast.py"; ###################################################
    #########################################################################################################

    my $miraPath = $assemblerPath."/mira-${miraversion}/bin/mira";

    #print out the input files received
    print "assembler.pl Input:\n";
    print "\tBAM File: $sampleName/$bamFile\n";
    print "\tSff Extract Path: $sffExtractPath\n";
    print "\tDownsample BAM Path: $downsampleBAMPath\n";
    print "\tFraction of Reads: $fractionOfReads\n";
    print "\tOutput Directory: $outputDirectory\n";
    print "\tURL Root: $urlRoot\n";
    print "\tMira Path: $miraPath\n";
    print "\tMIRA Reference: $reference\n";
    print "\tRun MIRA: $runMira\n";
    print "\tMIRA version: $miraversion\n";
    print "\tRun SPAdes: $runSpades\n";
    print "\tSPAdes options: $spadesOptions\n";
    print "\tSPAdes version: $spadesversion\n";
    print "\tQuast Only: $quastOnly\n";
    print "\tQuast reference: $qref\n\n";

    #see how many reads there are
    #what if not barcoded?
    my $qualityFile = $bamFile;
    $qualityFile =~ s/basecaller.bam/quality.summary/;
    print "\nExtracting read count\n";
    my $systemCall = "";
    my $numReads = "";

    #get the total number of reads          
    $systemCall = "samtools view $sampleName/$bamFile -c";
    $numReads = `$systemCall`;
    chomp($numReads);
    print "$numReads reads in $bamFile\n";

    if(defined $ENV{"PLUGINCONFIG__MIN_READS"} && $numReads<=$ENV{"PLUGINCONFIG__MIN_READS"}){
	print "\t*Does not have more than ", $ENV{"PLUGINCONFIG__MIN_READS"}, " reads. Skipping this file\n";
	exit(0);
    }

    #if fraction is not 1 then we will subsample. 
    my $bamRandomCall1 = "";
    my $bamRandomCall2 = "";
    if($fractionOfReads < 1){
	print "\nSubsampling using picard\n";
	my $newBamFile = "$outputDirectory/$sampleName/$bamFile\_scaled";
	$bamRandomCall1 = "java -Xmx2g -jar $downsampleBAMPath INPUT=$outputDirectory/$sampleName/$bamFile OUTPUT=$newBamFile PROBABILITY=$fractionOfReads";
	print "$bamRandomCall1\n";
	system($bamRandomCall1);
	$bamRandomCall2 = "mv $outputDirectory/$sampleName/$bamFile\_scaled $outputDirectory/$sampleName/$bamFile";
	print "$bamRandomCall2\n";
	system($bamRandomCall2);
    }

    #convert BAM to SFF
    print "\nRunning bam2sff\n";
    my $sffFile = $bamFile;
    $sffFile =~ s/bam/sff/;
    my $bam2sffCall = "bam2sff -o $outputDirectory/$sampleName/$sffFile $outputDirectory/$sampleName/$bamFile";
    print "$bam2sffCall\n";
    system($bam2sffCall);

    #run sff extract to generate fastq and trace xml
    my $fastq = $outputDirectory."/".$sampleName."/".$sffFile;
    my $xml = $outputDirectory."/".$sampleName."/".$sffFile;
    $fastq =~ s/.sff$/\_in.iontor.fastq/;
    $xml =~ s/.sff$/\_traceinfo_in.iontor.xml/;

    #my $command1 = "$sffExtractPath --min_left_clip=5 -s $fastq -x $xml $outputDirectory/$sampleName/$sffFile";
    my $command1 = "$sffExtractPath -s $fastq -x $xml $outputDirectory/$sampleName/$sffFile";
    print "\nRunning sff_extract\n$command1\n\n";
    system($command1);

    my $projectName = $sffFile;
    $projectName =~ s/\.sff$//;
    my $command2 = "";
    my $command4 = "";
    my $sampleID = $sampleName; $sampleID=~s/\.[ATCG]+$//i;

    if ($runSpades) {
      if  ( ! ($quastOnly && -s "$sampleName/spades/contigs.fasta")) {
        $spadesOptions .= " --tmp-dir /tmp/$$ -s $fastq -o $sampleName/spades";

        if ($spadesversion eq "2.5.1") {
                $spadesOptions .= " --rectangles";
        } else { # >= 3.0.0
                $spadesOptions .= " --iontorrent";
        } 
        $command2 = "$spadesPath $spadesOptions > /dev/null";
        print "\nRunning AssemblerPlus - SPAdes\n$command2\n";
        system($command2);
      }
      if($qref ne "None") {
        $command4 = "$quastPath -R $qref $sampleName/spades/contigs.fasta; rm -fr $sampleName/spades/quast_results; mv quast_results $sampleName/spades";
        print "\nRunning QUAST\n$command4\n";
        system($command4);
        my $report=`ls $sampleName/spades/quast_results/*/report.html`; chomp $report;
        if ($report)
          {`echo "<p align=left><a target=_blank href=$report>$sampleID SPAdes-QUAST report</a> <a target=_blank href=$sampleName/spades/contigs.fasta>contigs.fasta</a>" >> SPAdes_QUAST.html`
          }
        else {`echo "<p align=left><a target=_blank href=$sampleName/spades/spades.log>$sampleID failed</a>" >> SPAdes_QUAST.html`}
      }
    }

    if ($runMira) {

      if ( ! ($quastOnly && -s "$sampleName/$projectName\_assembly/$projectName\_d_results\/$projectName\_out.unpadded.fasta")) {

        #build mira executable command and run it
        #if mira 3.4.1.1
        if($miraversion eq "3.4.1.1"){
	  $command2 = "cd ./$sampleName; $miraPath --project=$projectName --job=denovo,genome,accurate,iontor -MI:sonfs=no IONTOR_SETTINGS -AS:mrpc=100 > $outputDirectory/$sampleName/mira.log";

	  #if they provided a genome we will use that here
	  if($reference ne "None"){
	     $command2 = "cd ./$sampleName; $miraPath --project=$projectName --job=denovo,genome,$type,iontor -MI:sonfs=no -SB:abnc=yes:bft=fasta -FN:bbin=$reference IONTOR_SETTINGS -AS:mrpc=100 > $outputDirectory/$sampleName/mira.log";
	  }

	  print "\nRunning AssemblerPlus - Mira\n$command2\n";
	  system($command2);
        }
        #otherwise 3.9.9 / 4.0.x
        else{
	  #need to create the  manfest file
	  my $manifestContent = "\n\#MIRA Manifest File\n\n";
	  $manifestContent .= "#Settings\n";
	  $manifestContent .= "project = $projectName\n";
	  $manifestContent .= "job = denovo,genome,$type\n";

	  #see if reference assisted or not
	  if($reference eq "None"){
            $manifestContent .= "parameters = -MI:sonfs=no IONTOR_SETTINGS -AS:mrpc=100 \n";
  	  }
	  else{
	    $manifestContent .= "parameters = -MI:sonfs=no IONTOR_SETTINGS -AS:mrpc=100 \n";
	    if($miraversion eq "3.9.9") {print "Reference assisted assemblies not currently supported in MIRA 3.9.9 development version\n";}
	    #$manifestContent .= "parameters = -MI:sonfs=no -SB:abnc=yes IONTOR_SETTINGS -AS:mrpc=100 TEXT_SETTINGS -LR:wqf=no\n\n";
	  }

          if($miraversion =~ /4\.0\./) {$manifestContent=~s/-MI:sonfs=no /-DI:trt=\/tmp -MI:/}

  	  $manifestContent .= "#Reads\n";
	  $manifestContent .= "readgroup = $sampleName\n";
	  $manifestContent .= "data = $fastq $xml\n";
	  $manifestContent .= "technology = iontor\n\n";

	  #add reference read group
	  if($reference ne "None"){
	    #$manifestContent .= "#Reference\n";
	    #$manifestContent .= "readgroup = Reference\n";
	    #$manifestContent .= "data = $reference\n";
	    #$manifestContent .= "technology = text\n";
	    #$manifestContent .= "as_reference\n";
	  }

	  #write the manfiest file
	  open(OUT, ">$outputDirectory/$sampleName/manfiest.txt") || die "Could not write manifest file for mira at $outputDirectory/$sampleName/manifest.txt\n";
	  print OUT "$manifestContent";
	  close(OUT);

	  #execute MIRA
	  $command2 = "cd ./$sampleName; $miraPath $outputDirectory/$sampleName/manfiest.txt > $outputDirectory/$sampleName/mira.log";
	  print "\nRunning AssemblerPlus - Mira\n$command2\n";
	  system($command2);

	  #append manifest content to command so it gets printed in logs
	  $command2 .= $manifestContent;
        }

      }

      #parse MIRA with QUAST
      if($qref ne "None") {
        $command4 = "$quastPath -R $qref $sampleName/$projectName\_assembly/$projectName\_d_results\/$projectName\_out.unpadded.fasta; rm -fr $sampleName/quast_results; mv quast_results $sampleName";
        print "\nRunning QUAST\n$command4\n";
        system($command4);
        my $report=`ls $sampleName/quast_results/*/report.html`; chomp $report;
        if ($report)
          {`echo "<p align=left><a target=_blank href=$report>$sampleID MIRA-QUAST report</a>&nbsp; &nbsp;<a target=_blank href=$sampleName/$projectName\_assembly/$projectName\_d_results\/$projectName\_out.unpadded.fasta>contigs.fasta</a>" >> MIRA_QUAST.html`
          }
        else {`echo "<p align=left><a target=_blank href=$sampleName/mira.log>$sampleID failed</a>" >> MIRA_QUAST.html`}
      }

#alternative parsing of MIRA output
   if (1) {

    #assembled reads total coverage, num contigs, N50, N90, N95 
    my $assembledReads = `grep "Num. reads assembled" $outputDirectory/$sampleName/$projectName\_assembly/$projectName\_d_info/$projectName\_info_assembly.txt | cut -f 2 -d ":" | sed -e 's/ //'`;
    my $coverage = `grep "Avg. total coverage" $outputDirectory/$sampleName/$projectName\_assembly/$projectName\_d_info/$projectName\_info_assembly.txt | cut -f 2 -d ":" | sed -e 's/ //'`;

    #these will all return two values: all contigs and using just the large contigs
    my $numContig = `grep "Number of contigs" $outputDirectory/$sampleName/$projectName\_assembly/$projectName\_d_info/$projectName\_info_assembly.txt | cut -f 2 -d ":" | sed -e 's/[ \t]//'`;
    my $consensusLength = `grep "Total consensus" $outputDirectory/$sampleName/$projectName\_assembly/$projectName\_d_info/$projectName\_info_assembly.txt | cut -f 2 -d ":" | sed -e 's/[ \t]//'`;
    my $largestContig = `grep "Largest contig" $outputDirectory/$sampleName/$projectName\_assembly/$projectName\_d_info/$projectName\_info_assembly.txt | cut -f 2 -d ":" | sed -e 's/[ \t]//'`;
    my $N50 = `grep "N50 contig size" $outputDirectory/$sampleName/$projectName\_assembly/$projectName\_d_info/$projectName\_info_assembly.txt | cut -f 2 -d ":" | sed -e 's/[ \t]//'`;
    my $N90 = `grep "N90 contig size" $outputDirectory/$sampleName/$projectName\_assembly/$projectName\_d_info/$projectName\_info_assembly.txt | cut -f 2 -d ":" | sed -e 's/[ \t]//'`;
    my $N95 = `grep "N95 contig size" $outputDirectory/$sampleName/$projectName\_assembly/$projectName\_d_info/$projectName\_info_assembly.txt | cut -f 2 -d ":" | sed -e 's/[ \t]//'`;
    
    chomp($assembledReads);
    chomp($coverage);
    chomp($numContig);
    chomp($consensusLength);
    chomp($largestContig);
    chomp($N50);
    chomp($N90);
    chomp($N95);

    my @numContigs = split(/\n/, $numContig);
    my @consensusLengths = split(/\n/, $consensusLength);
    my @largestContigs = split(/\n/, $largestContig);
    my @N50s = split(/\n/, $N50);
    my @N90s = split(/\n/, $N90);
    my @N95s = split(/\n/, $N95);

    #do some number formatting
    $assembledReads =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
    $coverage =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
    $numContigs[0] =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
    $numContigs[1] =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
    $consensusLengths[0] =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
    $consensusLengths[1] =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
    $largestContigs[0] =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
    $N50s[0] =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
    $N50s[1] =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
    $N90s[0] =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
    $N90s[1] =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
    $N95s[0] =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
    $N95s[1] =~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;

    #need to get total reads, API?
    print "\tAssembled Reads: $assembledReads\n";
    print "\tCoverage: $coverage\n";
    print "\tConsensus Lengths: $consensusLengths[0]\t$consensusLengths[1]\n";
    print "\tContigs: $numContigs[0]\t$numContigs[1]\n";
    print "\tLargest Contig: $largestContigs[0]\n";
    print "\tN50: $N50s[0]\t$N50s[1]\n";
    print "\tN90: $N90s[0]\t$N90s[1]\n";
    print "\tN95: $N95s[0]\t$N95s[1]\n";

    #print commands to txt file
    open(OUT, ">$outputDirectory/$sampleName/commands.txt") || die "Couldn't write command text file\n";

    print OUT "#OPTIONS\n";
    print OUT "SFF File: $sampleName/$sffFile\n";
    print OUT "Mira Path: $miraPath\n";
    print OUT "Sff Extract Path: $sffExtractPath\n";
    print OUT "Reference: $reference\n";
    print OUT "Fraction of Reads: $fractionOfReads\n";
    print OUT "Output Directory: $outputDirectory\n";
    print OUT "URL Root: $urlRoot\n\n";

    print OUT "#COMMANDS\n";

    if($bamRandomCall1 ne ""){
	print OUT "$bamRandomCall1\n";
	print OUT "$bamRandomCall2\n";
    }

    print OUT "$bam2sffCall\n";
    print OUT "$command1\n";
    print OUT "$command2\n";

    close(OUT);

    #print the data to a temp file that the report formater will pick up
    open(OUT, ">>$outputDirectory/assembly_summaries.txt") || die "Could not write assembly_summaries.txt file\n";    
    print OUT "$sampleName\t$assembledReads\t$coverage\t$consensusLengths[0]\t$consensusLengths[1]\t$numContigs[0]\t$numContigs[1]\t$largestContigs[0]\t$N50s[0]\t$N50s[1]\t$N90s[0]\t$N90s[1]\t$N95s[0]\t$N95s[1]\t$sampleName/commands.txt\t$sampleName\/$projectName\_assembly\/$projectName\_d_results\/$projectName\_out.unpadded.fasta\t$sampleName\/$projectName\_assembly\/$projectName\_d_results\/$projectName\_out.wig\t$sampleName\/$projectName\_assembly\/$projectName\_d_results\/$projectName\_out.ace\t$sampleName\/$projectName\_assembly\/$projectName\_d_info\/$projectName\_info_assembly.txt\t$sampleName\/mira.log\n";
    close(OUT);

   }

  }

}
