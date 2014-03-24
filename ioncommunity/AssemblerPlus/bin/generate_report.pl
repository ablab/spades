#!/usr/bin/perl -w
# Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved
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
    my ($outputDirectory) = @ARGV;
    
    #create th HTML doc and build the framework
    open(HTML, ">$outputDirectory/MIRA_classic.html") || die "Could not write MIRA_classic.html\n";

    print HTML "<!DOCTYPE html>\n";
    print HTML "<html>\n";
    print HTML "\t<head>\n";
    print HTML "\t\t<title>AssemblerPlus Plugin</title>\n";
    print HTML "\t\t<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n";
    print HTML "\t\t<link href=\"/pluginMedia/AssemblerPlus/css/kendo.common.min.css\" rel=\"stylesheet\" />\n";
    print HTML "\t\t<link href=\"/pluginMedia/AssemblerPlus/css/kendo.default.min.css\" rel=\"stylesheet\" />\n";
    print HTML "\t\t<link href=\"/pluginMedia/AssemblerPlus/css/kendo.ir.css\" rel=\"stylesheet\" />\n";
    print HTML "\t\t<link href=\"/pluginMedia/AssemblerPlus/css/ir.css\" rel=\"stylesheet\">\n";
    print HTML "\t\t<link href=\"/pluginMedia/AssemblerPlus/css/app.css\" rel=\"stylesheet\">\n";
    print HTML "\t\t<link href=\"/pluginMedia/AssemblerPlus/css/bootstrap.css\" rel=\"stylesheet\">\n";
    print HTML "\t\t<link href=\"/pluginMedia/AssemblerPlus/css/app.less\" rel=\"stylesheet/less\">\n";
    print HTML "\t\t<link href=\"/pluginMedia/AssemblerPlus/css/bootstrap-custom.css\" rel=\"stylesheet\">\n";
    print HTML "\t\t<link href=\"/pluginMedia/AssemblerPlus/css/bootstrap-select.min.css\" rel=\"stylesheet\">\n";
    print HTML "\t<head>\n";
    print HTML "\t<body>\n";
    print HTML "\t\t<script src=\"/pluginMedia/AssemblerPlus/js/less-1.4.1.min.js\"></script>\n";
    print HTML "\t\t<script src=\"/pluginMedia/AssemblerPlus/js/jquery-1.8.2.min.js\"></script>\n";
    print HTML "\t\t<script src=\"/pluginMedia/AssemblerPlus/js/bootstrap-select.min.js\"></script>\n";
    print HTML "\t\t<script src=\"/pluginMedia/AssemblerPlus/js/bootstrap.min.js\"></script>\n";
    print HTML "\t\t<script>\n";
    print HTML "\t\t\tfunction showSample(sample){\n";
    print HTML "\t\t\t\tvar divCollection = document.getElementsByTagName(\"div\");\n";
    print HTML "\t\t\t\tfor (var i=0; i<divCollection.length; i++) {\n";
    print HTML "\t\t\t\t\tif(divCollection[i].getAttribute(\"class\") == \"row-fluid sample\") {\n";
    print HTML "\t\t\t\t\t\tif(divCollection[i].getAttribute(\"id\") == sample) {\n";
    print HTML "\t\t\t\t\t\t\tdivCollection[i].style.display = 'block';\n";
    print HTML "\t\t\t\t\t\t}\n";
    print HTML "\t\t\t\t\t\telse{\n";
    print HTML "\t\t\t\t\t\t\tdivCollection[i].style.display = 'none';\n";
    print HTML "\t\t\t\t\t\t}\n";
    print HTML "\t\t\t\t\t}\n";
    print HTML "\t\t\t\t}\n";
    print HTML "\t\t\t}\n";
    print HTML "\t\t</script>\n";
    print HTML "\n";

    print HTML "\t\t<div class=\"main\">\n";
    print HTML "\t\t\t<div class=\"main-content clearfix\">\n";
    print HTML "\t\t\t\t<div class=\"container-fluid\">\n";
    print HTML "\t\t\t\t\t<div id=\"alertUser\" class=\"row-fluid\">\n";
    print HTML "\t\t\t\t\t\t<div class=\"span12\">\n";
    print HTML "\t\t\t\t\t\t\t<div class=\"alert alert-info\">\n";
    print HTML "\t\t\t\t\t\t\t\t<button class=\"close\" data-dismiss=\"alert\" type=\"button\">x</button>\n";
    print HTML "\t\t\t\t\t\t\t\tAssemblies were performed using <a href=\"http://mira-assembler.sourceforge.net\" target=\"_blank\">MIRA</a>. If you have questions on the quality of the assembly please refer to the <a href=\"http://mira-assembler.sourceforge.net\" target=\"_blank\">MIRA</a> site.\n";
    print HTML "\t\t\t\t\t\t\t</div>\n";
    print HTML "\t\t\t\t\t\t</div>\n";
    print HTML "\t\t\t\t\t</div>\n";
    print HTML "\n";

    #read in the summary file
    my $samples = {};
    open(IN, "$outputDirectory/assembly_summaries.txt") || die "Could not read assembly_summaries.txt file\n";    
    
    while(my $line = <IN>){
	chomp($line);
	my ($sampleName) = split(/\t/, $line);
	$samples->{$sampleName} = $line;
    }

    my $count = 0;
    
    if(scalar keys %{$samples} > 1){
	print HTML "\t\t\t\t\t<div class=\"row-fluid\">\n";
	print HTML "\t\t\t\t\t\t<div class=\"span12\" style=\"padding-top: 20px\">\n";
	print HTML "\t\t\t\t\t\t\t<div class=\"pull-left\">\n";
	print HTML "\t\t\t\t\t\t\t\t<strong>View Results</strong>:\n";
	print HTML "\t\t\t\t\t\t\t\t<select style=\"width:400px\" onchange=\"showSample(this.value)\">\n";
	
	foreach my $sample (keys %{$samples}){
	    if($count == 0){
		print HTML "\t\t\t\t\t\t\t\t\t<option value=\"$sample\" selected=\"selected\">$sample</option>\n";
		$count++;
	    }
	    else{
		print HTML "\t\t\t\t\t\t\t\t\t<option value=\"$sample\">$sample</option>\n";
	    }
	}
	
	print HTML "\t\t\t\t\t\t\t\t</select>\n";
	print HTML "\t\t\t\t\t\t\t</div>\n";
	print HTML "\t\t\t\t\t\t</div>\n";
	print HTML "\t\t\t\t\t</div>\n";
	print HTML "\n";
    }

    $count = 0;

    foreach my $sample (keys %{$samples}){
	my ($name, $assembledReads, $coverage, $consensusLength1, $consensusLength2, $numContig1, $numContig2, $largestContig, $n501, $n502, $n901, $n902, $n951, $n952, $commandsFile, $fastaFile, $wigFile, $aceFile, $statFile, $logFile) = split(/\t/, $samples->{$sample});

	if($count == 0){
	    print HTML "\t\t\t\t\t<div name=\"$sample\" id=\"$sample\" class=\"row-fluid sample\">\n";
	    $count++;
	}
	else{
	    print HTML "\t\t\t\t\t<div name=\"$sample\" id=\"$sample\" class=\"row-fluid sample\" style=\"display:none\">\n";
	}

	print HTML "\t\t\t\t\t\t<div class=\"span12\'\">\n";
	print HTML "\t\t\t\t\t\t\t<div style=\"padding:0 19px\">\n";
	print HTML "\n";

	#downloads
	print HTML "\t\t\t\t\t\t\t\t<div class=\"accordion-group\">\n";
	print HTML "\t\t\t\t\t\t\t\t\t<div class=\"accordion-heading\">\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t<a class=\"accordion-toggle\" data-toggle=\"collapse\" href=\"\#collapseDownload$sample\">Downloads</a>\n";
	print HTML "\t\t\t\t\t\t\t\t\t</div>\n";
	print HTML "\t\t\t\t\t\t\t\t\t<div id=\"collapseDownload$sample\" class=\"accordion-body collapse in\">\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t<div class=\"accordion-inner\">\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t<div class=\"row-fluid\" id=\"downloadPane\">\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t<div class=\"row-fluid\">\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t<div class=\"span12\">\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t<p>Download all your assembly result files.</p>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t<div style=\"padding-top:5px\">\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t<p><a target=\"_blank\" href=\"$fastaFile\"><i class=\"icon-download\"></i> Assembled Contigs (FASTA)</a>&nbsp;&nbsp;|&nbsp;&nbsp;<a target=\"_blank\" href=\"$wigFile\"><i class=\"icon-download\"></i> Coverage Information (WIG)</a>&nbsp;&nbsp;|&nbsp;&nbsp;<a target=\"_blank\" href=\"$statFile\"><i class=\"icon-download\"></i> Assembly Statistics (TXT)</a>&nbsp;&nbsp;|&nbsp;&nbsp;<a target=\"_blank\" href=\"$logFile\"><i class=\"icon-download\"></i> MIRA Log (TXT)&nbsp;&nbsp;|&nbsp;&nbsp;<a target=\"_blank\" href=\"$commandsFile\"><i class=\"icon-download\"></i> MIRA Commands (TXT)</a></p>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t</div>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t</div>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t</div>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t</div>\n";
	print HTML "\t\t\t\t\t\t\t\t\t</div>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t</div>\n";
	print HTML "\t\t\t\t\t\t\t\t</div>\n";
	print HTML "\n";

	#assembly stats
	print HTML "\t\t\t\t\t\t\t\t<div style=\"padding-top:10px\">\n";
	print HTML "\t\t\t\t\t\t\t\t\t<div class=\"accordion-group\">\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t<div class=\"accordion-heading\">\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t<a class=\"accordion-toggle\" data-toggle=\"collapse\" href=\"\#collapseStats$sample\">Assembly Statistics</a>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t</div>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t<div id=\"collapseStats$sample\" class=\"accordion-body collapse in\">\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t<div class=\"accordion-inner\">\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t<div class=\"row-fluid\" id=\"statPane\">\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t<div class=\"row-fluid\">\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t<div class=\"span12\">\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t<p>Assembly summary statistics for $sample.</p>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t<div style=\"padding-top:5px\">\n";
	print HTML "\n";

	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t<table class=\"table table-condensed\">\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t<tr><th>Parameter</th><th>Value</th></tr>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t<tr><td>Reference</td><td>".$ENV{"PLUGINCONFIG__AGENOME"}."</td></tr>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t<tr><td>MIRA Version</td><td>".$ENV{"PLUGINCONFIG__MIRAVERSION"}."</td></tr>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t<tr><td>Assembly Type</td><td>".$ENV{"PLUGINCONFIG__TYPE"}."</td></tr>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t<tr><td>Fraction of Reads Used</td><td>".$ENV{"PLUGINCONFIG__FRACTION_OF_READS"}."</td></tr>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t<tr><td>RAM</td><td>".$ENV{"PLUGINCONFIG__RAM"}."</td></tr>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t<tr><td>Barcode Minimum Read Cutoff</td><td>".$ENV{"PLUGINCONFIG__MIN_READS"}."</td></tr>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t</table>\n";
	print HTML "\n";

	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t<table class=\"table table-condensed\">\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t<tr><th>Metric</th><th>Large Contigs</th><th>All Contigs</th></tr>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t<tr><td>Assembled Reads</td><td colspan=\"2\">$assembledReads</td></tr>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t<tr><td>Coverage</td><td colspan=\"2\">$coverage</td></tr>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t<tr><td>Consensus Length</td><td>$consensusLength1</td><td>$consensusLength2</td></tr>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t<tr><td>Number of Contigs</td><td>$numContig1</td><td>$numContig2</td></tr>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t<tr><td>Largest Contig</td><td colspan=\"2\">$largestContig</td></tr>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t<tr><td>N50</td><td>$n501</td><td>$n502</td></tr>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t<tr><td>N90</td><td>$n901</td><td>$n902</td></tr>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t<tr><td>N95</td><td>$n951</td><td>$n952</td></tr>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t</table>\n";
	print HTML "\n";

	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t</div>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t\t</div>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t\t</div>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t\t</div>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t\t</div>\n";
	print HTML "\t\t\t\t\t\t\t\t\t\t</div>\n";
	print HTML "\t\t\t\t\t\t\t\t\t</div>\n";
	print HTML "\t\t\t\t\t\t\t\t</div>\n";
	print HTML "\n";

	print HTML "\t\t\t\t\t\t\t</div>\n";
	print HTML "\t\t\t\t\t\t</div>\n";
	print HTML "\t\t\t\t\t</div>\n";
    }

    #close it up
    print HTML "\t\t\t\t</div>\n";
    print HTML "\t\t\t</div>\n";
    print HTML "\t\t</div>\n";
    print HTML "\t</body>\n";
    print HTML "</html>\n";

    close(HTML);
}
