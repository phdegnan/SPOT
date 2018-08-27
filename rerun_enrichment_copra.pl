#!/usr/bin/perl


unless(@ARGV){die "\nUsage: $0 ncRNA_hIntaRNA.csv \n\n";}

print "\nRunning $0 on $ARGV[0]\n";

my $versiondirectory = "/software/CopraRNA/1.2.9/";


system $versiondirectory . "termClusterReport.pl $ARGV[0] 0.01"; #result goes into termClusterReport.txt


$check=`grep 'Annotation' termClusterReport.txt`;
if($check !~ /Annotation/){
	system "echo -e 'If you are reading this, then your prediction did not return an enrichment, your organism of interest is not in the DAVID database\nor the DAVID webservice is/was termporarily down. You can either rerun your CopraRNA\nprediction or create your enrichment manually at the DAVID homepage.' > termClusterReport.txt";
	print "Predicted sRNA-interacting mRNAs did not return valid functional enrichment results.\n";
	print "Check and rerun: $0 $ARGV[0]\n"; 
	print "However, there may not be any functional enrichment.\n"; 
}else{

## add enrichment visualization ## edit 1.2.5
print "cp $versiondirectory" . "copra_heatmap.html .\n";
system "cp $versiondirectory" . "copra_heatmap.html ."; ## edit 1.2.5 ## edit 1.2.7 (edited html file)

print "/usr/bin/R --slave -f " . $versiondirectory . "extract_functional_enriched.r\n"; #requires CopraRNA_result_all.csv
system "/usr/bin/R --slave -f " . $versiondirectory . "extract_functional_enriched.r"; ## edit 1.2.5 ## edit 1.2.7 (edited R code)

print $versiondirectory . "make_heatmap_json.pl enrichment.txt\n";
system $versiondirectory . "make_heatmap_json.pl enrichment.txt"; ##edit 1.2.5
system "cp $versiondirectory" . "index-thumb.html ."; ## edit 1.2.5
system "cp $versiondirectory" . "index-pdf.html ."; ## edit 1.2.6
#$versiondirectory . 
print "/usr/local/bin/phantomjs " . $versiondirectory . "rasterize.js " . "./index-thumb.html enriched_heatmap_big.png\n"; 
system "/usr/local/bin/phantomjs " . $versiondirectory . "rasterize.js " . "./index-thumb.html enriched_heatmap_big.png"; ## edit 1.2.5
system "/usr/local/bin/phantomjs " . $versiondirectory . "rasterize.js " . "./index-pdf.html enriched_heatmap_big.pdf"; ## edit 1.2.6
system "rm index-thumb.html"; ## edit 1.2.5
system "rm index-pdf.html"; ## edit 1.2.6

}
