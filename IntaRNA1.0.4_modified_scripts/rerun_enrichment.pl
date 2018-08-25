#!/usr/bin/perl

unless(@ARGV){die "\nUsage: $0 intarna_websrv_table_ncbi.csv \n\n";}

#$ARGV[0]

my $PATH_HOME="/data/software/IntaRNA1.0.4";
my $PATH_R="/usr/bin";

system "$PATH_HOME/termClusterReport.pl $ARGV[0]";
system "$PATH_R/R --slave -f $PATH_HOME/extract_enriched_IntaRNA.r"; ## edit 1.0.2  ## edit 1.0.3 changed R code for bugfix
system "cp $PATH_HOME/intarna_heatmap.html ."; ## edit 1.0.2 ## edit 1.0.3 bugfix steffen html
system "$PATH_HOME/make_heatmap_json.pl enrichment.txt"; ## edit 1.0.2
system "cp $PATH_HOME/index-thumb.html ."; ## edit 1.0.2
system "cp $PATH_HOME/index-pdf.html ."; ## edit 1.0.3
system "/usr/local/bin/phantomjs $PATH_HOME/rasterize.js ./index-thumb.html enriched_heatmap_big.png";
system "/usr/local/bin/phantomjs $PATH_HOME/rasterize.js ./index-pdf.html enriched_heatmap_big.pdf"; ## edit 1.0.3
system "rm index-thumb.html"; ## edit 1.0.2
system "rm index-pdf.html"; ## edit 1.0.3
