#!/usr/bin/perl


#################################################################################
###
### wraps IntaRNA for the webserver, producing automatic functional enrichment
### and regions plots
###
#################################################################################

use warnings;
use strict;

use Getopt::Long;
##added
unless($ARGV[0]){die "Usage IntaRNA_wrapper.pl -bindseq Yfr1.fasta -ncbiacc NC_005072 -ntupstream 80 -ntdownstream 30 -region 5utr
\n\t-n\t\tNCBI accession number
\t-tarseq\t\tsRNA targets to search  (default = target_RNAs.fa)
\t-bindseq\tsRNA fasta sequence (only 1 entry) (default = binding_RNAs.fa)
\t-e\t\tEnergy threshold (default = 0)
\t-subopres\tNumber of sub optimal results (default = 0)
\t-bpinseed\tBasepairs in seed (default=7)
\t-u\t\tNumber of unpaired sites in seed (default = 0)
\t-temp\t\tTemperature (default = 37)
\t-w\t\tWindow size (default = 140)
\t-m\t\tMax bp distance (default = 70)
\t-ntupstream\tNumber of nt upstream of start site (default = 30)
\t-ntdownstream\tNumber of nt downstream of start site (default = 30)
\t-r\t\tRegion 5utr, 3utr, cds (default = 5utr)
\t-seedncrna\tSeed region of sRNA (default ='')
\t-a\t\tCheck all replicons (chromosome + plasmid) Yes = 1, No = 0 (default = 1)
\t-g\tUse local GenBank file for IntaRNA? (default = N)
\n";}
## added -g  option PHD
## edit 1.0.4; added PATH variables
my $PATH_HOME="/software/IntaRNA1.0.4";
my $PATH_INTARNA_BIN="/software/intarna-1.2.5/src";
my $PATH_COPRARNA="/software/CopraRNA/1.2.9";
my $PATH_R="/usr/bin";

# -o will always be used on the IntaRNA call

my $ncbi_acc = ''; # NC_000913
my $tar_seqs = 'target_RNAs.fa'; # -t in intarna
my $bind_seqs = 'binding_RNAs.fa'; # -m in intarna
my $energy_threshold =  0.0; # -v
my $subopt_resultsnum = 0; # -s
my $bp_in_seed = 7; # -p
#my $unpair_seed1 = 0; # -u1
#my $unpair_seed2 = 0; # -u2
my $unpair_seed_both = 0; # -u
my $temprature = 37.0; # -T
my $window_size = 140; # -w
my $max_bp_dist = 70; # -L
my $nt_upstream = 30;
my $nt_downstream = 30;
my $region = "5utr";
my $seed_ncrna = "";
my $all_replicons = "";#1; # checkbox should be checked as standard option
my $local = "N"; ## added
my $kegg2refseqLoc = "$PATH_COPRARNA/kegg2refseqnew.csv";

# i=integet, f=float, s=string
GetOptions(
           # ncbi params
           'ncbiacc:s'		=> \$ncbi_acc, # empty unless specified
           'tarseq:s'		=> \$tar_seqs, # target sequences fasta file
           'bindseq:s'		=> \$bind_seqs,
           'enrgthr:f'		=> \$energy_threshold,	
	   'subopres:i'		=> \$subopt_resultsnum,
           'bpinseed:i'		=> \$bp_in_seed,
#	   'unpairseed1:i'	=> \$unpair_seed1, # targets
#	   'unpairseed2:i'	=> \$unpair_seed2, # binding RNAs
	   'unpairseedboth:i'	=> \$unpair_seed_both, # binding RNAs
           'temp:f'		=> \$temprature,                      
           'windowsize:i'	=> \$window_size,
           'maxbpdist:i'	=> \$max_bp_dist,
           'ntupstream:i'	=> \$nt_upstream,
	   'ntdownstream:i'	=> \$nt_downstream,
           'region:s'		=> \$region, # one of "5utr", "3utr", "cds"
           'seedncrna:s'	=> \$seed_ncrna, # example "1,10" for positions 1-10 in ncRNA
           'allreplicons:i'	=> \$all_replicons, # example either 1 or 0. If 1 then execute ncbi IntaRNA on all replicons else not
	   'genbank:s'          => \$local, #Yes/No
);

#override $unpair_seed_both
#$unpair_seed_both = -1;

my $switch = 0;
my $switch2 = 0;

# If user specified sRNA file is not called binding_RNAs.fa, copy it to that file name
unless($bind_seqs eq "binding_RNAs.fa"){
	`cp $bind_seqs binding_RNAs.fa`;
	$bind_seqs="binding_RNAs.fa";
}

my @files = <*>;
foreach(@files) {
    if ($_ eq $bind_seqs) { $switch = 1; }
    if ($_ eq $tar_seqs) { $switch2 = 1; }
   
}#print "Switch=[$switch]\n";

# add exception, only allow one sRNA on ncbi IntaRNA
if ($ncbi_acc) {
    system "grep -c '>' binding_RNAs.fa > binding_count.txt";
    open (MYDATA, "binding_count.txt") or die("Error: cannot open file binding_count.txt at IntaRNA_wrapper.pl\n");
        my @count = <MYDATA>;
    close MYDATA;
    chomp $count[0];
    #print "Count=[$count[0]]\n";
    system "rm binding_count.txt";
    if ($count[0] > 1) { die("only one sRNA sequence allowed in NCBI IntaRNA. $count[0] sequences are present.\n"); }
}

#print "Switch=[$switch]\n";
unless ($switch) { die("Error:No binding sequences supplied!\n"); }

unless ($ncbi_acc) {
    unless ($switch2) { die("Error:No target sequences supplied!\n"); }
}


if ($window_size < $max_bp_dist) { die("maxbpdist must be <= windowsize\n"); }


# if ncbi accession number is specified then run the ncbi intarna
if ($ncbi_acc) {
    my $replicon_count = 0;
    my $kegg2refseqLine = "";

    system "> termClusterReport.txt"; ## edit 1.0.2 make sure file exists in ncbi mode

    unless ($tar_seqs eq 'target_RNAs.fa') {
        print "overriding $tar_seqs with target_RNAs.fa\n";
        $tar_seqs = 'target_RNAs.fa';
    }
    $subopt_resultsnum = 0; # suboptimals not allowed for ncbi run

    #system "$PATH_HOME/get_refseq_from_ftp.sh $ncbi_acc" unless ($all_replicons);
    ## added control to allow flag to use local genbank file
    unless($local =~ /[Yy]/) {    
	system "$PATH_HOME/get_refseq_from_ftp.pl $ncbi_acc" unless ($all_replicons);
    }
    # add option for getting further replicons here if checkbox is on.
    if ($all_replicons) {
        open (MYDATA, $kegg2refseqLoc) or die("Error: cannot open file $kegg2refseqLoc at IntaRNA_wrapper.pl\n");
            my @kegg2refseq = <MYDATA>;
        close MYDATA;
        foreach(@kegg2refseq) {
            if ($_ =~ m/$ncbi_acc/) {
                $kegg2refseqLine = $_;
                last;
            }
        }
        chomp $kegg2refseqLine;
        my @split = split(/\s+|\t+/,$kegg2refseqLine);
        shift @split;
        foreach(@split) {
            system "$PATH_HOME/get_refseq_from_ftp.pl $_";
            print "$PATH_HOME/get_refseq_from_ftp.pl $_\n";
        }
        $replicon_count = scalar(@split); 
    }

    # manually added for karen kong eco plasmid NC_019063.gb ## edit 1.0.3
    if($ncbi_acc eq "NC_019063") {
        system "cp $PATH_HOME/NC_019063.gb .";
    }
    # manually added for joke lambrecht 
    if($ncbi_acc eq "NC_006882") {
        system "cp $PATH_HOME/NC_006882.gb .";
    }

    # check if RefSeq download worked otherwise retry
    my $repeats = 0;

    while($repeats <= 10) {
        my @RefSeq = <*gb>;
        $repeats++;
        unless($RefSeq[($replicon_count-1)]) {
            die("Issues downloading all RefSeq records from the NCBI. Please resubmit your request.\n");
            #print "unless\n";
            #sleep 30;
            #system "/scratch/rna/bisge001/Software/IntaRNA_wrapper/get_refseq_from_ftp.sh $ncbi_acc";
        #} else {
            #print "last\n";
         #   last;
        }
    }


    ## edit 1.0.1 adding the omission of the CONTIG line
    my @RefSeqFiles = ();
    @RefSeqFiles = <*gb>;

    foreach (@RefSeqFiles) {
        system "sed -i '/^CONTIG/d' $_"; ## d stands for delete
    }
    ## end CONTIG issue fix

    my @RefSeq = <*gb>;
    foreach (@RefSeq) {    
        system "$PATH_HOME/parse_region_from_genome.pl $_ $nt_upstream $nt_downstream $region >> target_RNAs.fa";
    }
}

# check if seq sizes are appropriate for windowsize and maxbplength
unless ($ncbi_acc) {
open (MYDATA, $tar_seqs) or die("Error: cannot open file $tar_seqs at IntaRNA_wrapper.pl\n");
    my @tar_seq_lines = <MYDATA>;
close MYDATA;

my $header = $tar_seq_lines[0];
chomp $header;

my $curr_seq = "";

for (my $i=1;$i<scalar(@tar_seq_lines);$i++) {
    unless($tar_seq_lines[$i] =~ m/>/) {
        my $curr_line = $tar_seq_lines[$i];
        chomp $curr_line;
        $curr_seq = $curr_seq . $curr_line;
    } else {
        if (length($curr_seq) < $window_size) {
            die("Sequence with name $header is to short with respect to the window size of $window_size. Remove it and rerun.\n"); 
        }
        $header = $tar_seq_lines[$i];
        $curr_seq = "";
        chomp $header;
    }
}
    # repeat for last sequence
    if (length($curr_seq) < $window_size) {
        die("Sequence with name $header is to short with respect to the window size of $window_size. Remove it and rerun.\n"); 
    }
}
# endcheck



# adjust $window_size and $max_bp_dist for sequences smaller than 140 nt
if ($ncbi_acc) {
    if (($nt_upstream+$nt_downstream) < $window_size) {
        $window_size = $nt_upstream+$nt_downstream;
        $max_bp_dist = $nt_upstream+$nt_downstream;
        print "overriding windowize and maximum bp distance window(-w)=$window_size, maxbp(-L)=$max_bp_dist\n";
    }
}


# run intarna
unless ($seed_ncrna) {
    system "$PATH_INTARNA_BIN/IntaRNA -v $energy_threshold -s $subopt_resultsnum -o -p $bp_in_seed -u $unpair_seed_both -T $temprature -w $window_size -L $max_bp_dist -t $tar_seqs -m $bind_seqs > predictions.intarna";
} else {
    system "$PATH_INTARNA_BIN/IntaRNA -v $energy_threshold -s $subopt_resultsnum -o -p $bp_in_seed -u $unpair_seed_both -T $temprature -w $window_size -L $max_bp_dist -t $tar_seqs -m $bind_seqs -f $seed_ncrna > predictions.intarna";
}

# process intarna output
system "$PATH_HOME/produce_semicolon_sep_results_from_intarna_out.pl --intarna-out-file predictions.intarna > predictions.csv";

# sort intarna output
system "$PATH_HOME/sort_intarna_csv_results.pl --intarna-csv-file predictions.csv --column 17 > predictions.sorted.csv";

# omit WARNING for 'N' from IntaRNA raw output. fixes interaction parsing warnings.
system "sed -i '/WARNING/{N;d;}' predictions.intarna"; ## d stands for delete N also deletes one line after. add another N to delete one more ##edit 1.0.1 

# add interactions to table
system "$PATH_HOME/add_interactions.pl predictions.sorted.csv predictions.intarna";

# if ncbi option is on then also add annatations and Entraz Gene IDs and gene name 
if ($ncbi_acc) {
    system "$PATH_HOME/add_GI_genename_annotation.pl";
}

# maka a functional enrichment for the top 50
print "Outside termClusterReport.pl\nncbi_acc[$ncbi_acc]\n";
if ($ncbi_acc) {
	my @enrich_lines=(-1,-1);
	#print "[$enrich_lines[0]]\n";
	my $counter=0;
	
	until($enrich_lines[0]  =~ /^Annotation/  || $counter == 10 ){## added PHD
    		print "$PATH_HOME/termClusterReport.pl intarna_websrv_table_ncbi.csv\n";
		system "$PATH_HOME/termClusterReport.pl intarna_websrv_table_ncbi.csv";
		open(MYDATA, "termClusterReport.txt") or die("Error: cannot open file termClusterReport.txt\n");
        	@enrich_lines = <MYDATA>;
    		close MYDATA;
		$counter++;
		sleep(5);
    	}
	
	unless ($enrich_lines[0]  =~ /^Annotation/ ) {
        	system "echo -e 'If you are reading this, then your prediction did not return an enrichment, your organism is not in the DAVID database\nor the DAVID webservice is/was termporarily down. You can either rerun your IntaRNA\nprediction or create your enrichment manually at the DAVID homepage.' > termClusterReport.txt";
	}
    
    	system "mv intarna_websrv_table_ncbi.csv intarna_websrv_table.csv";
}

system "head predictions.sorted.csv -n 101 > IntaRNA_result.csv";

# write file for regions plot script
open(WRITETO, '>intarna_pars.txt');
    print WRITETO "$region\n$nt_upstream\n$nt_downstream\n";
close(WRITETO);

# do length normalization for ncbiacc cds option
if ($ncbi_acc and ($region eq "cds")) {
    system "$PATH_HOME/length_normalization.pl";
}


# add pvalues and qvalues
if ($ncbi_acc) {
    system "$PATH_R/R --slave -f $PATH_HOME/add_pv_qv.r"; 
    system "paste PVQV.csv intarna_websrv_table.csv -d ';' > intarna_websrv_table_temp.csv";
    system "mv intarna_websrv_table_temp.csv intarna_websrv_table.csv";
}

system "head intarna_websrv_table.csv -n 100 > intarna_websrv_table_truncated.csv";

# make regions plots

if ($ncbi_acc) { ## edit 1.0.1
    system "$PATH_R/R --slave -f $PATH_HOME/plotting_script_intaRNA.r";

   # make heatmap
   system "$PATH_R/R --slave -f $PATH_HOME/extract_enriched_IntaRNA.r"; ## edit 1.0.2  ## edit 1.0.3 changed R code for bugfix
   system "cp $PATH_HOME/intarna_heatmap.html ."; ## edit 1.0.2 ## edit 1.0.3 bugfix steffen html
   system "$PATH_HOME/make_heatmap_json.pl enrichment.txt"; ## edit 1.0.2
   system "cp $PATH_HOME/index-thumb.html ."; ## edit 1.0.2
   system "cp $PATH_HOME/index-pdf.html ."; ## edit 1.0.3
   #$PATH_HOME/phantomjs
   system "$PATH_HOME/phantomjs $PATH_HOME/rasterize.js ./index-thumb.html enriched_heatmap_big.png";
   system "$PATH_HOME/phantomjs $PATH_HOME/rasterize.js ./index-pdf.html enriched_heatmap_big.pdf"; ## edit 1.0.3
   system "rm index-thumb.html"; ## edit 1.0.2
   system "rm index-pdf.html"; ## edit 1.0.3

    # thumbnails

    system "/usr/bin/convert -size 170x170 -resize 170x170 -flatten -rotate 90 sRNA_regions.ps thumbnail_sRNA.png";
    system "/usr/bin/convert -size 170x170 -resize 170x170 -flatten -rotate 90 mRNA_regions.ps thumbnail_mRNA.png";

    # blow up images

    system "/usr/bin/convert -density '300' -resize '700' -flatten -rotate 90 sRNA_regions.ps sRNA_regions.png";
    system "/usr/bin/convert -density '300' -resize '700' -flatten -rotate 90 mRNA_regions.ps mRNA_regions.png";
}

# write a README.txt

# remove not needed files from archive

#system "rm predictions.intarna";
#system "rm *gb PVQV.csv" if($ncbi_acc);
system "rm PVQV.csv" if($ncbi_acc); ##phd
