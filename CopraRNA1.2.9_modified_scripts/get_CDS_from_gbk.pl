#!/usr/bin/perl

use strict;
use warnings;

use Bio::SeqIO;

# ./get_CDS_from_gbk.pl NC_003197.gb > stm.fas

# gets a refseq file and creates a fasta
# with all CDS AA sequences in the following format:
# >eco:b0001
# ISLLGQCPDSINAALICRGEKMSIAIMAGVLEARGHNVTVIDPVEKLL

# >KEGGID:locus_tag
# Aminoacid sequence

my $versiondirectory = "/data/software/CopraRNA/1.2.9/";

my $kegg2refseq = $versiondirectory . "kegg2refseqnew.csv";

open(MYDATA, $kegg2refseq) or die("Error: cannot open file $kegg2refseq'\n");
    my @kegg2refseq = <MYDATA>;
close MYDATA;



my $refid = $ARGV[0];
my $refidchopped = $refid;
chop $refidchopped;
chop $refidchopped;
chop $refidchopped;

my $keggcode = "";

foreach(@kegg2refseq) {
    if($_ =~ m/$refidchopped/) {
        my @split = split(/\s/,$_);
        $keggcode = $split[0];
        last;
    }
}

##PHD
if($keggcode eq ""){die "$refid not included in kegg2refseqnew.csv\nPossible plasmid replicon or update to csv file required -- Skipped!\n";}



my $seqin = Bio::SeqIO->new( -format => 'genbank', -file => $refid);

while( (my $seq = $seqin->next_seq()) ) {
    foreach my $sf ( $seq->get_SeqFeatures() ) {
        if( $sf->primary_tag eq 'CDS' ) {
            my $ltag = "";
            my $protein = "";
            if ($sf->has_tag("locus_tag") and $sf->has_tag("translation")) {
                my @ltaglist = $sf->get_tag_values("locus_tag");
                $ltag = $ltaglist[0];
                my @translationlist = $sf->get_tag_values("translation");
                $protein = $translationlist[0];
            
                chomp $ltag;
                print ">$keggcode:$ltag\n";
                chomp $protein;
                print "$protein\n";
         }
      }
   }
}

