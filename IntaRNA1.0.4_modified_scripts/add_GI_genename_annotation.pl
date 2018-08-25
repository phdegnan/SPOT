#!/usr/bin/perl

use strict;
use warnings;

use Bio::SeqIO;

open(MYDATA, "intarna_websrv_table.csv") or die("Error: cannot open file intarna_websrv_table.csv\n");
    my @intarna_websrv_lines = <MYDATA>;
close MYDATA;

open(WRITETABLE, ">intarna_websrv_table_ncbi.csv");


my @gbks = <*gb>;
my %ltagToRestHash = ();

foreach(@gbks) {
my $gbkFile = $_;
chomp $gbkFile;

my $in  = Bio::SeqIO->new(-file => $gbkFile , '-format' => 'genbank');
    while ( my $seq = $in->next_seq() ) {
        foreach my $sf ( $seq->get_SeqFeatures() ) {
          if( $sf->primary_tag eq 'CDS' ) {
            my $ltag = "";
            my $annotation = "";
            my $GID = "";
            my $geneName = "";
            my $product = ""; # for hypothetical protein exception
            if ($sf->has_tag("locus_tag")) {
                my @ltaglist = $sf->get_tag_values("locus_tag");
                $ltag = $ltaglist[0];
                # add exception for "hypothetical protein" annotation
                if ($sf->has_tag("product")) {
                    my @products = $sf->get_tag_values("product");
                    $product = $products[0];
                }

                if ($sf->has_tag("product") and (not $product =~ m/hypothetical protein/) ) {
                    my @prodlist = $sf->get_tag_values("product");
                    $annotation = $prodlist[0];
                }
                elsif ($sf->has_tag("note")) {
                    my @notelist = $sf->get_tag_values("note");
                    $annotation = $notelist[0];
                }
                elsif ($sf->has_tag("function")) {
                    my @funclist = $sf->get_tag_values("function");
                    $annotation = $funclist[0];
                }
                if ($sf->has_tag("gene")) {
                    my @geneNameList = $sf->get_tag_values("gene");
                    $geneName = $geneNameList[0];
                }
                my @GIDlist = $sf->get_tag_values("db_xref") if ($sf->has_tag("db_xref"));
                foreach (@GIDlist) {
                    $GID = $_ if ($_ =~ m/GeneID:/);
                }
		if($GID){
                	$GID =~ s/GeneID://g;
                }else{ ## added PHD
			foreach (@GIDlist) {
                    		$GID = $_ if ($_ =~ m/GI:/);
               		}
			$GID =~ s/GI://g;
		}
		
		$annotation =~s/,//g;
                $annotation =~s/;//g;
                my $printLine = join(";",$geneName,$GID,$annotation); 
                $ltagToRestHash{$ltag} = $printLine;
            }
         }
    }
} 
}

foreach my $line (@intarna_websrv_lines) {
    chomp $line;
    print WRITETABLE $line;
    my @split = split(/;/,$line);
    my $ltag = $split[0];
    if(exists $ltagToRestHash{$ltag}) {
        print WRITETABLE ";" . $ltagToRestHash{$ltag} . "\n";
    } else {
        print WRITETABLE ";;;\n";
    }
}









