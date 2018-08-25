#!/usr/bin/perl 
##
## 	Patrick Degnan
##	fnaptt2gbk.pl
##	Convert GenBank ptt and fna files into a gb file
##	Does not use BioPerl
## 	v2 of script updated 30 Aug 2017
##

use Getopt::Long;
unless(@ARGV){
	print "\nUsage $0 :
\t-f\tfasta of genome sequence (multi sequence fasta file will not work)
\t-t\trnt, ptt or pttp formatted annotation table
\t-r\tRefseq number (if not used in file name)
\t-g\tRetreive GeneIDs? Will only work for bonafide GI numbers in ptt file (default = N)\n
$0 -f NC_016810.fna -t NC_016810.ptt \n\n";	
	exit;
}

GetOptions("fasta=s"=>\$fasta,"table=s"=>\$table,"refseq=s"=>\$genome_accession
	   ,"geneids=s"=>\$gids);

#\t-p\tPrint GI #s? If  (default = Y)

#NC_016810.ptt
#$fasta="NC_016810.fna";
#$genome_num=1000;

if(defined($genome_accession)){
	# do nothing
}else{
	if($fasta =~/(N.*?)\.f.*/){
		$genome_accession=$1;
	}elsif($table =~/(N.*?)\..*/){
		$genome_accession=$1;
	}else{
		die "Cannot determine Refseq number for output GenBank file\n";
	}
}
$file = $genome_accession . ".gb";

open(OUT, ">$file") or die "Cannot open file: $file\n";
open(FA, "<$fasta") or die "Cannot open file: $fasta\n";

$/=">";
$first=<FA>;
while ($line=<FA>) {
 	if($line=~ /^\S+?\|\s(.+?)\,.+\n/){ # if a refseq like name ending w/ , complete genome
		$genome_name_spaces=$1;
  	}else{
  		$line=~ /^(.+)\n/;
  		$genome_name_spaces=$1;
  	}
  	$line=~s/\|/\t/g;
  	$line=~s/[\]\[\)\(]//g;
  
  	if($line=~ /(.+)\n/){
  		$id=$1;
  		$line=~ s/$id//;
  		$line=~ s/[\>\n]//g;
  		#print "[$line]\n";
  		$genome_sequence=$line;
		$genome_length=length($line);
  }

}
close(FA);

$/="\n";

#$number_of_proteins=382;
#$genome_name_nospaces=$genome_name_spaces;
#$genome_name_nospaces=~s/\s/\_/g;
#$strain="Salmonella enterica subsp. enterica serovar Typhimurium str. SL1344";

$today=`date`;	#retrieve today's date
##Wed Aug  6 09:07:32 EDT 2003 to 080603
($day,$month,$numba,$time,$zone,$year)=split(" ",$today);	#break it down into its components
#make hash of Months and numerical counterparts
%months=qw(Jan 01 Feb 02 Mar 03 Apr 04 May 05 Jun 06 Jul 07 Aug 08 Sep 09 Oct 10 Nov 11 Dec 12);
if($numba < 10){$numba="0" . $numba;}		#if a single digit number place a 0 in front of it
$today=$numba . "-" . uc($month) . "-" . $year;		#concatenate new elements together

#Add PFAM, TIGRFAM annotations for future versions when present in pttp files
$pfam="";
$tigrfam="";
#$table="SL1344_sRNAs_pttp.txt";
open(TBL, "<$table") or die "Cannot open file: $table\n";
$head1=<TBL>;
$head1=~/.*?\t*(.+)$/;
$strain=$1;
$head2=<TBL>;
$head3=<TBL>;

$last=-1;  $annotations=""; $inregion="N"; $number_of_proteins=0;
while($l=<TBL>){
	chomp($l);
	@row=split(/\t/,$l);
	$gene="";$geneid="";$protid="";$printgeneid="";
	if($row[0]=~/(\d+)\.\.(\d+)/){ ## ptt 
		$start = $1 + 0;
		$stop  = $2 + 0;
		$locus=$row[5]; 
		$product=$row[8];
		$strand=$row[1];
		$type="CDS";
		$pid = $row[3];
		if($row[4] =~ /[\d\w]+/){
			$gene = "
                     /gene=\"$row[4]\"";
		}
	}elsif($row[3]=~/[\+\-]/){ ## pttp w/ first column containing gnum #
		$locus=$row[7]; 
		$start = $row[1] + 0;
		$stop  = $row[2] + 0;
		$product=$row[10];
		$strand=$row[3];
		$type=$row[11];
		$pid = $row[5];
		if($row[6] =~ /[\d\w]+/){
			$gene = "
                     /gene=\"$row[6]\"";
		}
	}elsif($row[2]=~/[\+\-]/){ ## pttp 
		$locus=$row[6]; 
		$start = $row[0] + 0;
		$stop  = $row[1] + 0;
		$product=$row[9];
		$strand=$row[2];	
		$type=$row[10];
		$pid = $row[4];
		if($row[5] =~ /[\d\w]+/){
			$gene = "
                     /gene=\"$row[5]\"";
		}
	}else{
		die "Input table does not conform to known format.\n";
	}	
	
	if($start < $stop){#determine order
		$begin=$start;
		$end=$stop;
	}else{
		$begin=$stop;
		$end=$start;
	}
	
	if($strand eq "+"){$coord="$begin..$end";}
	else{$coord="complement($begin..$end)";}
	@returna="";
	if($gids =~ /[Yy]/ && $pid =~ /\d+/){
		#print "[$pid]";
		$protid=`efetch -db protein -id $pid -format acc`;
		chomp($protid);#returns refseq protein accession
		$protid=~s/\.\d+//;
		#print "[$protid]";
		system("epost -db protein -id $protid >temp_a.txt");		
		until($returna[0] == 0 || sleep(60) ){
			@returna=system("elink -target gene -name protein_gene < temp_a.txt > temp_b.txt");
			print "[$pid][$returna[0]]\n";
		}
		$geneid=`efetch -db gene -format uid < temp_b.txt`;
		chomp($geneid);
		#print "[$geneid]\n";
		# all in one encounters difficulty with elink errors
		#system("epost -db protein -id $protid \| elink -target gene -name protein_gene \| efetch -db gene -format uid");
		if($geneid =~/\d+/){
			$printgeneid = "
                     /db_xref=\"GeneID=$geneid\"";
		}else{
			print "[$pid][$protid][$geneid]\n";
		}
	
	}
	#$pfam$tigrfam
	if($type =~ /CDS/){
		
		$annotations.="     gene            $coord$gene
                     /locus_tag=\"$locus\"$printgeneid
     CDS             $coord$gene
                     /locus_tag=\"$locus\"
                     /codon_start=1
                     /transl_table=11
                     /product=\"$product\"$printgeneid
                     /db_xref=\"GI:$pid\"
                     /note=\"$product\"
";	}else{
		
		if($row->[7] =~ /tRNA/ ){
			$type="tRNA ";
		}elsif($row->[7] =~ /rRNA/){
			$type="rRNA ";
		}else{
			$type="ncRNA";
		}
		$annotations.="     gene            $coord$gene
                     /locus_tag=\"$locus\"$printgeneid
     $type           $coord$gene
                     /locus_tag=\"$locus\"$printgeneid
                     /product=\"$product\"
                     /note=\"$product\"
";	
	
	}


	$last=$start;
	$number_of_proteins++;
}


$fragment_length = $genome_length;


## trim down accessions that are too long e.g., B. theta str 3731
$title=length($genome_accession) + length($fragment_length);
if($title >=28){
	$sub_len=$title - 27;
	$header=substr($genome_accession,0,(length($genome_accession)-$sublen-1));	
}else{
	$header=$genome_accession
}	
print OUT "LOCUS       $header";
## add spaces:
$temp= 28 - length($header) - length($fragment_length);
for $i (1..$temp){
	print OUT " ";
}
print OUT "$fragment_length bp    DNA     linear   PHG $today\n";
#LOCUS       KJ830768               11019 bp    DNA     linear   PHG 22-OCT-2014

print OUT "DEFINITION  $genome_name_spaces
ACCESSION   $genome_accession
REFERENCE   1  (bases 1 to $fragment_length)
  AUTHORS   Anonymous, A.
  TITLE     Direct Submission
  JOURNAL   Submitted ($today) HMP

FEATURES             Location/Qualifiers
     source          1..$fragment_length
                     /organism=\"$genome_name_spaces\"
                     /mol_type=\"genomic DNA\"
                     /strain=\"$genome_name_spaces\"
";

## PRINT ACCUMULATED ANNOTATIONS
print OUT $annotations;

#$copy=$genome_sequence;
#$as = $copy =~ s/[Aa]//g;
#$ts = $copy =~ s/[Tt]//g;
#$cs = $copy =~ s/[Cc]//g;
#$gs = $copy =~ s/[Gg]//g;
#$other = length($other) + 0;

## PRINT SEQUENCE

print OUT "ORIGIN      
";

$genome_sequence=~s/(.{10})/$1 /g;
#print "$genome_sequence\n";
$count=1; $extracounter=66;
while( $genome_sequence=~ /(.{66})/g ){ # || $genome_sequence=~ /(.+?)/
	$extra_spaces = 9 - length(scalar($count));
	for $i (1..$extra_spaces){
		print OUT " ";
	}	
	print OUT "$count $1\n";
	
	$extracounter+=66;
	$count+=60;
}

$length= length($genome_sequence) - ($extracounter - 66)  ;
#print length($genome_sequence) . "\t$extracounter\t$length\n";
$last=substr($genome_sequence,$extracounter-66,$length); 


#6260361
$extra_spaces = 9 - length(scalar($count));
for $i (1..$extra_spaces){
	print OUT " ";
}	
print OUT "$count $last\n//\n";

#$genome_sequence


__END__

=begin

=end
=cut


	unless(exists($PFAM{$cds})){$PFAM{$cds}="\n                     /db_xref=\"PFAM:$pfam\"";}
	else{$PFAM{$cds}=$PFAM{$cds} . "\n                     /db_xref=\"PFAM:$pfam\"";}
	unless(exists($PFAM_DESCR{$cds})){$PFAM_DESCR{$cds}=$pfam_descr;}
	else{$PFAM_DESCR{$cds}=$PFAM_DESCR{$cds} . ", $pfam_descr";}

	unless(exists($TIGRFAM{$cds})){$TIGRFAM{$cds}="\n                     /db_xref=\"TIGRFAM:$tigrfam\"";}
	else{$TIGRFAM{$cds}=$TIGRFAM{$cds} . "\n                     /db_xref=\"TIGRFAM:$tigrfam\"";}
	unless(exists($TIGRFAM_DESCR{$cds})){$TIGRFAM_DESCR{$cds}=$tigrfam_descr;}
	else{$TIGRFAM_DESCR{$cds}=$TIGRFAM_DESCR{$cds} . ", $tigrfam_descr";}
