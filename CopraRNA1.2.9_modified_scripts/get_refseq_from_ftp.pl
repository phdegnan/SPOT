#!/usr/bin/perl

unless(@ARGV){die "\nUsage: get_refseq_from_ftp.pl NC_XXXXXXXX\n\n";}

## 25 Aug 2017 Edited IntaRNA copy. 

$text=`efetch -db nucleotide -id $ARGV[0] -mode text |head -2000`;
$text=~s/\n//g;
if($text=~/BioSample\".+?\"(SAM.+?)\"/){
	$biosample=$1;
	print "$ARGV[0]\t[$biosample]\n";
}if($text=~/BioProject\".+?\"(PRJ.+?)\"/){
	$bioproject=$1;
	print "$ARGV[0]\t[$bioproject]\n";
}else{
	die "Cannot retrieve BioProject or BioSample ID for accession $ARGV[0]";
}

$assembly=`curl ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt`;

if($biosample && $assembly=~/\n(.+$biosample.+)\n/){
	$line=$1;
	print "[$biosample][$line]\n";
}elsif($bioproject && $assembly=~/\n(.+$bioproject.+)\n/){
	$line=$1;
	print "[$bioproject][$line]\n";
}else{
	die "Cannot match BioSample/BioProject [$biosample][$bioproject] to RefSeq assembly summary\n";
}

@cols=split(/\t/,$line);

$ftpfile=$cols[19] . "/" . $cols[0] . "_" . $cols[15] . "_genomic.gbff.gz";

$command="wget --quiet $ftpfile";
`$command`;

#sleep(10);

$command="gunzip " . $cols[0] . "_" . $cols[15] . "_genomic.gbff";
`$command`;

#$command="mv " . $cols[0] . "_" . $cols[15] . "_genomic.gbff $ARGV[0].gb";
$command="gbk_split2.pl " . $cols[0] . "_" . $cols[15] . "_genomic.gbff";
@files=`$command`;

#delete gbff
$command="rm " . $cols[0] . "_" . $cols[15] . "_genomic.gbff";
system($command);


__END__
