#!/usr/bin/perl
##
## 	Patrick Degnan
##	gbk_strip2.pl
##	Convert gb/gbk/gbff files into GenBank ptt, rnt or combined annotation files and fna files
##	Does not use BioPerl
##

$/="     gene            ";
@types=("CDS","rRNA","RNA","tRNA","tmRNA","misc_RNA","ncRNA");
unless(defined($ARGV[0])){die "Usage gbk_strip2.pl <1=ptt/2=ptt+rnt/3=rnt/4=pttp> <GenBank formatted file> \n";}
$outfmt=$ARGV[0];
$file=$ARGV[1];
$seq="N";


open(IN,"<$file") or die "Cannot open $file\n";

$outtbl=$file;
@extensions=qw(null .ptt .txt .rnt _pttp.txt);
$outtbl=~s/(\..+)/$extensions[$outfmt]/;
open(OUT,">$outtbl") or die "Cannot open $outtbl\n";
$outfna=$file;
$outfna=~s/(\..+)/.fna/;
open(FA,">$outfna") or die "Cannot open $outfna\n";

$first=<IN>;
$first=~/LOCUS\s+(.+?)\s+(\d+)\sbp/;
$accession=$1;
$genome_length=$2;
$first=~s/\n//g;
$first=~/DEFINITION\s+(.+?)\.\s*ACCESSION\s+$accession/;
$organism=$1;
$organism=~ s/\s+/ /g;

$head="";$body="";$count=0;	
if($outfmt =~/1|2|3/){
	$head="$organism - 1..$genome_length\nXXXX proteins\n";
	$head=$head . "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n";
}else{$head="Start\tStop\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\tType\tNT\tOld_locus\n";}


while($line=<IN>){

	if($line=~ /\/locus\_tag\=\"(.*?)\"/){
 		$name=$1;
 		#$entire= ">gnl|NAM|$name ";
 		#print "$name\t";
 	}else{$name=""; }#print "\t";
	
	
 	if($line=~ /gene\=\"(.*)\"/){
		$gene=$1;
		#print "$gene\t";
		#$entire= $entire . "$gene ";
	}else{$gene="";}#print "\t";
	
	if($line=~ /product\=\"(.*?)\"\s+.+?/s){
		$product=$1;
		$product=~ s/\n\s+/ /g;
		#$entire= $entire . "$product ";
		#print "$product\t";
		if($line=~ /note\=\"(.*?)\"\s+.+?/s){
			$added=$1;
			$added=~ s/\n\s+/ /g;
			$product=$product . "; $added";
		}
	}elsif($line=~ /note\=\"(.*?)\"\s+.+?/s){
		$product=$1;
		$product=~ s/\n\s+/ /g;
		#print "$note\t";
	}else{$product="";}#print "\t";
	
	if($line=~/\/db\_xref\=\"GI\:(\d+)\"/){
		$gi=$1;
	}else{$gi="";}

	if($line=~/\/db\_xref\=\"COG\:(COG\d+)\"/){
		$cog=$1;
	}else{$cog="-";}	

	if($line=~/EcoGene\:(EG\d+)\"/){
		$eg=$1;
	}elsif($line=~ /\/matching\_locus\_tag\=\"(.*?)\"/){
 		$eg=$1;
 	}else{$eg=""; }#print "\t";
 	
		
	if($line=~ /\/pseudo\n/){
		if($line=~ /<?(\d+)\.\.>?(\d+)/){
			$start=$1;
			$stop=$2;
			$strand="+";
		}elsif($line=~ /complement\(<?(\d+)\.\.>?(\d+)\)/){
			$start=$1;
			$stop=$2;
			$strand="-";	
		}
		$type="pseudogene";	
	}else{		
		foreach $cat (@types){
			if($line=~ /\s+$cat\s+<?(\d+)\.\.>?(\d+)/){
				$start=$1;
				$stop=$2;
				$strand="+";
				$type=$cat;
				last;
			}elsif($line=~ /\s+$cat\s+complement\(<?(\d+)\.\.>?(\d+)\)/){
				$start=$1;
				$stop=$2;
				$strand="-";
				$type=$cat;	
				last;
			}
		}
	}
	if($type eq "CDS"){$aa=((($stop-$start+1)/3) - 1);}
	else{$aa="";}
	$nt=$stop-$start+1;
	$product=~ s/\n\s+/ /g;
	if(($outfmt == 1 || $outfmt == 2 ) && $type eq "CDS"){
		$body=$body . "$start..$stop\t$strand\t$aa\t$gi\t$gene\t$name\t-\t$cog\t$product\n";
		$count++;
	}elsif(($outfmt == 2 || $outfmt == 3 ) && $type =~ /RNA/){
		$body=$body . "$start..$stop\t$strand\t$nt\t$gi\t$gene\t$name\t-\t$cog\t$product\n";
		$count++;
	}elsif($outfmt == 4){
		$body=$body . "$start\t$stop\t$strand\t$aa\t$gi\t$gene\t$name\t-\t$cog\t$product\t$type\t$nt\t$eg\n";
		$count++;
	}
	#print "$start..$stop\t$strand\t$aa\t$gi\t$gene\t$name\t-\t-\t$product\t$type\t$eg\n";
	#print "$name\t$gene\t$eg\n";
	#if($line=~ /translation\=\"(.*?)\"\s+.+?/s){
	#	$seq=$1;
	#	$seq=~s/[\s\n]//g;
	#	print "$entire\n$seq\n";
	#}
	
	if($line=~/ORIGIN(.*)\/\//s){
		#print " found the origin!\n";
		$seq=$1;
		print FA ">$accession $organism\n";
	#}elsif($seq eq "Y" && $line =~ /\s+\d+/){
		$seq=~s/\d+//g;
		$seq=~s/\s+//g;
		$seq=~s/(.{60})/$1\n/g;
		print FA "$seq\n";	
	}
 
}
if($outfmt == 1 || $outfmt == 2){$head=~s/XXXX\sproteins/$count proteins/;}
elsif($outfmt == 3){$head=~s/XXXX\sproteins/$count RNAs/;}

#unless($outfmt == 4){}
print OUT "$head";
print OUT $body;

__END__
