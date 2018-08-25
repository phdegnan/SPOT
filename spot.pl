#!/usr/bin/perl
##
## 	Patrick Degnan
##	spot.pl
##	SPOT = sRNA-target Prediction Organizing Tool 
##
## v0.1 Patrick Degnan
## v0.2 Jan 12 2015
## v0.3 Jan 7  2016
## v0.4 Jan    2017
## v0.5 Aug    2017
## v0.6 Sept   2017
## v0.7 Sept   2017
## v1.0 Aug    2018
##

use Getopt::Long;
unless(@ARGV){
	print "\nUsage $0 
Input parameters:
\t-r\tFasta file of sRNA query
\t-a\tRefSeq Accession number (assumes any local files have RefSeq number as their prefix)
\t-o\toutput file prefix (default = TEST)
\t-g\tUse local GenBank or PTT&FNA files for all Programs? (default = N use latest from GenBank,
\t\tCopraRNA cannot use local files)
\t-n\tOther genome RefSeq ids for CopraRNA listed in quotes '' , current max is 5 genomes (default ='')
\t-m\tMultisequence sRNA file for each genome in CopraRNA list (default ='')
\t-x\tEmail address for job completion notification (default ='')
\nAlgorithm parameters:
\t-u\tNumber of nt upstream of start site to search (default = 60)
\t-d\tNumber of nt downstream of start site to search (default = 60)	
\t-s\tseed sizes for I, T, S e.g., '6 7 6' (defaults TargetRNA = 7, IntaRNA & Starpicker = 6)
\t-c\tP/Threshold value Cutoff for T, S, I e.g., '0.5 .001 un'
\t\t(defaults Target = 0.05, Starpicker = 0.5, IntaRNA = top)
\nResults Filters:
\t-b\tNumber of nt upstream of start site to filter results (default = -20)
\t-e\tNumber of nt downstream of start site to filter results (default = 20)
\t\tNote: -b and -e ignored if using a list (-l) or Rockhopper results (-t)
\t\tNote: Set -b and -e to -u and -d to get all possible matches in results
\t-l\tList of up and/or down regulated genes, include binding coord if known e.g., 
\t\tb1101\\tdown\\n
\t\tb3826\\tup\\tsRNA_start\\tsRNA_stop\\tmRNA_start\\tmRNA_stop\\n
\t\tOR
\t-t\ttranscriptome expression file from Rockhopper *_transcripts.txt
\t-f\tRockhopper fold change cutoff (default = 1.5)
\t-q\tRockhopper q value cutoff (default = 0.01)
\t-k\tRockhopper RPKM cutoff value (default = 100)
\t-p\tOperon file from DOOR-2 (http://csbl.bmb.uga.edu/DOOR/index.php) (optional)
\t-w\tReport all genes even if List or Rockhopper provided ? (default = No)
\t-y\tExclude target predictions by only 1 method ? (default = Yes)
\t\tNote: Does not apply to genes on List or significantly expressed from Rockhopper
\t-z\tSkip sRNA-mRNA detection steps, and just re-analyze data [Yy]es (default = No)
\t\t(Run in the same directory & requires original results files from each program)
\n\nExample run to examine entire genome target mRNAs:
$0 -r sgrS.fasta -o stringent -a NC_000913 \n
Example of running for select target mRNAs:
$0 -r sgrS.fasta -l sgrS_diff.txt -u 150 -d 100 -o relaxed -a NC_000913\n
Example of running for reanalysis of data with different filters:
$0 -r sgrS.fasta -l sgrS_diff.txt -a NC_000913 -o changed_50_30 -b -50 -e 30 -z Y\n
Example run using Rockhopper results and CopraRNA :
$0 -r sgrS.fasta -t NC_000913_SgrS_transcripts.txt -o express -a NC_000913 -m sgrS_homologs.fasta -n 'NC_002695 NC_008563' \n\n\n";	
	exit;
}

## unused -h -i -j -v 
GetOptions("rna=s"=>\$rna,"up=s"=>\$upstream,"down=s"=>\$downstream
	   ,"out=s"=>\$out,"seed=s"=>\$seed,"acc=s"=>\$acc
	   ,"begin=s"=>\$begin,"end=s"=>\$end,"l=s"=>\$list
	   ,"z=s"=>\$skip,"t=s"=>\$expression,"p=s"=>\$operon,"g=s"=>\$local
	   ,"q=s"=>\$qval,"k=s"=>\$rpkm,"f=s"=>\$fold, "c=s"=>\$cutoff
	   ,"n=s"=>\$copra_genomes,"m=s"=>\$rnas,"x=s"=>\$email
	   ,"y=s"=>\$num_predictions,"w=s"=>\$include_all);

#set default parameters

unless($upstream=~/\d/){$upstream=60;}# targetRNA default
if($upstream == 0){$upstream=1;print "Warning: US flank cannot equal 0 - set to 1\n";}
unless($downstream=~/\d/){$downstream=60;}
if($downstream == 0){$downstream=1;print "Warning: DS flank cannot equal 0 - set to 1\n";}
unless($begin=~/\d/){$begin=-20;}
unless($end=~/\d/){$end=20;}
unless($local){$local="N";}

unless($num_predictions){$num_predictions="Y";}
unless($include_all){$include_all="N";}

if($seed){($seedi,$seedt,$seeds)=split(" ",$seed);}
else{$seedi=6;$seedt=7;$seeds=6;}
#TargetRNA2 default=7
#IntRNA default=5
#StarPicker default=7
if($cutoff){($cutofft,$cutoffs,$cutoffi)=split(" ",$cutoff);}
else{$cutofft=0.05;$cutoffs=0.5;$cutoffi="top";}# program defaults

unless($out){$out="TEST";}
unless($skip){$skip="N";}
unless($list){$list="";}
unless($expression){$expression="";}
unless($operon){$operon="";}
unless($qval){$qval=0.01;}
unless($rpkm){$rpkm=100;}
unless($fold){$upfold=1.5;$downfold= -1.5;}#$downfold=2/3;
else{$upfold=$fold;$downfold= -1 * $fold;}#$downfold=1/$fold;
unless($email){$email="";}

#set Paths

$PWD=`pwd`; chomp($PWD);
$STARDIR="/software/starpicker";
$TARGDIR="/software/TargetRNA2";
$INTADIR="/software/IntaRNA1.0.4";
$COPRADIR="/software/CopraRNA/1.2.9";


unless($skip =~ /[Yy]/){ # skip search, go straight to analysis
print "\n=========Prepping RefSeq Files====================\n";
$d=`date`;chomp($d);
print "[$d]\n\n";

##find directory that has matching accession number
#$acc_dir="/data/DB/Genomes/Escherichia_coli_K_12_substr__MG1655_uid57779";
#TargetRNA2 will only work on refseq files in /data/DB/Genomes ?
#Make it use local copies?


## Using local files
if($local =~ /[Yy]/){

	$file1="$acc.ptt";
	$file2="$acc.fna";
	$file3="$acc.gb";
	$file4="$acc.gbk";

	if(-e $file1){# if file1 exists
		$acc_dir= $PWD;
		if(-e $file2 && (-e $file3 || -e $file4)  ){
			# do nothing
		}elsif(-e $file2 && !(-e $file3) && !(-e $file4)){
			system("/scripts/fnaptt2gbk.pl -f $file2 -t $file1");
		}else{# 
			die "*PTT file present locally, but *FNA file [$file2] and/or but *GBK file [$file3] are not\n";
		}
			
	}elsif(-e $file2 && !(-e $file1) && !(-e $file3) && !(-e $file4) ){
		die "*FNA file present locally, but *PTT file [$file1] and/or but *GBK file [$file3] are not\n";
	
	}elsif(-e $file3  && !(-e $file1) && !(-e $file2) ){
		$acc_dir= $PWD;
		system("/scripts/gbk_strip2.pl 1 $file3");
	
	}elsif(-e $file4 && !(-e $file1) && !(-e $file2) ){
		$acc_dir= $PWD;
		system("/scripts/gbk_strip2.pl 1 $file4");
	
	}else{ ## Use cached genome copy (old locus numbers?)
		$total_dir=`ls /data/DB/Genomes/*/$acc.ptt`;
		if($total_dir=~/(.+?)\/$acc/){$acc_dir=$1;}
		else{ die "Cannot find $acc files in current directory or in /data/DB/Genomes/\n";}
		#Copy *fna and *ptt file to local directory, possibly delete later - starpicker cannot write to that directory.
		`cp $acc_dir/$acc.fna ./`;
		`cp $acc_dir/$acc.ptt ./`;
		system("/scripts/fnaptt2gbk.pl -f $file2 -t $file1");
		
	}
}else{
	$acc_dir= $PWD;
	#Makes sure that all three programs are using the same data - ignores old refseq that TargetRNA downloaded.
	system("/software/IntaRNA1.0.4/get_refseq_from_ftp.pl $acc");
	system("/scripts/gbk_strip2.pl 1 $acc.gb");
	
}
#print "$acc_dir\n";
#} ## debug of file input test
#__END__
print "\n=========Starting Searches========================\n";
$d=`date`;chomp($d);
print "[$d]\n\n";

## Run TargetRNA2
#=begin
chdir($TARGDIR);
#$path=`pwd`;
#print $path;
$pid="";
defined($pid = fork()) or die "unable to fork: $!\n";
if ($pid == 0) { # child
	#print "[$pid]\n";
	# Change P value threshold? P= 0.5? 0.99 ? (TargetRNA2 default = 0.05)
	$command="java -jar TargetRNA2.jar -a $upstream -b $downstream -l $seedt -y -g $acc_dir/$acc -s $PWD/$rna -d Bacteria -o $PWD/TargetRNA2_$out.txt -p $cutofft";
	print "$command\n\n";
	exec("java","-jar","TargetRNA2.jar","-a","$upstream","-b","$downstream","-l","$seedt","-y","-g","$acc_dir/$acc","-s","$PWD/$rna","-d","Bacteria","-o","$PWD/TargetRNA2_$out.txt","-p","$cutofft");
	die "unable to exec: $!\n";
}else{
	push(@IDS,$pid);
}
chdir($PWD);
#$path=`pwd`;
#print $path

## Run sTarPicker_global2.pl

$pid="";
defined($pid = fork()) or die "unable to fork: $!\n";
if ($pid == 0) { # child
	#print "[$pid]\n";
	# Change P value threshold? -th 0.001 ? 0.01 ? (sTarPicker default = 0.5)
	$command="$STARDIR/sTarPicker_global2.pl  -srna $rna -gen $acc.fna -ptt $acc.ptt -th $cutoffs -u $upstream -d $downstream -seed $seeds &";
	print "$command\n\n";
	exec("$STARDIR/sTarPicker_global2.pl","-srna","$rna","-gen","$acc.fna","-ptt","$acc.ptt","-th","$cutoffs","-u","$upstream","-d","$downstream","-seed","$seeds","&");
	die "unable to exec: $!\n";
}else{
	push(@IDS,$pid);
}
sleep(2);

## Run IntaRNA_wrapper.pl

$pid="";
defined($pid = fork()) or die "unable to fork: $!\n";
if ($pid == 0) { # child
	#print "[$pid]\n";
	# -genbank $local changed to -genbank Y because first step above is to download if not local. No need to overwrite/repeat step
	$command="$INTADIR/IntaRNA_wrapper.pl -bindseq $rna -ncbiacc $acc -ntupstream $upstream -ntdownstream $downstream -region 5utr -bpinseed $seedi -genbank Y &";
	print "$command\n\n";
	exec("$INTADIR/IntaRNA_wrapper.pl","-bindseq","$rna","-ncbiacc","$acc","-ntupstream","$upstream","-ntdownstream","$downstream","-region","5utr","-bpinseed","$seedi","-genbank","Y","&");
	die "unable to exec: $!\n";
}else{
	push(@IDS,$pid);
}
sleep(2);
#=begin
#=end
#=cut
## Run CopraRNA.pl
#print "$copra_genomes $rnas\n";
if($copra_genomes =~/.+/ && -e $rnas ){
	unless(-d COPRASUB){`mkdir COPRASUB`;}
	`cp $rnas COPRASUB`;
	chdir("$PWD/COPRASUB");
	$pid=""; 
	#print "pid = [$pid]\n";
	defined($pid = fork()) or die "unable to fork: $!\n";
	#print "pid = [$pid]\n";
	if ($pid == 0) { # child
		#print "[$pid]\n";
		$command="$COPRADIR/homology_intaRNA.pl $rnas $upstream $downstream 5utr $acc $copra_genomes &";
		print "$command\n\n";
		@cgenomes=split(/\s/,$copra_genomes);
		#print "\$copra_genomes = '$copra_genomes'\n";
		#print scalar(@cgenomes) ."[$cgenomes[0]][$cgenomes[1]]\n";
		if(@cgenomes == 2){
			exec("$COPRADIR/homology_intaRNA.pl","$rnas","$upstream","$downstream","5utr","$acc","$cgenomes[0]","$cgenomes[1]");
			die "unable to exec: $!\n";
		}elsif(@cgenomes == 3){
			exec("$COPRADIR/homology_intaRNA.pl","$rnas","$upstream","$downstream","5utr","$acc","$cgenomes[0]","$cgenomes[1]","$cgenomes[2]");
			die "unable to exec: $!\n";
		}elsif(@cgenomes == 4){
			exec("$COPRADIR/homology_intaRNA.pl","$rnas","$upstream","$downstream","5utr","$acc","$cgenomes[0]","$cgenomes[1]","$cgenomes[2]","$cgenomes[3]");
			die "unable to exec: $!\n";
		}elsif(@cgenomes == 5){
			exec("$COPRADIR/homology_intaRNA.pl","$rnas","$upstream","$downstream","5utr","$acc","$cgenomes[0]","$cgenomes[1]","$cgenomes[2]","$cgenomes[3]","$cgenomes[4]");
			die "unable to exec: $!\n";
		}
	}else{
		push(@IDS,$pid);
	}
	#chdir($PWD);
	sleep(2);
}elsif($rnas){
	print "**Missing genome RefSeq numbers to analyze**\n\tCopraRNA not Run!\n";
}elsif($copra_genomes){
	print "**Missing multisequence sRNA file to analyze**\n\tCopraRNA not Run!\n";
}# if neither is true, simply skip to next task
#=end
#=cut
#wait until jobs done

foreach $p (@IDS){
	$d=`date`;chomp($d);
	print "RUNNING[$p][$d]\n";
	waitpid($p,0); ## found code online - debugged in test_fork.pl and worked
	$d=`date`;chomp($d);
	print "COMPLETED[$p][$d]\n";
}

chdir($PWD); #return to pwd if left for CopraRNA

##Check for problem with first run by IntaRNA - should be fixed GeneID vs. GI issue...
$check=`grep 'Annotation' termClusterReport.txt`;
unless($check =~ /Annotation/){	
	print "\n$INTADIR/rerun_enrichment.pl intarna_websrv_table.csv\n\n";
	system("$INTADIR/rerun_enrichment.pl intarna_websrv_table.csv");
}

print "\n=========Jobs Complete============================\n";

}## end of skip
#__END__

## Process/Analyze results
## Process IntaRNA results and SAVE IN HASH
print "\n=========Reading/Saving Raw Results===============\n";
$d=`date`;chomp($d);
print "[$d]\n";
$rna_name=`grep ">" $rna`;
chomp($rna_name);
$rna_name=~s/\>//;

$underscore="N";
if($cutoffi=~/top/i){
	$inta_results="intarna_websrv_table_truncated.csv";
}else{
	$inta_results="intarna_websrv_table.csv";
}
open(IN,"<$inta_results") or die "cannot open $inta_results\n";
$ranked="N";$structure="N";
while($l=<IN>){
	chomp($l);
	@cols=split(/;/,$l);
	if($cutoffi == "top" || $cutoffi == "un" || ($cutoffi=~/\d/ && $cols[0] <= $cutoffi)){
		if($cols[20] eq ""){
			$gene=$cols[2];
		}else{
			$gene=$cols[20];
		}
		$INVERSE{$gene}=$cols[2];
		
		if($cols[2] =~ /\_/){$underscore="Y";}
		
		$TOTALLOCUS{$cols[2]}=$gene;
		$INTA{$cols[2]}{GENE}=$gene; # if no gene value is locus
		$INTA{$cols[2]}{EVAL}=$cols[18];
		$INTA{$cols[2]}{PVAL}=$cols[0];
		$INTA{$cols[2]}{SSTA}=$cols[12];
		$INTA{$cols[2]}{SEND}=$cols[13];
		$INTA{$cols[2]}{MSTA}=$cols[5] - $upstream - 1;
		$INTA{$cols[2]}{MEND}=$cols[6] - $upstream - 1;
		
		$GENEID{$cols[2]}=$cols[21];
	
		@lines=split(/\\n/,$cols[19]);
		foreach $r (@lines){
			$INTA{$cols[2]}{STRUCTURE}.="$r\r";
			#\r is carriage return for Macs. Excel will use and wrap text in cell.
			#Not helpful in text editor
		}
	
	}
}close(IN);
## Debug:
#print "$INTA{b1101}{GENE}\t$INTA{b1101}{EVAL}\t$INTA{b1101}{PVAL}\t$INTA{b1101}{SSTA}\t$INTA{b1101}{SEND}\t$INTA{b1101}{MSTA}\t$INTA{b1101}{MEND}\n";
#print  $INTA{b1101}{STRUCTURE};

## Process sTarpicker_global.pl results and SAVE IN HASH

$star_results=`ls *.output`;
chomp($star_results);
#print "$star_results\n";
open(IN,"<$star_results") or die "cannot open $star_results\n";
while($l=<IN>){
	chomp($l);
	
	if($l=~/^\d/){
		@cols=split(/\t/,$l);
		#handle exceptions where locus numbers have underscores in their names
		if($underscore eq "Y"){
			@names=split(/\_/,$cols[1]);
			if(@names == 3){
				$gene=$names[0];
				$locus=$names[1] ."_". $names[2];	
			}else{
				$gene=$names[0] ."_". $names[1];
				$locus=$names[0] ."_". $names[1];
			}
			#DEBUG
			#if($cols[1] eq "BT_4680"){print "[$cols[1]]\t[$names[0]]\t[$names[1]]\n[$gene]\t[$locus]\n";}
	
		}elsif($cols[1]=~/(\S+?)\_(\S+)/){
			#print "$cols[1]\n";
			$gene=$1;
			$locus=$2;
		}else{
			$gene=$cols[1];
			$locus=$cols[1];
		}
		$INVERSE{$gene}=$locus;
		#print "[$gene]\t[$locus]\n";
		$TOTALLOCUS{$locus}=$gene;
		$STAR{$locus}{GENE}=$gene;
		$first="N"
	}else{
		@cols=split(/\s+/,$l);
		unless($first eq "Y"){
			if($l=~/Probability/){
				$STAR{$locus}{EVAL}=$cols[4];
				$STAR{$locus}{PVAL}=1-$cols[3];
			}elsif($l=~/sRNA/){
				$STAR{$locus}{SSTA}=$cols[2];
				$STAR{$locus}{SEND}=$cols[4];
				$l=~s/^\t//;
				$STAR{$locus}{STRUCTURE}="";
				$STAR{$locus}{STRUCTURE}.=$l . "\r"; 
			}elsif($l=~/\|/){
				$l=~s/^\t//;
				$STAR{$locus}{STRUCTURE}.=$l . "\r";
			}elsif($l=~/Target/){
				$STAR{$locus}{MSTA}=$cols[4] ;
				$STAR{$locus}{MEND}=$cols[2] ;
				$l=~s/^\t//;
				$STAR{$locus}{STRUCTURE}.=$l . "\r";
				$first="Y";
			}
		}
	}
	
}close(IN);
##Debug
#print "$STAR{b1101}{GENE}\t$STAR{b1101}{EVAL}\t$STAR{b1101}{PVAL}\t$STAR{b1101}{SSTA}\t$STAR{b1101}{SEND}\t$STAR{b1101}{MSTA}\t$STAR{b1101}{MEND}\n";
#print  $STAR{b1101}{STRUCTURE};

## Process TargetRNA2 results and SAVE IN HASH
$targ_results=`ls TargetRNA2_*.txt`;
chomp($targ_results);
open(IN,"<$targ_results") or die "cannot open $targ_results\n";
$ranked="N";$structure="N";
while($l=<IN>){
	chomp($l);
	if($l=~/^Rank/){
		$ranked="Y";
	}elsif($l=~/^\%/){
		$ranked="N";
		$structure="Y";
		$NEXT=<IN>;
		$NEXT=~/^(\S+)/;
		$locus=$INVERSE{$1};
		#print "$1\t$GENE\n";
	}elsif($ranked eq "Y"){
		@cols=split(/\s+/,$l);
		if($cols[1] eq "" || $cols[1] eq "-"){
			$gene=$cols[2];
		}else{
			$gene=$cols[1];
		}
		$TOTALLOCUS{$cols[2]}=$gene;
		$TARGET{$cols[2]}{GENE}=$gene;
		$TARGET{$cols[2]}{EVAL}=$cols[3];
		$TARGET{$cols[2]}{PVAL}=$cols[4];
		$TARGET{$cols[2]}{SSTA}=$cols[5];
		$TARGET{$cols[2]}{SEND}=$cols[6];
		$TARGET{$cols[2]}{MSTA}=$cols[7];
		$TARGET{$cols[2]}{MEND}=$cols[8];
		$INVERSE{$gene}=$cols[2];
		
	}elsif($structure eq "Y"){
		if($l=~/^\s+.+/){
			$l=~s/^\t//;
			## Parse output to remove \t which have issues when viewed in Excel
			($gene,$coord1,$seq,$coord2)=split(/\t/,$l);
			#12 char, 4 char, whatever, whatever
			$lg=12-length($gene);
			$string=$gene;
			for $i (1..$lg){$string.=" ";}
			$lg=4-length($coord1);
			$string.=$coord1;
			for $i (1..$lg){$string.=" ";}
			$string.=" $seq   $coord2";
			$TARGET{$locus}{STRUCTURE}.="$string\r";
			#print "$string\n";
		}
	}
	
}close(IN);
## Debug
#print "$TARGET{b1101}{GENE}\t$TARGET{b1101}{EVAL}\t$TARGET{b1101}{PVAL}\t$TARGET{b1101}{SSTA}\t$TARGET{b1101}{SEND}\t$TARGET{b1101}{MSTA}\t$TARGET{b1101}{MEND}\n";
#print  $TARGET{b1101}{STRUCTURE};

if($operon){
	open(IN,"<$operon") or die "cannot open $operon\n";
	while($l=<IN>){
		chomp($l);
		@cols=split(/\t/,$l);
		$UNITS{$cols[2]}=$cols[0];
		$OPERONS{$cols[0]}.="$cols[2],";
	}close(IN);
}else{
	%OPERONS=();
	%UNITS=();
}


## Read in of Up or Down regulated genes
$yes="N";
%ALL=();
if($list){
	print "\n=========Sig. different expression================\n";
	$d=`date`;chomp($d);
	print "[$d]\n";
	open(IN,"<$list") or die "cannot open $list\n";
	while($l=<IN>){
		chomp($l);
		@cols=split(/\t/,$l);
		$EXPRESS{$cols[0]}="$cols[1]\t\t";
		print "$cols[0]\t$cols[1]\t";
		#if known coordinate data are provided in LIST
		if($cols[2] =~ /\d/ && $cols[3] =~ /\d/ && $cols[4] =~ /\d/ && $cols[5] =~ /\d/){
			$KNOWN{$cols[0]}{FLAG}=1;
			$KNOWN{$cols[0]}{SSTA}=$cols[2];
			$KNOWN{$cols[0]}{SEND}=$cols[3];
			$KNOWN{$cols[0]}{MSTA}=$cols[4];
			$KNOWN{$cols[0]}{MEND}=$cols[5];
			print "$cols[2]\t$cols[3]\t$cols[4]\t$cols[5]\t";
		}
		print "\n";
		
	}close(IN);
	$EXPRESS{BLANK}="UP/DOWN\tQ_value\t";
	$yes="Y";
}
if($expression){## read expression output from Rockhopper *_transcripts.txt file
	open(IN,"<$expression") or die "cannot open $expression\n";
	print "\n=========Sig. different expression================\n";
	$d=`date`;chomp($d);
	print "[$d]\n";
	while($l=<IN>){
		chomp($l);
		@cols=split(/\t/,$l);
		if($cols[8] == 0){$cols[8] =1;}
		if($cols[9] == 0){$cols[9] =1;}
		
		$change=$cols[8] / $cols[9]; 
		if($change < 1){
			$invert= -1*(1/$change);
			$change=$invert;
		}
		$ALL{$cols[6]}="$change\t$cols[10]\t";
		
		if( ($cols[8]+$cols[9]) >= $rpkm && $cols[10] <= $qval && $change >= $upfold){
			$EXPRESS{$cols[6]}="$change\t$cols[10]\t";
			print "$cols[6]\tup\n";
		}elsif( ($cols[8]+$cols[9]) >= $rpkm && $cols[10] <= $qval && $change <= $downfold ){
			$EXPRESS{$cols[6]}="$change\t$cols[10]\t";
			print "$cols[6]\tdown\n";
		}

	}close(IN);
	$EXPRESS{BLANK}="UP/DOWN\tQ_value\t";
	$yes="Y";
}
if($yes eq "N"){# no list provided
	$EXPRESS{BLANK}="";
}

##EXPAND EXPRESSION RESULTS TO GENES IN SAME OPERONS
if($yes eq "Y" && -e $operon){## if $list OR $expression was true AND $operon
	print "\n=========Added genes in operons===================\n";
	$d=`date`;chomp($d);
	print "[$d]\n";
	foreach $locus (keys(%EXPRESS)){#go through sig. results
		if($UNITS{$locus}=~/(\d+)/){# get operon number if true
			$op=$1;
			$OPERONS{$op}=~s/\,$//; #remove final comma
			@inoperon=split(/\,/,$OPERONS{$op});#split on commas
			foreach $i (@inoperon){ # for each locus give that new locus the expression up/down of original locus
				unless($EXPRESS{$i}=~/\w/){
					if($ALL{$i}=~/\w/){
						$added_value=$ALL{$i}; # if Rockhopper
					}else{
						$added_value="n.a.";   # if just a list and using operon
					}
					$EXPRESS{$i}=$added_value;
					print "$i\t$added_value\n";
				}
			}
		}
	}
}

#=begin

## Extract Copra Results... if it was run - Columns in final table will be bl
%COPRA=();
if($copra_genomes && $rnas ){
	#Add structures to the predictions of the target genome ($acc)
	chdir("$PWD/COPRASUB");
	unless(-e "intarna_websrv_table_ncbi.csv"){
		$initial_int=`ls $acc*intarna`;
		chomp($initial_int);
		system "cp $initial_int predictions.intarna";
		# Code from IntaRNA_wrapper.pl
		# process intarna output
		system "$INTADIR/produce_semicolon_sep_results_from_intarna_out.pl --intarna-out-file predictions.intarna > predictions.csv";
		# sort intarna output
		system "$INTADIR/sort_intarna_csv_results.pl --intarna-csv-file predictions.csv --column 17 > predictions.sorted.csv";
		# omit WARNING for 'N' from IntaRNA raw output. fixes interaction parsing warnings.
		system "sed -i '/WARNING/{N;d;}' predictions.intarna"; ## d stands for delete N also deletes one line after. add another N to delete one more ##edit 1.0.1 
		# add interactions to table
		system "$INTADIR/add_interactions.pl predictions.sorted.csv predictions.intarna"; #writes to: intarna_websrv_table.csv
		# if ncbi option is on then also add annatations and Entraz Gene IDs and gene name 
		system "$INTADIR/add_GI_genename_annotation.pl"; ## writes to: intarna_websrv_table_ncbi.csv
	}
	chdir("$PWD");
	
	#Data in $PWD/COPRASUB
	$copra_results="$PWD/COPRASUB/CopraRNA_result.csv";
## Only has 100 top results / not a P value cutoff ?	
## data format:
#0	1	2		3		4		5		6		7	
#fdr	p-value	NC_000913	NC_009792	NC_013716	NC_011740	Annotation	Additional homologs	
#Each entry
#	locus(gene|c-energy|c-p-value|target_start|target_end|sRNA_start|sRNA_end|GeneID:XXXXX)	
	
	open(IN,"<$copra_results") or die "cannot open $copra_results\n";
	$head=<IN>;
	@cols=split(/\,/,$head);
	$last_column=@cols - 3;
	$anno=@cols - 2;
	while($l=<IN>){
		chomp($l);
		@cols=split(/\,/,$l);
		$cols[2]=~s/\(/\|/;
		@parts=split(/\|/,$cols[2]);
		$locus=$parts[0];
		if($parts[1] eq "N/A"){
			$COPRA{$locus}{GENE}=$locus; # if no gene value is locus
		}else{
			$COPRA{$locus}{GENE}=$parts[1];
		}
		$COPRA{$locus}{EVAL}=$parts[2];
		$COPRA{$locus}{PVAL}=$parts[3];	#Individual IntaRNA p Value
		$COPRA{$locus}{SSTA}=$parts[6];
		$COPRA{$locus}{SEND}=$parts[7];
		$COPRA{$locus}{MSTA}=$parts[4] - $upstream - 1;
		$COPRA{$locus}{MEND}=$parts[5] - $upstream - 1;
		$COPRA{$locus}{CPVAL}=$cols[1];	# CopraRNA composite P value
		$COPRA{$locus}{ANNO}=$cols[$anno];	# annotation
		$COPRA{$locus}{ANNO}=~s/\;//g;
		#print "[3][$last_column]\n";
		for $i (3..$last_column){ # not all genomes will have a match...
			#print "[$cols[$i]]\n";
			if($cols[$i] =~/(.+?)\(/){
				$COPRA{$locus}{HOMOLOGS}.="$1,";
			}
		}
		$COPRA{$locus}{HOMOLOGS}=~s/\,$//;# delete final comma

	}close(IN);
	## get structures
	$copra_acc_results="$PWD/COPRASUB/intarna_websrv_table_ncbi.csv";
	open(IN,"<$copra_acc_results") or die "cannot open $copra_acc_results\n";
	while($l=<IN>){
		chomp($l);
		@cols=split(/\;/,$l);
		$locus=$cols[0];
		#print "[$cols[0]]\t[$cols[17]]\t[$cols[19]]\n";
		@lines=split(/\\n/,$cols[17]);
		foreach $r (@lines){
			$COPRA{$locus}{STRUCTURE}.="$r\r";
			#\r is carriage return for Macs. Excel will use and wrap text in cell.
			#Not helpful in text editor
		}
		#print "[ $COPRA{$locus}{STRUCTURE} ]\n";
	}
}
#=end
#=cut

## How to present results of up and down? Simply add column? Generate seperate out files? Only print those with matches?

## Condense results to summary of all matches AND complete data table for interactions -20/+20
## No P value or E value cutoffs employed here, could be added.
print "\n=========Collating Results========================\n";
$d=`date`;chomp($d);
print "[$d]\n";
$sub="COLLATED_RESULTS_" . $out;
unless(-d $sub){mkdir($sub);}

## ADD CODE TO NOT PRINT CopraRNA columns or headers when there are no results for it...
## Currently empty columns will be printed.
if($copra_genomes && $rnas ){
	$coprahead1="\tCopraRNA";
	$coprahead2="\tCopraRNA_$begin\_$end";
	$coprahead3="\tCopraRNA_M";
}

# Open out files
open(SUM,">$out\_summary.txt") or die "Cannot open $out\_summary.txt\n";
# Print headers
print SUM "Locus\tGene\t$EXPRESS{BLANK}TargetRNA2\tStar\tIntaRNA$coprahead1\tCount_ALL\tTargetRNA2_$begin\_$end\tStar_$begin\_$end\tIntaRNA_$begin\_$end$coprahead2\tCount_$begin\_$end";
print SUM "\tTargetRNA2_M\tStar_M\tIntaRNA_M$coprahead3\tRank\n";

#open(SOM,">$out\_complete.txt") or die "Cannot open $out\_complete.txt\n";
open(COM,">$out\_complete.txt") or die "Cannot open $out\_complete.txt\n";
# Print headers
print COM "Locus\tGene\t$EXPRESS{BLANK}T-Energy\tT-Pvalue\tT-sRNA_start\tT-sRNA_stop\tT-mRNA_start\tT-mRNA_stop\tT-Structure\t";
print COM "S-Energy\tS-Pvalue\tS-sRNA_start\tS-sRNA_stop\tS-mRNA_start\tS-mRNA_stop\tS-Structure\t";
print COM "I-Energy\tI-Pvalue\tI-sRNA_start\tI-sRNA_stop\tI-mRNA_start\tI-mRNA_stop\tI-Structure";
if($copra_genomes && $rnas ){
	print COM "\tC-Pvalue\tC-Energy\tC-sRNA_start\tC-sRNA_stop\tC-mRNA_start\tC-mRNA_stop\tC-Structure\tC-IPvalue\tC-Homologs\tC-Annotation";
}
print COM "\n";

open(CSV,">$out\_intarna_websrv_tbl.csv") or die "Cannot open $out\_intarna_websrv_tbl.csv\n";

open(FUN,">$sub/$out\_composite_websrv_tbl.csv") or die "Cannot open $out\_composite_websrv_tbl.csv\n";
# Print headers
print FUN "p-value,Method1,Method2,Method3,Annotation,Additional homologs\n";
open(FDR,">$sub/CopraRNA_result_all.csv") or die "Cannot open CopraRNA_result_all.csv\n";
# Print headers
print FDR "fdr,p-value,Method1,Method2,Method3,Annotation,Additional homologs\n";

@KEYS=("EVAL","PVAL","SSTA","SEND","MSTA","MEND");


foreach $locus (sort {lc($a) cmp lc($b)}  keys %TOTALLOCUS) { 
	$com=""; 
	if($locus=~/\w/){ 
	@MBRNA=();
	@MERNA=();
	@SBRNA=();
	@SERNA=();
	$icount=0;$scount=0;$ccount=0;$tcount=0;
	$csv_line="";
	$sum_line="";
	$rank="";
	## if $locus in list of differentially expressed targets

	    if(	($EXPRESS{$locus}=~/\w/ && defined($list)) || ($EXPRESS{$locus}=~/\w/ && defined($expression)) 
	         || ($list eq "" && $expression eq "") || $include_all =~ /[Yy]/ ){

		$sum_line="$locus\t$TOTALLOCUS{$locus}\t";
		if($EXPRESS{$locus}!~/\t/ && $include_all =~ /[Yy]/ && (defined($list) || defined($expression))){
			$sum_line.="\t\t";
		}else{
			$sum_line.=$EXPRESS{$locus};
		}
		$count1=0;
		$count2=0;
		$region=$upstream+$downstream;
		$largewindow="";
		$smallwindow="";
		$overlapwindow="";
		if($TARGET{$locus}{GENE}){
	
			#print SUM "1\t";
			$count1++;
			$largewindow="1\t";
			if(($TARGET{$locus}{MSTA} <= $end && $TARGET{$locus}{MSTA} >= $begin) ||
			   ($TARGET{$locus}{MEND} <= $end && $TARGET{$locus}{MEND} >= $begin) ||
			   ($EXPRESS{$locus}=~/\w/ && (defined($list) || defined($expression)) ) ){
				#print SUM "1\t";
				$count2++;
				$smallwindow.="1\t";
				foreach $col (@KEYS){
					$com.="$TARGET{$locus}{$col}\t";
				}
				$TARGET{$locus}{STRUCTURE}=~s/\r$//;
				$com.="\"$TARGET{$locus}{STRUCTURE}\"\t";
				$mlen=$TARGET{$locus}{MEND} - $TARGET{$locus}{MSTA} + 1;
				$slen=$TARGET{$locus}{SEND} - $TARGET{$locus}{SSTA} + 1;
				$msta=$TARGET{$locus}{MSTA} + $upstream + 1;
				$mend=$TARGET{$locus}{MEND} + $upstream + 1;
				$stru=$TARGET{$locus}{STRUCTURE};
				$stru=~s/\r/\\n/g;
				$stru=~s/\t/\\n/g;
				$stru=~s/\;/ /g;
				$hybrid=$TARGET{$locus}{EVAL}*2;
				unless($TARGET{$locus}{GENE}=~/\w/){$gene=$locus;}
				else{$gene=$TARGET{$locus}{GENE};}
				$csv_line.="$TARGET{$locus}{PVAL};$TARGET{$locus}{PVAL};$gene\_t;$region;";
				$csv_line.="$mlen;$msta;$mend;$msta;$mend;na;$rna_name;$slen;$TARGET{$locus}{SSTA};";
				$csv_line.="$TARGET{$locus}{SEND};$TARGET{$locus}{SSTA};$TARGET{$locus}{SEND};na;";
				$csv_line.="$hybrid;$TARGET{$locus}{EVAL};$stru;$locus;na;na\n";
				
				$smallrna="";$rnatarget="";
				if($KNOWN{$locus}{FLAG} == 1){
					$smallrna=&MATCH($TARGET{$locus}{SSTA},$TARGET{$locus}{SEND},$KNOWN{$locus}{SSTA},$KNOWN{$locus}{SEND});
					$rnatarget=&MATCH($TARGET{$locus}{MSTA},$TARGET{$locus}{MEND},$KNOWN{$locus}{MSTA},$KNOWN{$locus}{MEND});
					if($smallrna > 0 && $rnatarget > 0){$overlapwindow.="A\t";}#print "A\n";}
					else{$overlapwindow.="B\t";}#print "B\n";}
				}else{$overlapwindow.="F\t";}#print "C\n";}
				push(@MBRNA,$TARGET{$locus}{MSTA});
				push(@MERNA,$TARGET{$locus}{MEND});
				push(@SBRNA,$TARGET{$locus}{SSTA});
				push(@SERNA,$TARGET{$locus}{SEND});
				$tcount=1;
			}else{
				#print SUM "\t";
				$smallwindow.="\t";
				$overlapwindow.="\t";
				$com.="\t\t\t\t\t\t\t";
				push(@MBRNA,"");push(@MERNA,"");push(@SBRNA,"");push(@SERNA,"");
			}   
			
		}else{
			#print SUM "\t\t";
			$largewindow.="\t";
			$smallwindow.="\t";
			$overlapwindow.="\t";
			$com.="\t\t\t\t\t\t\t";
			push(@MBRNA,"");push(@MERNA,"");push(@SBRNA,"");push(@SERNA,"");
		}
		
		if($STAR{$locus}{GENE}){
			#print SUM "1\t";
			$count1++;
			$largewindow.="1\t";
			if(($STAR{$locus}{MSTA} <= $end && $STAR{$locus}{MSTA} >= $begin) ||
			   ($STAR{$locus}{MEND} <= $end && $STAR{$locus}{MEND} >= $begin) ||
			   ($EXPRESS{$locus}=~/\w/ && (defined($list) || defined($expression))) ){
				#print SUM "1\t";
				$count2++;
				$smallwindow.="1\t";
				foreach $col (@KEYS){
					$com.="$STAR{$locus}{$col}\t";
				}
				$STAR{$locus}{STRUCTURE}=~s/\r$//;
				$com.="\"$STAR{$locus}{STRUCTURE}\"\t";
				$mlen=$STAR{$locus}{MEND} - $STAR{$locus}{MSTA} + 1;
				$slen=$STAR{$locus}{SEND} - $STAR{$locus}{SSTA} + 1;
				$msta=$STAR{$locus}{MSTA} + $upstream + 1;
				$mend=$STAR{$locus}{MEND} + $upstream + 1;
				$stru=$STAR{$locus}{STRUCTURE};
				$hybrid=$STAR{$locus}{EVAL}*2;
				unless($STAR{$locus}{GENE}=~/\w/){$gene=$locus;}
				else{$gene=$STAR{$locus}{GENE};}
				$stru=~s/\r/\\n/g;
				$csv_line.="$STAR{$locus}{PVAL};$STAR{$locus}{PVAL};$gene\_s;$region;$mlen;$msta;";
				$csv_line.="$mend;$msta;$mend;na;$rna_name;$slen;$STAR{$locus}{SSTA};$STAR{$locus}{SEND};";
				$csv_line.="$STAR{$locus}{SSTA};$STAR{$locus}{SEND};na;$hybrid;$STAR{$locus}{EVAL};$stru;$locus;na;na\n";
				
				$smallrna="";$rnatarget="";
				if($KNOWN{$locus}{FLAG} == 1){
					$smallrna=&MATCH($STAR{$locus}{SSTA},$STAR{$locus}{SEND},$KNOWN{$locus}{SSTA},$KNOWN{$locus}{SEND});
					$rnatarget=&MATCH($STAR{$locus}{MSTA},$STAR{$locus}{MEND},$KNOWN{$locus}{MSTA},$KNOWN{$locus}{MEND});
					if($smallrna > 0 && $rnatarget > 0){$overlapwindow.="A\t";}#print "A\n";}
					else{$overlapwindow.="B\t";}#print "B\n";}
				}else{$overlapwindow.="F\t";}#print "C\n";}
				push(@MBRNA,$STAR{$locus}{MSTA});
				push(@MERNA,$STAR{$locus}{MEND});
				push(@SBRNA,$STAR{$locus}{SSTA});
				push(@SERNA,$STAR{$locus}{SEND});
				$scount=1;
			}else{
				#print SUM "\t";
				$smallwindow.="\t";
				$overlapwindow.="\t";
				$com.="\t\t\t\t\t\t\t";
				push(@MBRNA,"");push(@MERNA,"");push(@SBRNA,"");push(@SERNA,"");
			}
		}else{
			#print SUM "\t\t";
			$largewindow.="\t";
			$smallwindow.="\t";
			$overlapwindow.="\t";
			$com.="\t\t\t\t\t\t\t";
			push(@MBRNA,"");push(@MERNA,"");push(@SBRNA,"");push(@SERNA,"");			
		}
				
		if($INTA{$locus}{GENE}){
			#print SUM "1\t";
			$count1++;
			$largewindow.="1\t";
			if(($INTA{$locus}{MSTA} <= $end && $INTA{$locus}{MSTA} >= $begin) ||
			   ($INTA{$locus}{MEND} <= $end && $INTA{$locus}{MEND} >= $begin) ||
			   ($EXPRESS{$locus}=~/\w/ && (defined($list) || defined($expression))) ){
				#print SUM "1\t";
				$count2++;
				$smallwindow.="1\t";
				foreach $col (@KEYS){
					$com.="$INTA{$locus}{$col}\t";
				}
				$INTA{$locus}{STRUCTURE}=~s/\r$//;
				$com.="\"$INTA{$locus}{STRUCTURE}\"\t";
				$mlen=$INTA{$locus}{MEND} - $INTA{$locus}{MSTA} + 1;
				$slen=$INTA{$locus}{SEND} - $INTA{$locus}{SSTA} + 1;
				$msta=$INTA{$locus}{MSTA} + $upstream + 1;
				$mend=$INTA{$locus}{MEND} + $upstream + 1;
				$stru=$INTA{$locus}{STRUCTURE};
				$hybrid=$INTA{$locus}{EVAL}*2;
				unless($INTA{$locus}{GENE}=~/\w/){$gene=$locus;}
				else{$gene=$INTA{$locus}{GENE};}
				$stru=~s/\r/\\n/g;
				$csv_line.="$INTA{$locus}{PVAL};$INTA{$locus}{PVAL};$gene\_i;$region;$mlen;$msta;";
				$csv_line.="$mend;$msta;$mend;na;$rna_name;$slen;$INTA{$locus}{SSTA};$INTA{$locus}{SEND};";
				$csv_line.="$INTA{$locus}{SSTA};$INTA{$locus}{SEND};na;$hybrid;$INTA{$locus}{EVAL};$stru;$locus;na;na\n";
				
				$smallrna="";$rnatarget="";
				if($KNOWN{$locus}{FLAG} == 1){
					$smallrna=&MATCH($INTA{$locus}{SSTA},$INTA{$locus}{SEND},$KNOWN{$locus}{SSTA},$KNOWN{$locus}{SEND});
					$rnatarget=&MATCH($INTA{$locus}{MSTA},$INTA{$locus}{MEND},$KNOWN{$locus}{MSTA},$KNOWN{$locus}{MEND});
					if($smallrna > 0 && $rnatarget > 0){$overlapwindow.="A\t";}#print "A\n";}
					else{$overlapwindow.="B\t";}#print "B\n";}
				}else{$overlapwindow.="F\t";}#print "C\n";}
				push(@MBRNA,$INTA{$locus}{MSTA});
				push(@MERNA,$INTA{$locus}{MEND});
				push(@SBRNA,$INTA{$locus}{SSTA});
				push(@SERNA,$INTA{$locus}{SEND});
				$icount=1;
			}else{
				#print SUM "\t";
				$smallwindow.="\t";
				$overlapwindow.="\t";
				$com.="\t\t\t\t\t\t\t";
				push(@MBRNA,"");push(@MERNA,"");push(@SBRNA,"");push(@SERNA,"");				
			}
		}else{
			#print SUM "\t\t";
			$largewindow.="\t";
			$smallwindow.="\t";
			$overlapwindow.="\t";
			$com.="\t\t\t\t\t\t\t";
			push(@MBRNA,"");push(@MERNA,"");push(@SBRNA,"");push(@SERNA,"");
		}

		
		if($COPRA{$locus}{GENE} && $copra_genomes && $rnas){ # Don't do it if not running CopraRNA
			#print SUM "1\t";
			$count1++;
			$largewindow.="1\t";
			if(($COPRA{$locus}{MSTA} <= $end && $COPRA{$locus}{MSTA} >= $begin) ||
			   ($COPRA{$locus}{MEND} <= $end && $COPRA{$locus}{MEND} >= $begin) ||
			   ($EXPRESS{$locus}=~/\w/ && (defined($list) || defined($expression))) ){
				#print SUM "1\t";
				$count2++;
				$smallwindow.="1\t";
				foreach $col (@KEYS){
					$com.="$COPRA{$locus}{$col}\t";
				}
				$COPRA{$locus}{STRUCTURE}=~s/\r$//;
				$com.="\"$COPRA{$locus}{STRUCTURE}\"\t";
				$com.="$COPRA{$locus}{CPVAL}\t$COPRA{$locus}{HOMOLOGS}\t$COPRA{$locus}{ANNO}\t";
				$mlen=$COPRA{$locus}{MEND} - $COPRA{$locus}{MSTA} + 1;
				$slen=$COPRA{$locus}{SEND} - $COPRA{$locus}{SSTA} + 1;
				$msta=$COPRA{$locus}{MSTA} + $upstream + 1;
				$mend=$COPRA{$locus}{MEND} + $upstream + 1;
				$stru=$COPRA{$locus}{STRUCTURE};
				$stru=~s/\r/\\n/g;
				$hybrid=$COPRA{$locus}{EVAL}*2;
				unless($COPRA{$locus}{GENE}=~/\w/){$gene=$locus;}
				else{$gene=$COPRA{$locus}{GENE};}
				
				$csv_line.="$COPRA{$locus}{PVAL};$COPRA{$locus}{PVAL};$gene\_c;$region;$mlen;$msta;";
				$csv_line.="$mend;$msta;$mend;na;$rna_name;$slen;$COPRA{$locus}{SSTA};$COPRA{$locus}{SEND};";
				$csv_line.="$COPRA{$locus}{SSTA};$COPRA{$locus}{SEND};na;$hybrid;$COPRA{$locus}{EVAL};$stru;$locus;na;na\n";
				
				$smallrna="";$rnatarget="";
				if($KNOWN{$locus}{FLAG} == 1){
					$smallrna=&MATCH($COPRA{$locus}{SSTA},$COPRA{$locus}{SEND},$KNOWN{$locus}{SSTA},$KNOWN{$locus}{SEND});
					$rnatarget=&MATCH($COPRA{$locus}{MSTA},$COPRA{$locus}{MEND},$KNOWN{$locus}{MSTA},$KNOWN{$locus}{MEND});
					if($smallrna > 0 && $rnatarget > 0){$overlapwindow.="A\t";}#print "A\n";}
					else{$overlapwindow.="B\t";}#print "B\n";}
				}else{$overlapwindow.="F\t";}#print "C\n";}
				push(@MBRNA,$COPRA{$locus}{MSTA});
				push(@MERNA,$COPRA{$locus}{MEND});
				push(@SBRNA,$COPRA{$locus}{SSTA});
				push(@SERNA,$COPRA{$locus}{SEND});
				$ccount=1;
			}else{
				#print SUM "\t";
				$smallwindow.="\t";
				$overlapwindow.="\t";
				$com.="\t\t\t\t\t\t\t\t\t";
				push(@MBRNA,"");push(@MERNA,"");push(@SBRNA,"");push(@SERNA,"");
			}
		}elsif($copra_genomes && $rnas){
			#print SUM "\t\t";
			$largewindow.="\t";
			$smallwindow.="\t";
			$overlapwindow.="\t";
			$com.="\t\t\t\t\t\t\t\t\t";
			push(@MBRNA,"");push(@MERNA,"");push(@SBRNA,"");push(@SERNA,"");			
		}
		# Cross reference all non-A matches assign letter groups: 
		# B->E are predictions not coincident with known, as many as 4 bins (forward compatable)
		# F->I are predictions when no known region, as many as 4 bins
		@LETTERS=split("\t",$overlapwindow);
		($better,$bin) = &SHARED(\@LETTERS,\@MBRNA,\@MERNA,\@SBRNA,\@SERNA);
		#print "$locus\t[$better][$bin]"; ##debug
		if($bin <= 4 || $bin==10){$rank=$bin;}
		elsif($ccount == 1 ){$rank=5;}
		elsif($scount == 1 || $tcount == 1){$rank=6;}
		elsif($icount == 1 ){$rank=7;}
		#print "[$rank]\n"; ##debug
		$sum_line.="$largewindow$count1\t$smallwindow$count2\t$better$rank\n";
		$PRINT_ORDER{$locus} = $rank;
		$PRINT_SUM{$locus}=$sum_line;
		
		#print SOM "$locus\t$TOTALLOCUS{$locus}\t$EXPRESS{$locus}$com\n";
		## Already inside loop that accounted for $include_all
		if(($EXPRESS{$locus}=~/\w/ && (defined($list) || defined($expression))) 
			|| ($num_predictions =~ /[Nn]/ && $rank != 10)
			|| ($num_predictions =~ /[Yy]/ && $rank == 4 ) ){
			#IF ON A LIST
			#IF WANT ALL genes - NO BLANKS
			#IF ONLY WANT genes that HAVE 2 or more matching predictions
			if($EXPRESS{$locus}!~/\t/ && $include_all =~ /[Yy]/ && (defined($list) || defined($expression))){
				$spacer="\t\t";
			}else{
				$spacer=$EXPRESS{$locus};
			}#$EXPRESS{$locus};
			$PRINT_COM{$locus}="$locus\t$TOTALLOCUS{$locus}\t$spacer$com\n";
			$PRINT_CSV{$locus}=$csv_line;
		}#
		
		
	    }## end if defined in list OR no list used
	    
	    	#write comma seperated file for CopraRNA version of termClusterReport.pl
		# only 1 line per gene unlike intarna_websrv_tbl.csv used for diagrams
	    if( $EXPRESS{$locus}=~/\w/ && (defined($list) || defined($expression))  && $count1 > 0 ) {
			# included if listed and has >=1 match or no list and 
			print FUN "0.001,";
			print FDR "0.001,0.001,";
	    }elsif( $list eq "" && $expression eq "" && $rank == 4){ # Rank will be "" if not processed 4 = 2 or more agreeing algorithms
			# no list or expression but 2 or more AGREEING predictions within region -b to -e
			# Genes with disagreeing predictions or only 1 prediction 
			print FUN "0.001,";
			print FDR "0.001,0.001,";
	    }else{	## All other genes for background functional distribution
			print FUN "0.5,";
			print FDR "0.5,0.5,";
	    } # no matches or fewer than required
	    print FUN "$locus($TOTALLOCUS{$locus}|-100|0|1|83|1|83|";
	    print FDR "$locus($TOTALLOCUS{$locus}|-100|0|1|83|1|83|";
		#if(length($GENEID{$locus}) >=8){
		#	print FUN "GI:$GENEID{$locus}),,,annotation,\n";
		#	print FDR "GI:$GENEID{$locus}),,,annotation,\n";
		#}else{	
	    print FUN "GeneID:$GENEID{$locus}),,,annotation,\n";
	    print FDR "GeneID:$GENEID{$locus}),,,annotation,\n";
		#}
		
			
	}## end if locus not blank
	
}## end foreach


## Print to files based on rank 1 --> 10
foreach $locus (sort {$PRINT_ORDER{$a} <=> $PRINT_ORDER{$b}} keys %PRINT_ORDER){
	
	if($PRINT_COM{$locus} =~/\w/){
		print SUM $PRINT_SUM{$locus};
		print COM $PRINT_COM{$locus}; 
		print CSV $PRINT_CSV{$locus};
	}

}
close(SUM); close(CSV); close(SUM); close(FUN); close(FDR);

print "\n=========Final summary Plots======================\n";
$d=`date`;chomp($d);
print "[$d]\n";
##Final Region and Heatmap Plots
system("cp binding_RNAs.fa intarna_pars.txt $sub/");
system("cp $out\_intarna_websrv_tbl.csv $sub/intarna_websrv_table.csv");
chdir($sub);
$path=`pwd`;
print $path;

$PATH_HOME="/software/IntaRNA1.0.4";
$PATH_R="/usr/bin";

# make regions plots - R code stolen from IntaRNA1.0.4 wrapper script
print "\n$PATH_R/R --slave -f $PATH_HOME/plotting_script_intaRNA2.r \n\n";
system("$PATH_R/R --slave -f $PATH_HOME/plotting_script_intaRNA2.r");
system("/usr/bin/convert -density '300' -resize '700' -flatten -rotate 90 sRNA_regions.ps sRNA_regions.png");
system("/usr/bin/convert -density '300' -resize '700' -flatten -rotate 90 mRNA_regions.ps mRNA_regions.png");


# Run Copra enrichment on collated/selected results

print "\n$INTADIR/rerun_enrichment.pl $out\_composite_websrv_tbl.csv \n\n";
system("$INTADIR/rerun_enrichment.pl $out\_composite_websrv_tbl.csv ");

chdir($PWD);

unless($skip =~ /[Yy]/){ 
	#Cleanup unused not to be reused results/intermediate files
	$folder="STARPICKER";
	unless(-d $folder){mkdir($folder);}
	system("mv *binding_site* *mfe_GC_strict_seed.txt *seed_selected.txt *target_up* $folder");
	$folder="INTARNA";
	unless(-d $folder){mkdir($folder);}
	system("mv predictions.* termClusterReport.txt target_RNAs.fa $folder");
	#binding_RNAs.fa intarna_pars.txt intarna_websrv_table_truncated.csv
}

print "\n/scripts/make_sRNA_heatmap.pl $out $out\_summary.txt \n\n";
system( "/scripts/make_sRNA_heatmap.pl $out $out\_summary.txt ");


print "\n/scripts/make_sRNA_xlsx.pl $out $out\_complete.txt $out\_summary.txt \n\n";
system( "/scripts/make_sRNA_xlsx.pl $out $out\_complete.txt $out\_summary.txt ");



##send mail
if($email=~/\w/){
	@program=split("/",$0);
	system("/scripts/sendmail.pl $out $email $program[-1]");
}

print "\n=========Complete=================================\n";
$d=`date`;chomp($d);
print "[$d]\n";

## ========= SUBROUTINES ========= ##

#Subroutine to check for overlap of two sets of coordinates
#Works on + and - integers series as long as ascending
sub MATCH{
	
	my ($start,$stop,$begin,$end)= @_;
	my $flag="0";

	if($start >= $begin && $stop <= $end){
		#within region	
		$flag=1;
	}elsif($begin >= $start && $end <= $stop){
		#Spans entire region
		$flag=2;
	}elsif($start < $begin && $stop >= $begin && $stop <= $end){
		#overlaps 5' region
		$flag=3;	
	}elsif($start >= $begin && $start <= $end && $stop > $end ){
		#overlaps 3' region
		$flag=4;			
	}
	#print "$start,$stop and $begin,$end = $flag\n"; ##debug
	return($flag);	
		
}	

#Subroutine to check for overlaps of non-KNOWN predictions (As)
sub SHARED{
	my $params = shift;
	my @LETTERS= @$params;
	$params = shift;
	my @MBRNA = @$params;
	$params = shift;
	my @MERNA = @$params;
	$params = shift;
	my @SBRNA = @$params;
	$params = shift;
	my @SERNA = @$params;
	my @BS=();
	my @FS=();
	my $m="";my $s="";my $x="";my $W="";my $y="";my $z="";
	for $i (0..3){
		if($LETTERS[$i] eq "B"){
			push(@BS,$i);
		}
		if($LETTERS[$i] eq "F"){
			push(@FS,$i);
		}
	}
#Check for all category Bs that did not match known set of coordinates
#BCDE
	#print @BS . "\n";
	if(@BS >= 2){
		$m=&MATCH($MBRNA[$BS[0]],$MERNA[$BS[0]],$MBRNA[$BS[1]],$MERNA[$BS[1]]);
		$s=&MATCH($SBRNA[$BS[0]],$SERNA[$BS[0]],$SBRNA[$BS[1]],$SERNA[$BS[1]]);
		#print "0-1	[$m][$s]\n";
		if($m > 0 && $s >0){
			#don't do anything already marked the same	
		}else{
			$LETTERS[$BS[1]]="C";	
		}
	}
	if(@BS >= 3){
		$m=&MATCH($MBRNA[$BS[0]],$MERNA[$BS[0]],$MBRNA[$BS[2]],$MERNA[$BS[2]]);
		$s=&MATCH($SBRNA[$BS[0]],$SERNA[$BS[0]],$SBRNA[$BS[2]],$SERNA[$BS[2]]);
		#print "0-2	[$m][$s]\n";
		if($m > 0 && $s >0){
			#don't do anything already marked the same	
		}else{
			$w=&MATCH($MBRNA[$BS[1]],$MERNA[$BS[1]],$MBRNA[$BS[2]],$MERNA[$BS[2]]);
			$x=&MATCH($SBRNA[$BS[1]],$SERNA[$BS[1]],$SBRNA[$BS[2]],$SERNA[$BS[2]]);
			#print "1-2	[$w][$x]\n";
			if($w > 0 && $x >0){
				$LETTERS[$BS[2]]="C";	
			}else{
				$LETTERS[$BS[2]]="D";	
			}
		}
	}
	if(@BS == 4){
		$m=&MATCH($MBRNA[$BS[0]],$MERNA[$BS[0]],$MBRNA[$BS[3]],$MERNA[$BS[3]]);
		$s=&MATCH($SBRNA[$BS[0]],$SERNA[$BS[0]],$SBRNA[$BS[3]],$SERNA[$BS[3]]);
		#print "0-3	[$m][$s]\n";
		if($m > 0 && $s >0){
			#don't do anything already marked the same	
		}else{
			$w=&MATCH($MBRNA[$BS[1]],$MERNA[$BS[1]],$MBRNA[$BS[3]],$MERNA[$BS[3]]);
			$x=&MATCH($SBRNA[$BS[1]],$SERNA[$BS[1]],$SBRNA[$BS[3]],$SERNA[$BS[3]]);
			#print "1-3	[$w][$x]\n";
			if($w > 0 && $x > 0){
				$LETTERS[$BS[3]]="C";	
			}else{
				$y=&MATCH($MBRNA[$BS[2]],$MERNA[$BS[2]],$MBRNA[$BS[3]],$MERNA[$BS[3]]);
				$z=&MATCH($SBRNA[$BS[2]],$SERNA[$BS[2]],$SBRNA[$BS[3]],$SERNA[$BS[3]]);
				#print "2-3	[$w][$x]\n";
				if($y > 0 && $z >0){
					$LETTERS[$BS[3]]="D";
				}else{
					$LETTERS[$BS[3]]="E";
				}	
			}
		}
	}
#Check for all category Fs that did not match known set of coordinates
#FGHI
	if(@FS >= 2){
		$m=&MATCH($MBRNA[$FS[0]],$MERNA[$FS[0]],$MBRNA[$FS[1]],$MERNA[$FS[1]]);
		$s=&MATCH($SBRNA[$FS[0]],$SERNA[$FS[0]],$SBRNA[$FS[1]],$SERNA[$FS[1]]);
		#print "0-1	[$m][$s]\n";
		if($m > 0 && $s >0){
			#don't do anything already marked the same	
		}else{
			$LETTERS[$FS[1]]="G";	
		}
	}
	if(@FS >= 3){
		$m=&MATCH($MBRNA[$FS[0]],$MERNA[$FS[0]],$MBRNA[$FS[2]],$MERNA[$FS[2]]);
		$s=&MATCH($SBRNA[$FS[0]],$SERNA[$FS[0]],$SBRNA[$FS[2]],$SERNA[$FS[2]]);
		#print "0-2	[$m][$s]\n";
		if($m > 0 && $s >0){
			#don't do anything already marked the same	
		}else{
			$w=&MATCH($MBRNA[$FS[1]],$MERNA[$FS[1]],$MBRNA[$FS[2]],$MERNA[$FS[2]]);
			$x=&MATCH($SBRNA[$FS[1]],$SERNA[$FS[1]],$SBRNA[$FS[2]],$SERNA[$FS[2]]);
			#print "1-2	[$w][$x]\n";
			if($w > 0 && $x >0){
				$LETTERS[$FS[2]]="G";	
			}else{
				$LETTERS[$FS[2]]="H";	
			}
		}
	}	
	if(@FS == 4){
		$m=&MATCH($MBRNA[$FS[0]],$MERNA[$FS[0]],$MBRNA[$FS[3]],$MERNA[$FS[3]]);
		$s=&MATCH($SBRNA[$FS[0]],$SERNA[$FS[0]],$SBRNA[$FS[3]],$SERNA[$FS[3]]);
		#print "0-3	[$m][$s]\n";
		if($m > 0 && $s >0){
			#don't do anything already marked the same	
		}else{
			$w=&MATCH($MBRNA[$FS[1]],$MERNA[$FS[1]],$MBRNA[$FS[3]],$MERNA[$FS[3]]);
			$x=&MATCH($SBRNA[$FS[1]],$SERNA[$FS[1]],$SBRNA[$FS[3]],$SERNA[$FS[3]]);
			#print "1-3	[$w][$x]\n";
			if($w > 0 && $x > 0){
				$LETTERS[$FS[3]]="G";	
			}else{
				$y=&MATCH($MBRNA[$FS[2]],$MERNA[$FS[2]],$MBRNA[$FS[3]],$MERNA[$FS[3]]);
				$z=&MATCH($SBRNA[$FS[2]],$SERNA[$FS[2]],$SBRNA[$FS[3]],$SERNA[$FS[3]]);
				#print "2-3	[$w][$x]\n";
				if($y > 0 && $z >0){
					$LETTERS[$FS[3]]="H";
				}else{
					$LETTERS[$FS[3]]="I";
				}	
			}
		}
	}
	my $return_list="";
	for my $i (0..3){
		$return_list.="$LETTERS[$i]\t";
	}
	my $twoormore=10;
	#print "[$return_list]\n";
	my $counter = () = $return_list =~ m/A/g;
	if($counter >= 2 ){$twoormore="1";}
	elsif($counter == 1 ){$twoormore="2";}
	elsif($return_list=~/[BCDE]/){$twoormore="3";}
	else{
		my @POSSIBLE=("F","G","H","I");
		foreach $c (@POSSIBLE){
	 		my $counter1 = () = $return_list =~ m/$c/g;
	 		if($counter1 >= 2 && 4 < $twoormore){$twoormore="4";}
	 		elsif($counter1 == 1 && 5 < $twoormore){$twoormore="5";}
	 	}
	}
	return($return_list,$twoormore);
	
}#end subroutine

__END__


