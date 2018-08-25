#!/usr/bin/perl
##
## 	Patrick Degnan
##	gbk_split2.pl
##	Split multi gb/gbk/gbff files intoindividual *gb files
##	Does not use BioPerl
##
unless($ARGV[0]){print "Convert concatenated *.gbk/*.gbff/*.gb files into single files\n";}

$/="";

if(@ARGV){                     #read in all *.br files one at a time 
	foreach $f (@ARGV){
		open(IN, "<$f") or die "Cannot open $f\n";
		$all=<IN>;
		@entries=split(/\n\/\//,$all);
	#	$name=$ARGV[0];
	#	$name=~ s/.gbk//;
		$count=1;
		foreach $e (@entries){
			if($e=~ /LOCUS\s+(\S+)\s+/){
				#$fileout=$1 . ".gbk";
				$name=$1;
				$fileout="$name.gb";
				open (OUT,">$fileout");
				$count++;
				$e=~s/^\n+//;
				#$e=~s/\/\///;
				print OUT "$e\n//\n";
				close(OUT);
				print "$fileout\n";
			}
		}
		#$command="mv $f ORIGINAL";
		#system($command);
	}	
}
 
