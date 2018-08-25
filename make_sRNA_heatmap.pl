#!/usr/bin/perl
##
## 	Patrick Degnan
##	make_sRNA_heatmap.pl
##	Convert text output from spot.pl and make into R heatmap summary of sRNA-mRNA matches
## 

unless(@ARGV){die "\nUsage: make_sRNA_heatmap.pl output_prefix *summary.txt\n\n";}

print "Make heatmap pdf(s) of summary.txt file for up to the 1st 96 entries\n";

## basic reading in of file, and splitting lines
$file=$ARGV[1];
$prefix=$ARGV[0];
open(IN,"<$file") or die "cannot open $file\n";
@lines=<IN>;
$header=shift(@lines);
#print "[$header]\n";
@cols=split(/\t/,$header);
$total_columns=$#cols;
$total_rows=$#lines;
#print "[$total_columns][$total_rows]\n";
if($total_columns == 13 || $total_columns ==15){
	# only three columns
	@set=(-4,-3,-2);
	$head="\tt\ts\ti\n";
}else{
	# four columns
	@set=(-5,-4,-3,-2);
	$head="\tt\ts\ti\tc\n";
}

if($total_rows < 96){$last=$total_rows;}
else{$last=95;}
#print "\$last = [$last]\n";
$count=0;
$merge_command="gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$prefix\_summary.pdf ";

for $i (0..$last){
	$l=$lines[$i];
	@cols=split(/\t/,$l);
	##add code
	if($i==0 || $i==48){
		$count++;
		open(OUT,">temp_table$count.txt") or die "Cannot open temp_table$count.txt\n";
		print OUT $head;
	}
	print OUT "'$cols[0] $cols[1]'";
	foreach $j (@set){
		if($cols[$j] eq ""){
			print OUT "\t0";
		}else{
			$cols[$j]=~tr/ABCDEFGHI/123334555/;
			print OUT "\t$cols[$j]";
		}
 	}
 	print OUT "\n";
 	if($i == $last || $i ==47){
 		$columns=@set;
		if($columns==4){print OUT "blank\t0\t0\t0\t5\n";}
		else{print OUT "blank\t0\t0\t5\n";}
		close(OUT);
		
		if($i > 47){$row=$i+2-48;}
		else{$row=$i+2;}
		#print "[$i][$row]\n";
	 	&WRITER("temp_table$count.txt","temp_table$count.pdf",$row,$columns);
	 	if($last > 47){
	 		$merge_command.="temp_table$count.pdf ";
	 	}
	 	
 	}
 	
}close(IN);

## File renaming/merging and deletion of intermediates
if($count > 1){
	`$merge_command`;
}else{
	`mv temp_table$count.pdf $prefix\_summary.pdf`;
}
`rm temp_table* temp.r`;

## write r script for heatmap
sub WRITER {

	my ($table,$output,$rows,$cols)=@_;

	my $cvar=54-$rows;
	my $rvar;
	if($cols==4){$rvar=35}
	else{$rvar=36.75}
	#print "[$rows][$cvar][$rvar]\n";
	open(ROUT,">temp.r") or die "Cannot open temp.r\n";
	print ROUT "library(gplots)
x=read.table(\"$table\", header=TRUE)
mat=data.matrix(x)

pdf(\"$output\", height=11, width=8.5)
heatmap.2(mat,
Rowv=NULL,
Colv=NULL,
dendrogram= c(\"none\"),
distfun = dist,
hclustfun = hclust,
key=TRUE,
keysize=.9,
trace=\"none\",
density.info=c(\"none\"),
margins=c($cvar, $rvar), 
col=c(
'0' = \"#FFFFFF\",
'1' = \"#FF9999\",
'2' = \"#99CCFF\",
'3' = \"#78C98E\",
'4' = \"#FFCCFF\",
'5' = \"#FFCC99\"),
sepwidth=c(0.005,0.005),
sepcolor=\"#A0A0A0\",
colsep=1:ncol(x),
rowsep=1:nrow(x),
srtCol=0,
adjCol = c(1,1),
)

dev.off()
";

	close(ROUT);
	` /usr/bin/R --vanilla <temp.r`;
	return();
}
