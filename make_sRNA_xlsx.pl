#!/usr/bin/perl 
##
## 	Patrick Degnan
##	make_sRNA_xlsx.pl
##	Convert text output from spot.pl and make into xlsx file
## 

# http://search.cpan.org/~jmcnamara/Excel-Writer-XLSX/lib/Excel/Writer/XLSX.pm
use Excel::Writer::XLSX;

unless(@ARGV){die "\nUsage: $0  output_prefix *complete.txt *summary.txt\n\n";}

$prefix=$ARGV[0];
$out="$prefix.xlsx";
# Create a new Excel workbook
$workbook = Excel::Writer::XLSX->new($out);

$file=$ARGV[1];
$file2=$ARGV[2];
# Add a worksheet
$worksheet = $workbook->add_worksheet($file);
$worksheet2 = $workbook->add_worksheet($file2); 

$format_white= $workbook->add_format();
$format_white->set_bg_color('#FFFFFF');
$format_white->set_border(1);
$format_white->set_border_color('#A0A0A0');

$format_red= $workbook->add_format();
$format_red->set_bg_color('#FF9999');
$format_red->set_border(1);
$format_red->set_border_color('#A0A0A0');

$format_blue= $workbook->add_format();
$format_blue->set_bg_color('#99CCFF');
$format_blue->set_border(1);
$format_blue->set_border_color('#A0A0A0');

$format_green= $workbook->add_format();
$format_green->set_bg_color('#78C98E');
$format_green->set_border(1);
$format_green->set_border_color('#A0A0A0');

$format_purple= $workbook->add_format();
$format_purple->set_bg_color('#FFCCFF');
$format_purple->set_border(1);
$format_purple->set_border_color('#A0A0A0');

$format_orange= $workbook->add_format();
$format_orange->set_bg_color('#FFCC99');
$format_orange->set_border(1);
$format_orange->set_border_color('#A0A0A0');

 
#  Add and define a format(s)	
$format_bold = $workbook->add_format();
$format_bold->set_bold();

$courier = $workbook->add_format();
$courier->set_font("Courier");
$courier->set_size(10);
$courier->set_text_wrap();	

#set_row( $row, $height, $format, $hidden, $level, $collapsed )
#set_column( $first_col, $last_col, $width, $format, $hidden, $level, $collapsed )
	
#If you wish to set the format without changing the height you can pass undef as the height parameter:
$worksheet->set_row( 0, undef, $format_bold );

# Row and column are zero indexed
$row = 0;
$total=0;

open(IN,"<$file") or die "cannot open $file\n";
while ( <IN> ) {
	chomp;
	# Split on single tab
	@fields = split( '\t', $_ );

	if($row==0){$total=@fields;}

	$col = 0;
	foreach $token ( @fields ) {
		$token=~s/\"//g;
		$worksheet->write( $row, $col, $token );
		$col++;
	}
	$row++;
}close(IN);   

if($total >= 35 || $total == 25){
	@COLS=(10,17,24,31);
}else{
	@COLS=(8,15,22,29);
}

foreach $c (@COLS){
	$worksheet->set_column( $c, $c, 70, $courier );
}


$worksheet2->set_row( 0, undef, $format_bold );    
# Row and column are zero indexed
$row = 0;
$total=0;

$header=`head -1 $file2`;
@cols=split(/\t/,$header);
$total_columns=$#cols;
#print "[$total_columns]\n";
if($total_columns == 13 ){
	# 3 columns, no list
	%set=(10,1,11,1,12,1);
}elsif($total_columns ==15){
	# 3 columns, list
	%set=(12,1,13,1,14,1);
}elsif($total_columns ==16){
	# 4 columns, no list
	%set=(12,1,13,1,14,1,15,1);
}else{
	# 4 columns, list
	%set=(14,1,15,1,16,1,17,1);
}

%COLORS=("",$format_white,"A",$format_red,"B",$format_blue,"C",$format_green,"D",$format_green,"E",$format_green,"F",$format_purple,"G",$format_orange,"H",$format_orange,"I",$format_orange);


open(IN,"<$file2") or die "cannot open $file2\n";
while ( <IN> ) {
	chomp;
	# Split on single tab
	@fields = split( '\t', $_ );

	if($row==0){$total=@fields;}

	$col = 0;
	foreach $token ( @fields ) {
		$token=~s/\"//g;
		if($set{$col} == 1 && $row != 0){
			#print "[$row,$col][$set{$col}][$token][$COLORS{$token}]\n";
			$format=$COLORS{$token};
			$worksheet2->write( $row, $col, $token, $format );
		}else{
			$worksheet2->write( $row, $col, $token );
		}
		
		
		$col++;
	}
	$row++;
}close(IN);   

$workbook->close();
