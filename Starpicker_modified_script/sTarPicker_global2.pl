#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;

sub extract_startcodon_with_up_down;
sub seed_accessibility_mfe_GC_strict_up350_global;
sub select_seed_with_hybrid_open_global;
sub search_binding_site_up350_global;
sub compute_feature_global;
sub classify_global;
sub output_with_pairwise_alignment_global;

sub energy_open_seed;
sub unstructured_energy_seed;
sub energy_open_binding_site;
sub unstructured_energy_binding_site;
sub basecontent;
sub pairwise_alignment;
sub in_array;
my $STARDIR="/software/starpicker";##
use constant USAGE=><<END;
   This is sTarPicker program for global target prediction of bacterial small RNAs.
   
   Usage: 
   $0 -srna sRNA sequences in fasta format -gen genome sequences in fasta format -ptt ptt files in ptt format -th the threshold for sTarPicker
     
   Parameters:
   -srna         sRNA sequences in fasta format, supporting multiple sRNA sequences;   
   -gen          genome sequences in fasta format, which can be downloaded from GenBank;
   -ptt          ptt files in ptt format, which can be downloaded from GenBank;
   -th    the threshold for sTarPicker, ranging from 0.001 to 1.
   -u		number of nt upstream of start site (default = 150)
   -d		number of nt downstream of start site (default = 100)
   -seed	seed size (default = 5)
   Example:
   $0 -srna Yfr1.fasta -gen NC_005072.fna -ptt NC_005072.ptt -th 0.5

END

die USAGE if scalar @ARGV < 8;
my ($srna_fasta, $genome_file, $ptt_file, $threshold,$upstream,$downstream, $seed_size);##
GetOptions("srna=s"=>\$srna_fasta, "gen=s"=>\$genome_file, "ptt=s"=>\$ptt_file, "th=f"=>\$threshold,"up=s"=>\$upstream,"down=s"=>\$downstream, "seed=s"=>\$seed_size);

unless($upstream){$upstream=150;}##
unless($downstream){$downstream=100;}##
unless($seed_size){$seed_size=5;}##

#my $genome_file_prefix = substr $genome_file, 0, (index $genome_file, '.'); 
$genome_file=~/(N\S+)\./;##
my $genome_file_prefix =$1;##

# make hash table of gene name and gene information
my %gene_product = ();
open IN, "<$ptt_file";
my @lines = <IN>;
close IN;
shift @lines; 
shift @lines;
shift @lines; # the first three lines are expelaination.
foreach my $line (@lines) {
	chomp $line;
	my ($position, $strand, $aa_len, $pid, $gene_name, $Synonym, $code, $cog, $product) = split /\t/, $line;
	my ($gene_full_name);
	if ($gene_name eq '-') {
		$gene_full_name = $Synonym;
	}
	else {	
		$gene_full_name = $gene_name . '_' . $Synonym;
	}
	$gene_product{$gene_full_name} = $product;
}

my $seq_in = Bio::SeqIO->new(-file => "<$srna_fasta", -format => 'fasta');
while (my $seqobj = $seq_in->next_seq()) {
	my $srna_name = $seqobj->display_id();
	my $srna_seq = $seqobj->seq();
	$srna_seq =~ s/\s//g;
	my $ups=350;##
	my $downs=300;##	
	my $prefix = $srna_name . '_' . $genome_file_prefix; 
	
	# Step 1: Input target genome and ptt file, output two datasets with upstream 150nt and downstream 100nt of targets and 350nt and 300nt of targets
	my $target_up150_down100 = $prefix . "_target_up" . $upstream . "_down" . $downstream .".fa"; 
	my $target_up350_down300 = $prefix . "_target_up" . $ups . "_down" . $downs .".fa"; 
	&extract_startcodon_with_up_down($genome_file, $ptt_file, $upstream, $downstream, $target_up150_down100);##

	&extract_startcodon_with_up_down($genome_file, $ptt_file, $ups, $downs, $target_up350_down300);##
	
	# Step 2: search seeds, preliminary filter seeds, and calculate energy of seeds, output seed file with energy
	my $seed_file = $prefix . '_mfe_GC_strict_seed.txt';
	&seed_accessibility_mfe_GC_strict_up350_global($srna_name, $srna_seq, $target_up150_down100, $target_up350_down300, $seed_file);
	
	# Step 3:  filter seeds according to the threshold, output final seeds
	my $seed_final = $prefix . '_seed_selected.txt';
	&select_seed_with_hybrid_open_global($seed_file, -2, $seed_final);
	
	# Step 4: find binding sites according to seed
	my $binding_site = $prefix . '_binding_site.txt';
	&search_binding_site_up350_global($srna_name, $srna_seq, $target_up150_down100, $target_up350_down300, $seed_final, $binding_site);
	
	# Step 5: calculate features of binding site, output feature matrix
	my $feature = $prefix . '_binding_site.feature';
	&compute_feature_global($srna_name, $srna_seq, $binding_site, $feature);
	
	# Step 6: Classify binding site with the model and select the binding site with the probability > 0.
	my $final_binding_site = $prefix . '_binding_site_selected.txt';
	&classify_global($feature, $binding_site, $threshold, $final_binding_site);
	
	# Step 7: Output final result with pairwise alignment
	my $final_binding_site_sorted = $prefix . '_binding_site_selected_sorted.txt';
	my $final_output = $prefix . '.output';
	&output_with_pairwise_alignment_global($srna_name, $srna_seq, $target_up150_down100, $target_up350_down300, $final_binding_site, $final_binding_site_sorted, $final_output);
}
$seq_in->close();

exit;


###############################################################################################
#                                                                                             #
#                                           Subroutines                                       #
#                                                                                             #
###############################################################################################
sub extract_startcodon_with_up_down {
	my ($sub_genome, $sub_ptt, $sub_upstream, $sub_downstream, $sub_output_target_up_down) = @_;

	open IN, "<$sub_ptt";
	my @sub_lines = <IN>;
	close IN;
	my $info = shift @sub_lines; # contain genome length
	shift @sub_lines;
	shift @sub_lines; # the first three lines are expelaination.
	
	# extract genome length from the $info line, i.e. first line
	chomp $info;
	# my ($genome_len) = $info =~ /\.\.(\d+)$/;
	my ($temp_temp_str, $genome_len) = split /\.\./, $info;

	open OUT, ">>$sub_output_target_up_down";
	my $seqin = Bio::SeqIO->new(-file => "<$sub_genome", -format => 'largefasta');
	my $seqobj = $seqin->next_seq();
	foreach my $line (@sub_lines) {
		my ($position, $strand, $aa_len, $pid, $gene_name, $Synonym) = split /\t/, $line;
		my ($gene_start, $gene_stop) = split /\.+/, $position;
		if ($gene_stop > $gene_start) {
			my $gene_len = $gene_stop - $gene_start + 1;
			my ($new_start, $new_stop);
			my ($real_upstream, $real_downstream);
			if ($strand eq '+') {
				$new_start = (($gene_start - $sub_upstream) > 0) ? ($gene_start - $sub_upstream) : 1;
				$new_stop = ($gene_len >= $sub_downstream) ? ($gene_start + $sub_downstream - 1) : $gene_stop;
				$real_upstream = $gene_start - $new_start;
				$real_downstream = $new_stop - $gene_start + 1;
			}
			else {
				$new_start = ($gene_len >= $sub_downstream) ? ($gene_stop - $sub_downstream + 1) : $gene_start;
				$new_stop = (($gene_stop + $sub_upstream) <= $genome_len) ? ($gene_stop + $sub_upstream) : $genome_len;
				$real_upstream = $new_stop - $gene_stop;
				$real_downstream = $gene_stop - $new_start + 1;
			}
			my $updown = $seqobj->subseq($new_start, $new_stop);
			if ($strand eq '-') {
				$updown = reverse $updown;
				$updown =~ tr/AGCTagct/TCGAtcga/;
			}
			
			my ($gene_name_full);
			if ($gene_name eq '-') {
				$gene_name_full = $Synonym . '_-' . $real_upstream . '_+' . $real_downstream;
			}
			else {
				$gene_name_full = $gene_name . '_' . $Synonym . '_-' . $real_upstream . '_+' . $real_downstream;
			}
			print OUT ">$gene_name_full\n";
			print OUT "$updown\n";
		}
	}
	$seqin->close();
	close OUT;
	
	return 0;
}

sub seed_accessibility_mfe_GC_strict_up350_global {
	my ($sub_srna_name, $sub_srna_seq, $up150_down100, $up350_down300, $sub_output) = @_;
	
	my $flanking = 100; # extract flanking 100nt sequence around seed

	# make hash table of target_up350_down300.
	# make hash table of target_name and upstream downstream for upstream350 and downstream 300.
	my (%target_up350_down300_hash);
	my (%target_name_upstream_up350, %target_name_downstream_up350);
	my $seqin = Bio::SeqIO->new(-file => "<$up350_down300", -format => 'fasta');
	while (my $seqobj = $seqin->next_seq()) {
		my $id = $seqobj->display_id();
		my ($new_id) = $id =~ /^(.*)_-\d+/;if($new_id eq ""){print "[$id]\n";}
		my ($temp_upstream) = $id =~ /_-(\d+)_/;
		my ($temp_downstream) = $id =~ /_\+(\d+)$/;
		my $seq = $seqobj->seq();
		$target_up350_down300_hash{$new_id} = $seq;
		$target_name_upstream_up350{$new_id} = $temp_upstream;
		$target_name_downstream_up350{$new_id} = $temp_downstream;
	}
	$seqin->close();
	
	# make hash table of target_name and upstream downstream for upstream350 and downstream 300.
	my (%target_name_upstream_up150, %target_name_downstream_up150);
	$seqin = Bio::SeqIO->new(-file => "<$up150_down100", -format => 'fasta');
	while (my $seqobj = $seqin->next_seq()) {
		my $id = $seqobj->display_id();
		my ($new_id) = $id =~ /^(.*)_-\d+/;
		my ($temp_upstream) = $id =~ /_-(\d+)_/;
		my ($temp_downstream) = $id =~ /_\+(\d+)$/;
		$target_name_upstream_up150{$new_id} = $temp_upstream;
		$target_name_downstream_up150{$new_id} = $temp_downstream;
	}
	$seqin->close(); 

	$sub_srna_name =~ s/\s//g;
	# calculate the structured energy
	my (@rnafold_results, $energy_start, $energy_end, $srna_energy_structured);
	@rnafold_results = `echo $sub_srna_seq|RNAfold -noPS`;
	$energy_start = rindex $rnafold_results[1], '(';
	$energy_end = rindex $rnafold_results[1], ')';
	$srna_energy_structured = substr $rnafold_results[1], ($energy_start + 1), ($energy_end - 1 - $energy_start);
	$srna_energy_structured =~ s/\s//g;

	my (%srna_start_end_open_energy);
	my $temp_target_for_guugle_dir = 'temp_target_for_guugle.fa';
	my $temp_query_for_guugle_dir = 'temp_query_for_guugle.fa';	
	my $temp_out_dir = 'temp_out.guugle';
	# find seed using guugle and calculate energy
	open OUTTMP, ">$temp_query_for_guugle_dir";
	print OUTTMP ">$sub_srna_name\n";
	print OUTTMP "$sub_srna_seq\n";
	close OUTTMP;

	open OUT, ">$sub_output";
	print OUT "target_name\tsRNA_seed_start\tsRNA_seed\ttarget_seed_start\ttarget_seed\tsRNA_seed_open\ttarget_seed_open\thybrid\thybrid_open\n";
	$seqin = Bio::SeqIO->new(-file => "<$up150_down100", -format => 'fasta');
	while (my $seqobj = $seqin->next_seq()) {
		my $sub_target_name_temp = $seqobj->display_id();
		my ($sub_target_name) = $sub_target_name_temp =~ /^(.*)_-\d+/;
		my $sub_target_seq_up150_down100 = $seqobj->seq();
		
		my (%target_start_end_open_energy);
		
		# run guugle to find seed
		open OUTTMP, ">$temp_target_for_guugle_dir";
		print OUTTMP ">$sub_target_name\n";
		print OUTTMP "$sub_target_seq_up150_down100\n";
		close OUTTMP;
		system ("$STARDIR/guugle -d $seed_size $temp_target_for_guugle_dir $temp_query_for_guugle_dir >$temp_out_dir");
		# filter seed
		open IN, "<$temp_out_dir";
		my @guugle_lines = <IN>;
		close IN;
		my ($seed_length, $target_start, $srna_start, $target_seed, $srna_seed);
		foreach my $guugle_line (@guugle_lines) {
			if ($guugle_line =~ /^MatchLength/) {
				($seed_length) = $guugle_line =~ /:\s*(\d+)\s*/;
				($target_start) = $guugle_line =~ /at\s*(\d+)\s*vs/;
				($srna_start) = $guugle_line =~ /at\s*(\d+)\s*\n/;
			}
			elsif ($guugle_line =~ /^5/) {
				($target_seed) = $guugle_line =~ /5(.*)3/;
			}
			elsif ($guugle_line =~ /^3/) {
				($srna_seed) = $guugle_line =~ /3(.*)5/;
			
				# count GU pairs in 5' end
				my $GU_5_primer = 0;
				for (my $i = 0; $i < $seed_length; $i++) {
					my $temp_target_nucleotide = substr $target_seed, $i, 1;
					my $temp_srna_nucleotide = substr $srna_seed, $i, 1;
					if ((($temp_srna_nucleotide eq 'g') and ($temp_target_nucleotide eq 'u')) or (($temp_srna_nucleotide eq 'u') and ($temp_target_nucleotide eq 'g'))) {
						$GU_5_primer++;
					}
					else {
						last;
					}
				}
				# count GU pairs in 3' end
				my $GU_3_primer = 0;
				for (my $i = ($seed_length - 1); $i >= 0; $i--) {
					my $temp_target_nucleotide = substr $target_seed, $i, 1;
					my $temp_srna_nucleotide = substr $srna_seed, $i, 1;
					if ((($temp_srna_nucleotide eq 'g') and ($temp_target_nucleotide eq 'u')) or (($temp_srna_nucleotide eq 'u') and ($temp_target_nucleotide eq 'g'))) {
						$GU_3_primer++;
					}
					else {
						last;
					}
				}
				# remove ending GU pairs
				my ($new_seed_length, $new_target_seed, $new_srna_seed, $new_target_start, $new_srna_start);
				$new_seed_length = $seed_length - $GU_5_primer - $GU_3_primer;
				$new_target_seed = substr $target_seed, $GU_5_primer, $new_seed_length;
				$new_srna_seed = substr $srna_seed, $GU_5_primer, $new_seed_length;
				$new_target_start = $target_start + $GU_5_primer;
				$new_srna_start = $srna_start + $GU_3_primer;
				
				# filter seed
				my $flag = 0;
				 
				if ($new_seed_length >= 5) {
					# count GU pairs and GC pairs
					my $GU_pairs = 0;
					my $GC_pairs = 0;
					my $AU_pairs = 0;
					my (@Cs_srna_seed) = $new_srna_seed =~ /(c)/g;
					my (@Cs_target_seed) = $new_target_seed =~ /(c)/g;
					$GC_pairs = (scalar @Cs_srna_seed) + (scalar @Cs_target_seed);
					my (@As_srna_seed) = $new_srna_seed =~ /(a)/g;
					my (@As_target_seed) = $new_target_seed =~ /(a)/g;
					$AU_pairs = (scalar @As_srna_seed) + (scalar @As_target_seed);
					$GU_pairs = $new_seed_length - $GC_pairs - $AU_pairs;
				
					if ($new_seed_length == 5) {
						if ($GU_pairs == 0) {
							if ($GC_pairs >= 4) {
								$flag = 1;
							}
							elsif ($GC_pairs >= 3) {
								if ($new_srna_seed =~ /.*[cg][cg][cg].*/) {
									$flag = 1;
								}
							}
						}
					}
					elsif (($new_seed_length > 5) and ($new_seed_length < 8)) {
						if (($GU_pairs <= 1) and ($GC_pairs >= 3)) {
							$flag = 1;
						}
					}
					elsif (($new_seed_length >= 8) and ($new_seed_length < 10)) {
						if (($GU_pairs <= 1) and ($GC_pairs >= 2)) {
							$flag = 1;
						}
					}
					elsif ($new_seed_length >= 10) {
						if (($GU_pairs <= 4) and ($GC_pairs >= 2)) {
							$flag = 1;
						}
					}
				}	
				# calculate energy
				if ($flag) {
					# calculate srna seed opening energy
					my $srna_seed_start = $new_srna_start;
					my $srna_seed_end = $new_srna_start + $new_seed_length - 1;
					
					my ($srna_energy_unstructured);
					my $srna_open_energy_key = $srna_seed_start . '_' . $srna_seed_end;
					if (defined $srna_start_end_open_energy{$srna_open_energy_key}) {
						$srna_energy_unstructured = $srna_start_end_open_energy{$srna_open_energy_key};
					}
					else {
						$srna_energy_unstructured = &unstructured_energy_seed($sub_srna_seq, $srna_seed_start, $srna_seed_end);
						$srna_start_end_open_energy{$srna_open_energy_key} = $srna_energy_unstructured;
					}
					my $open_srna_energy = $srna_energy_unstructured - $srna_energy_structured;
					$open_srna_energy = sprintf("%6.2f",$open_srna_energy);
				
					# calculate target seed opening energy		
					my ($open_target_energy);
					my $new_target_seed_end = $new_target_start + $new_seed_length - 1;
					my $target_open_energy_key = $new_target_start . '_' . $new_target_seed_end;
					if (defined $target_start_end_open_energy{$target_open_energy_key}) {
						$open_target_energy = $target_start_end_open_energy{$target_open_energy_key};
					}
					else {										
						my $target_up350_down300_seq = $target_up350_down300_hash{$sub_target_name};
						my $target_up350_down300_len = length $target_up350_down300_seq;
						my $new_target_start_in_up350 = $new_target_start + ($target_name_upstream_up350{$sub_target_name} - $target_name_upstream_up150{$sub_target_name}); # up350_down300 and up150_down100		
						my $new_target_start_flanking = ($new_target_start_in_up350 > $flanking) ? ($new_target_start_in_up350 - $flanking) : 1; # extract upstream 100 and downstream 100 around seed 
						my $new_target_end_flanking = (($new_target_start_in_up350 + $new_seed_length - 1 + $flanking) < $target_up350_down300_len) ? ($new_target_start_in_up350 + $new_seed_length - 1 + $flanking) : $target_up350_down300_len;
						my $flanking_seq = substr $target_up350_down300_seq, ($new_target_start_flanking - 1), ($new_target_end_flanking - $new_target_start_flanking + 1);
						my $temp_target_seed_start = $new_target_start_in_up350 - $new_target_start_flanking + 1;
						my $temp_target_seed_end = $temp_target_seed_start + $new_seed_length - 1;
						$open_target_energy = &energy_open_seed($flanking_seq, $temp_target_seed_start, $temp_target_seed_end);
						$target_start_end_open_energy{$target_open_energy_key} = $open_target_energy;
					}	
					$open_target_energy = sprintf("%6.2f",$open_target_energy);
				
					# calculate hybrid energy
					my $new_srna_seed_reverse = reverse $new_srna_seed;
					my @templines = ` echo "$new_srna_seed_reverse\n$new_target_seed\n"|RNAduplex `;
					my $energy_start = rindex $templines[0], '('; 
					my $energy_end = rindex $templines[0], ')'; 
					my $hybrid = substr $templines[0], ($energy_start + 1), ($energy_end - 1 - $energy_start); 
					$hybrid =~ s/\s//g;

					# calculate the energy difference
					my $hybrid_open = $hybrid + $open_srna_energy + $open_target_energy;
					$hybrid_open = sprintf("%6.2f",$hybrid_open);
					# output
					print OUT "$sub_target_name\t$new_srna_start\t$new_srna_seed\t$new_target_start\t$new_target_seed\t$open_srna_energy\t$open_target_energy\t$hybrid\t$hybrid_open\n";
				}
			}
		}
	}
	unlink $temp_out_dir;
	unlink $temp_target_for_guugle_dir;
	unlink $temp_query_for_guugle_dir;
	close OUT;
	
	return 0;
}

sub select_seed_with_hybrid_open_global {
	# $threshold is -1 if using -p0; $threshold is -2 if using mfe; 
	my ($in_seed_global, $sub_threshold, $out_seed_global) = @_;
	
	open IN, "<$in_seed_global";
	my @sub_lines = <IN>;
	close IN;
	shift @sub_lines; # the first line is title

	open OUT, ">$out_seed_global";
	print OUT "target_name\tsRNA_seed_start\tsRNA_seed\ttarget_seed_start\ttarget_seed\tsRNA_seed_open\ttarget_seed_open\thybrid\thybrid_open\n";
	foreach (@sub_lines) {
		chomp;
		my ($sub_target_name, $new_srna_start, $new_srna_seed, $new_target_start, $new_target_seed, $open_srna_energy, $open_target_energy, $hybrid, $hybrid_open) = split /\t/;
		if ($hybrid_open < $sub_threshold) {
			print OUT "$_\n";
		}
	}
	close OUT;
	
	return 0;
}

sub search_binding_site_up350_global {
	# output: Target_Name	Target_seq_up350_down300	sRNA_binding_site	Target_binding_site	sRNA_binding_seq	Target_binding_seq	Structure	Energy
	# binding start and end are delimited by comma
	# using target up350_down300 sequences to extract up and downstream 100nt arount seeds. 
	# However, the position of seeds provided by the seed_file is using target up150_down100.
	# So remember to convert the positions of seeds.
	my ($sub_srna_name, $sub_srna_seq, $up150_down100, $up350_down300, $seed_file, $out_sample) = @_;
	my $flanking = 100;

	# make hash table of target_up350_down300.
	# make hash table of target_name and upstream downstream for upstream350 and downstream 300.
	my (%target_up350_down300_hash);
	my (%target_name_upstream_up350, %target_name_downstream_up350);
	my $seqin = Bio::SeqIO->new(-file => "<$up350_down300", -format => 'fasta');
	while (my $seqobj = $seqin->next_seq()) {
		my $id = $seqobj->display_id();
		my ($new_id) = $id =~ /^(.*)_-\d+/;
		my ($temp_upstream) = $id =~ /_-(\d+)_/;
		my ($temp_downstream) = $id =~ /_\+(\d+)$/;
		my $seq = $seqobj->seq();
		$target_up350_down300_hash{$new_id} = $seq;
		$target_name_upstream_up350{$new_id} = $temp_upstream;
		$target_name_downstream_up350{$new_id} = $temp_downstream;
	}
	$seqin->close();
	
	# make hash table of target_name and upstream downstream for upstream350 and downstream 300.
	my (%target_name_upstream_up150, %target_name_downstream_up150);
	$seqin = Bio::SeqIO->new(-file => "<$up150_down100", -format => 'fasta');
	while (my $seqobj = $seqin->next_seq()) {
		my $id = $seqobj->display_id();
		my ($new_id) = $id =~ /^(.*)_-\d+/;
		my ($temp_upstream) = $id =~ /_-(\d+)_/;
		my ($temp_downstream) = $id =~ /_\+(\d+)$/;
		$target_name_upstream_up150{$new_id} = $temp_upstream;
		$target_name_downstream_up150{$new_id} = $temp_downstream;
	}
	$seqin->close(); 

	open IN, "<$seed_file";
	my @seed_lines = <IN>;
	close IN;
	shift @seed_lines; # the first line is title

	open OUT, ">$out_sample";
	print OUT "Target_Name\tTarget_seq_up350_down300\tsRNA_binding_site\tTarget_binding_site\tsRNA_binding_seq\tTarget_binding_seq\tStructure\tEnergy\tsRNA_seed_position\ttarget_seed_position_up150\tsRNA_seed\ttarget_seed\tOpen_sRNA_seed\tOpen_target_seed\tHybrid_seed\tHybrid_open_seed\n";
	foreach my $seed_line (@seed_lines) {
		chomp $seed_line;
		my ($sub_target_name, $srna_seed_start, $srna_seed, $target_seed_start_up150, $target_seed, $open_srna_seed_energy, $open_target_seed_energy, $seed_hybrid, $seed_hybrid_open) = split /\t/, $seed_line;
		my $seed_len = length $srna_seed;
	
		my $srna_seed_end = $srna_seed_start + $seed_len - 1;
		my $srna_len = length $sub_srna_seq;
		my $srna_restriction = "." x ($srna_seed_start - 1) . "(" x $seed_len . "." x ($srna_len - $srna_seed_end); 
		
		my $target_seed_start_up350 = ($target_name_upstream_up350{$sub_target_name} - $target_name_upstream_up150{$sub_target_name}) + $target_seed_start_up150;
		my $target_seed_end_up350 = $target_seed_start_up350 + $seed_len - 1;
		my $target_up350_down300_seq = $target_up350_down300_hash{$sub_target_name};
		my $target_up350_down300_seq_len = length $target_up350_down300_seq;
		my $target_fold_start = (($target_seed_start_up350 - $flanking) >= 0) ? ($target_seed_start_up350 - $flanking) : 1;
		my $target_fold_end = (($target_seed_end_up350 + $flanking) < $target_up350_down300_seq_len) ? ($target_seed_end_up350 + $flanking) : $target_up350_down300_seq_len;
		my $target_fold_len = $target_fold_end - $target_fold_start + 1;
		my $target_fold_seq = substr $target_up350_down300_seq, ($target_fold_start - 1), $target_fold_len;
		my $target_restriction = "." x ($target_seed_start_up350 - $target_fold_start) . ")" x $seed_len . "." x ($target_fold_end - $target_seed_end_up350); 
	
		my $seq = $sub_srna_seq . '&' . $target_fold_seq;
		my $restriction = $srna_restriction . '&' . $target_restriction;
		#print "$sub_target_name\n$seq\n$restriction\n";##
		my @ct_lines = ` echo "$seq\n$restriction\n"|RNAcofold -C -noLP -noCloseGU -noPS |$STARDIR/b2ct `;
		shift @ct_lines; # the first line is energy
		# make paring hash table
		my (%pair_hash);
		foreach my $ct_line (@ct_lines) {
			chomp $ct_line;
			#print "[$ct_line]\n";##
			my (@array) = split /\s+/, $ct_line;
			#print @array;##
			$pair_hash{$array[-1]} = $array[-2];
		}
	
		my $target_seed_start_ct = ($target_seed_start_up350 - $target_fold_start) + 1 + $srna_len + 1; # '&' in $srna_len + 1
		my $target_seed_end_ct = $target_seed_start_ct + $seed_len - 1;
		#print "[$target_seed_start_ct][$target_seed_end_ct]\n";	##
		my ($i_start, $i_end);
		# find sRNA binding start
		my $srna_binding_start = $srna_seed_start;
		$i_start = $srna_seed_start - 1;
		#print "$seed_line\n[$seq]\n[$restriction]\n";
		for (my $i = $i_start; $i >= 1; $i--) {
			my $value = $pair_hash{$i};
			if ($value > 0) {
				if ($value <= $srna_len) { 
					last; #form intra-structure, stop
				}
				else {
					$srna_binding_start = $i;
				}
			}
		}
		#die "Past first find sRNA binding start\n";##
		# find sRNA binding end
		my $srna_binding_end = $srna_seed_end;
		$i_start = $srna_seed_end + 1;
		for (my $i = $i_start; $i <= $srna_len; $i++) {
			my $value = $pair_hash{$i};
			if ($value > 0) {
				if ($value <= $srna_len) { 
					last; #form intra-structure, stop
				}
				else {
					$srna_binding_end = $i;
				}
			}
		}
	
		# find target binding start
		my $target_binding_start_ct = $target_seed_start_ct;
		$i_start = $target_seed_start_ct - 1;
		$i_end = $srna_len + 1;
		my $target_offset = $srna_len + 1;
		for (my $i = $i_start; $i > $i_end; $i--) {
			my $value = $pair_hash{$i};
			if ($value > 0) {
				if ($value > $target_offset) { 
					last; #form intra-structure, stop
				}
				else {
					$target_binding_start_ct = $i;
				}
			}
		}
	
		# find target binding end
		my $target_binding_end_ct = $target_seed_end_ct;
		$i_start = $target_binding_end_ct + 1;
		$i_end = $srna_len + $target_fold_len + 1;	
		for (my $i = $i_start; $i <= $i_end; $i++) {
			my $value = $pair_hash{$i};
			if ($value > 0) {
				if ($value > $target_offset) { 
					last; #form intra-structure, stop
				}
				else {
					$target_binding_end_ct = $i;
				}
			}
		}
	
		# trim unpaired nucleotides from binding_start_ct and binding_end_ct
		for (my $i = $srna_binding_end; $i >= $srna_binding_start; $i--) {
			my $value = $pair_hash{$i};
			if ($value > 0) {
				$srna_binding_end = $i;
				last;
			}
		}
		for (my $i = $srna_binding_start; $i <= $srna_binding_end; $i++) {
			my $value = $pair_hash{$i};
			if ($value > 0) {
				$srna_binding_start = $i;
				last;
			}
		}
		for (my $i = $target_binding_end_ct; $i >= $target_binding_start_ct; $i--) {
			my $value = $pair_hash{$i};
			if ($value > 0) {
				$target_binding_end_ct = $i;
				last;
			}
		}
		for (my $i = $target_binding_start_ct; $i <= $target_binding_end_ct; $i++) {
			my $value = $pair_hash{$i};
			if ($value > 0) {
				$target_binding_start_ct = $i;
				last;
			}
		}	
	
		# get the final sRNA and target binding site
		my $target_binding_start_pairing = $pair_hash{$srna_binding_end};
		my $target_binding_end_pairing = $pair_hash{$srna_binding_start};

		my $target_binding_start_final = ($target_binding_start_ct > $target_binding_start_pairing) ? $target_binding_start_ct : $target_binding_start_pairing;
		my $target_binding_end_final = ($target_binding_end_ct < $target_binding_end_pairing) ? $target_binding_end_ct : $target_binding_end_pairing;
		
		# final srna binding site
		my $srna_binding_start_final = $pair_hash{$target_binding_end_final};	
		my $srna_binding_end_final = $pair_hash{$target_binding_start_final};	
		
		# minus offset
		$target_binding_start_final = $target_binding_start_final - $target_offset;
		$target_binding_end_final = $target_binding_end_final - $target_offset;		
		
		my $target_binding_start_up350 = $target_seed_start_up350 + $target_binding_start_final - ($target_seed_start_up350 - $target_fold_start + 1);
		my $target_binding_end_up350 = $target_binding_start_up350 + ($target_binding_end_final - $target_binding_start_final);
	
		# calculate hybrid energy
		my $srna_binding_seq = substr $sub_srna_seq, ($srna_binding_start_final - 1), ($srna_binding_end_final - $srna_binding_start_final + 1);
		my $target_binding_seq = substr $target_up350_down300_seq, ($target_binding_start_up350 - 1), ($target_binding_end_up350 - $target_binding_start_up350 + 1);
		my @templines = ` echo "$srna_binding_seq\n$target_binding_seq\n"|RNAduplex `;
		chomp $templines[0];
		my ($structure, $srna_position, $colon, $target_position, $energy) = split /\s+/, $templines[0];
		my $energy_start = index $energy, '(';
		my $energy_end = index $energy, ')';
		my $hybrid = substr $energy, ($energy_start + 1), ($energy_end - 1 - $energy_start);
		$hybrid =~ s/\s//g;
		$hybrid = sprintf("%6.2f",$hybrid);
	
		my $srna_seed_reverse = reverse $srna_seed;
		my $target_seed_end_up150 = $target_seed_start_up150 + $seed_len - 1;
		print OUT "$sub_target_name\t$target_up350_down300_seq\t$srna_binding_start_final,$srna_binding_end_final\t$target_binding_start_up350,$target_binding_end_up350\t$srna_binding_seq\t$target_binding_seq\t$structure\t$hybrid\t$srna_seed_start,$srna_seed_end\t$target_seed_start_up150,$target_seed_end_up150\t$srna_seed_reverse\t$target_seed\t$open_srna_seed_energy\t$open_target_seed_energy\t$seed_hybrid\t$seed_hybrid_open\n";
	}
	close OUT;
	
	return 0;
}		

sub compute_feature_global {
	my ($sub_srna_name, $sub_srna_seq, $samples, $output) = @_;

	my $flanking = 100;

	# calculate the structured energy
	my (@rnafold_results, $energy_start, $energy_end, $srna_energy_structured);
	@rnafold_results = `echo $sub_srna_seq|RNAfold -p0 -noPS`;
	$srna_energy_structured = substr($rnafold_results[2],index($rnafold_results[2],"-"),index($rnafold_results[2],"k") - index($rnafold_results[2],"-") - 1);
	#print @rnafold_results;##
	
	open IN, "<$samples";
	my @sub_lines = <IN>;
	close IN;
	shift @sub_lines; # the first line is title
	open OUT, ">$output";
	foreach my $line (@sub_lines) {
		chomp $line;
		my ($sub_target_name, $target_seq_up350_down300, $srna_binding_site, $target_binding_site, $srna_binding_seq, $target_binding_seq, $structure, $energy, $srna_seed_position, $target_seed_position_up150, $srna_seed, $target_seed, $open_srna_seed_energy, $open_target_seed_energy, $seed_hybrid, $seed_hybrid_open) = split /\t/, $line;
		# calculate base content of hybridization
		my $hybrid_seq = $srna_binding_seq . $target_binding_seq;
		my ($content_array_ref) = &basecontent($hybrid_seq);
		my @content_array = @$content_array_ref;
		
		# calculate the energies of opening the structure of sRNA binding site and target binding site. Input sRNA and target sequences and binding site
		my ($srna_binding_start, $srna_binding_end) = split /,/, $srna_binding_site; 
		my $srna_energy_unstructured = &unstructured_energy_binding_site($sub_srna_seq, $srna_binding_start, $srna_binding_end);
		my $open_srna_energy = $srna_energy_unstructured - $srna_energy_structured;
		$open_srna_energy = sprintf("%6.2f",$open_srna_energy);
	
		my ($target_binding_start, $target_binding_end) = split /,/, $target_binding_site; 
		
		my $flanking_start_real = (($target_binding_start - $flanking) >= 1) ? $flanking : ($target_binding_start - 1);
		my $target_binding_start_flanking = $target_binding_start - $flanking_start_real; # extract upstream 100 and downstream 100 around seed 
		
		my $target_seq_up350_down300_len = length $target_seq_up350_down300;
		my $flanking_end_real = (($target_binding_end + $flanking) <= $target_seq_up350_down300_len) ? $flanking : ($target_seq_up350_down300_len - $target_binding_end);
		my $target_binding_end_flanking = $target_binding_end + $flanking_end_real;
		
		my $flanking_seq = substr $target_seq_up350_down300, ($target_binding_start_flanking - 1), ($target_binding_end_flanking - $target_binding_start_flanking + 1);
		my $temp_target_binding_start = $flanking_start_real + 1;
		my $temp_target_binding_end = $flanking_start_real + $target_binding_end - $target_binding_start + 1;
		#print "[$flanking_seq]\t[$temp_target_binding_start]\t[$temp_target_binding_end]\n";##
		my $open_target_energy = &energy_open_binding_site($flanking_seq, $temp_target_binding_start, $temp_target_binding_end);
		$open_target_energy = sprintf("%6.2f",$open_target_energy);

		my $energy_plus_srna_target = $energy + $open_target_energy + $open_srna_energy;
		$energy_plus_srna_target = sprintf("%6.2f",$energy_plus_srna_target);
	
		# ouput
		my (@final_features, $num);
		$num = scalar @content_array;
		# base content
		for (my $i = 0; $i < $num; $i++) {
			my $temp_feature = sprintf("%6.2f",$content_array[$i]);
			push @final_features, $temp_feature;
		}
		# energy and delta delta G of binding site
		push @final_features, $energy;
		push @final_features, $energy_plus_srna_target;
		# seed length
		my $real_seed_len = length $srna_seed;
		push @final_features, $real_seed_len;
		
		my $final_str = join "\t", @final_features;
		print OUT "$final_str\n";
	}
	close OUT;
	#die;##flag
	return 0;
}
	
	
sub classify_global {
	my ($feature_file, $binding_site_file, $sub_threshold, $output_binding_file) = @_;
	
	open IN, "<$feature_file";
	my @feature_lines = <IN>;
	close IN;
	my $line_num = scalar @feature_lines;
	
	open IN, "<$binding_site_file";
	my @binding_lines = <IN>;
	close IN;
	my $binding_title_line = shift @binding_lines; # The first line is title line.
	chomp $binding_title_line;
	
	my $classifyfile = "$STARDIR/model16m.m";
	open IN, "<$classifyfile" || die "Could not open file $classifyfile!\n"; 
	my @class_lines = <IN>;
	close IN;
	splice @class_lines, 0, 4; # The first four lines are not used.
	
	open OUT, ">$output_binding_file";
	print OUT "$binding_title_line\tHybrid_plus_open_both\tProbability\n";
		
	for (my $i = 0; $i < $line_num; $i++) {
		my $votes = 0;
		my @features = split /\t/, $feature_lines[$i];
		my $feature_num = scalar @features;
		
		for (my $j = 0; $j < 1999; $j = $j + 2) {
			my @pos_weights = split /\s+/, $class_lines[$j];
			shift @pos_weights; # the first element is blank;
			my @neg_weights = split /\s+/, $class_lines[$j + 1];
			shift @neg_weights; # the first element is blank;
			my $sum_pos = shift @pos_weights; # the first element is constant
			my $sum_neg = shift @neg_weights; # the first element is constant
			for (my $k = 0; $k < $feature_num; $k++) {
				$sum_pos += $pos_weights[$k] * $features[$k];
				$sum_neg += $neg_weights[$k] * $features[$k];
			}
			$votes++ if ($sum_pos >= $sum_neg);
		}
		$votes = $votes/1000;
		
		if ($votes >= $sub_threshold) {
			chomp $binding_lines[$i];
			print OUT "$binding_lines[$i]\t$features[-2]\t$votes\n";
		}
	}
	close OUT;
	
	return 0;
}

sub output_with_pairwise_alignment_global {
	my ($sub_srna_name, $sub_srna_seq, $sub_up150_down100, $sub_up350_down300, $sub_final_binding_site, $sub_final_binding_site_sorted, $sub_final_output) = @_;
	
	# make hash table of target_name and upstream downstream for upstream350 and downstream 300.
	my (%target_name_upstream_up350, %target_name_downstream_up350);
	my $seqin = Bio::SeqIO->new(-file => "<$sub_up350_down300", -format => 'fasta');
	while (my $seqobj = $seqin->next_seq()) {
		my $id = $seqobj->display_id();
		my ($new_id) = $id =~ /^(.*)_-\d+/;
		my ($temp_upstream) = $id =~ /_-(\d+)_/;
		my ($temp_downstream) = $id =~ /_\+(\d+)$/;
		$target_name_upstream_up350{$new_id} = $temp_upstream;
		$target_name_downstream_up350{$new_id} = $temp_downstream;
	}
	$seqin->close();
	
	# make hash table of target_name and upstream downstream for upstream350 and downstream 300.
	my (%target_name_upstream_up150, %target_name_downstream_up150);
	$seqin = Bio::SeqIO->new(-file => "<$sub_up150_down100", -format => 'fasta');
	while (my $seqobj = $seqin->next_seq()) {
		my $id = $seqobj->display_id();
		my ($new_id) = $id =~ /^(.*)_-\d+/;
		my ($temp_upstream) = $id =~ /_-(\d+)_/;
		my ($temp_downstream) = $id =~ /_\+(\d+)$/;
		$target_name_upstream_up150{$new_id} = $temp_upstream;
		$target_name_downstream_up150{$new_id} = $temp_downstream;
	}
	$seqin->close(); 
	
	open IN, "<$sub_final_binding_site";
	my @lines = <IN>;
	close IN;
	shift @lines; #The first line is title
	#Target_Name	Target_seq_up350_down300	sRNA_binding_site	Target_binding_site	sRNA_binding_seq	Target_binding_seq	Structure	Energy	sRNA_seed_position	target_seed_position_up150	sRNA_seed	target_seed	Open_sRNA_seed	Open_target_seed	Hybrid_seed	Hybrid_open_seed	Hybrid_plus_open_both	Probability
	
	my (%target_lines, %target_best_line);
	my $temp_target_name = '';
	my $line_num = scalar @lines;
	for (my $i = 0; $i < $line_num; $i++) {
		chomp $lines[$i];
		my ($Target_Name, $Target_seq_up350_down300, $sRNA_binding_site, $Target_binding_site, $sRNA_binding_seq, $Target_binding_seq, $Structure, $Energy, $sRNA_seed_position, $target_seed_position_up150, $sRNA_seed, $target_seed, $Open_sRNA_seed, $Open_target_seed, $Hybrid_seed, $Hybrid_open_seed, $Hybrid_plus_open_both, $Probability) = split /\t/, $lines[$i];
		if ($Target_Name eq $temp_target_name) {
			my $sub_array_ref = $target_lines{$Target_Name};
			my @sub_array = @$sub_array_ref;
			my $sub_flag = 0;
			my $sub_array_num = scalar @sub_array;
			for (my $j = 0; $j < $sub_array_num; $j++) {
				chomp $lines[$sub_array[$j]];
				my (@temp_array) = split /\t/, $lines[$sub_array[$j]];
				if ($Probability > $temp_array[-1]) {
					splice @sub_array, $j, 0, $i;
					$target_lines{$Target_Name} = \@sub_array;
					$sub_flag = 1;
					last;
				}
				elsif ($Probability == $temp_array[-1]) {
					if ($Hybrid_plus_open_both < $temp_array[-2]) {
						splice @sub_array, $j, 0, $i;
						$target_lines{$Target_Name} = \@sub_array;
						$sub_flag = 1;
						last;
					}
					elsif ($Hybrid_plus_open_both == $temp_array[-2]) {
						if ($Hybrid_open_seed < $temp_array[-3]) {
							splice @sub_array, $j, 0, $i;
							$target_lines{$Target_Name} = \@sub_array;
							$sub_flag = 1;
							last;
						}
					}
				}
			}
			if (!$sub_flag) {
				push @sub_array, $i;
				$target_lines{$Target_Name} = \@sub_array;
			}
			$target_best_line{$Target_Name} = $sub_array[0]; # The first one is the best binding site	
		}
		else {
			my @sub_array = ($i);
			$target_lines{$Target_Name} = \@sub_array;
			$target_best_line{$Target_Name} = $i;
			$temp_target_name = $Target_Name;
		}
	}
	
	my $temp_binding_site = 'temp_binding_site.txt';
	open OUT, ">$temp_binding_site";
	my @sub_keys = keys %target_best_line;
	foreach my $sub_key (@sub_keys) {
		my $temp_line_num = $target_best_line{$sub_key};
		print OUT "$lines[$temp_line_num]\n";
	}
	close OUT;
	
	my $temp_binding_site_sorted = 'temp_binding_site_sorted.txt';
	my $command = "cat $temp_binding_site | sort -k 18nr -k 17n -k 16n -o $temp_binding_site_sorted";
	system ($command);

	open IN, "<$temp_binding_site_sorted";
	my @best_lines = <IN>;
	close IN;
	open OUT, ">$sub_final_output";
	open OUTBIND, ">$sub_final_binding_site_sorted";
	my $i = 0;
	foreach my $best_line (@best_lines) {
		$i++;
		my ($sub_Target_Name) = split /\t/, $best_line;
		print OUT "$i\t$sub_Target_Name\t$gene_product{$sub_Target_Name}\n";
		my $sub_array_ref = $target_lines{$sub_Target_Name};
		my @sub_array = @$sub_array_ref;
		foreach my $sub_array_ele (@sub_array) {
			print OUTBIND "$lines[$sub_array_ele]\n";
			chomp $lines[$sub_array_ele];
			my ($Target_Name, $Target_seq_up350_down300, $sRNA_binding_site, $Target_binding_site, $sRNA_binding_seq, $Target_binding_seq, $Structure, $Energy, $sRNA_seed_position, $target_seed_position_up150, $sRNA_seed, $target_seed, $Open_sRNA_seed, $Open_target_seed, $Hybrid_seed, $Hybrid_open_seed, $Hybrid_plus_open_both, $Probability) = split /\t/, $lines[$sub_array_ele];
			$sRNA_binding_seq =~ s/T/U/g;
			$sRNA_binding_seq =~ s/t/u/g;
			$Target_binding_seq =~ s/T/U/g;
			$Target_binding_seq =~ s/t/u/g;
			my ($srna_binding_start, $srna_binding_end) = split /,/, $sRNA_binding_site;
			my ($target_binding_start_up350, $target_binding_end_up350) = split /,/, $Target_binding_site;
			my ($srna_seed_start, $srna_seed_end) = split /,/, $sRNA_seed_position;
			my ($target_seed_start_up150, $target_seed_end_up150) = split /,/,  $target_seed_position_up150;
			my $target_seed_start_up350 = ($target_name_upstream_up350{$Target_Name} - $target_name_upstream_up150{$Target_Name}) + $target_seed_start_up150;
			my $target_seed_end_up350 = ($target_name_upstream_up350{$Target_Name} - $target_name_upstream_up150{$Target_Name}) + $target_seed_end_up150;
			my $srna_len = length $sub_srna_seq;
			my $target_up350_down300_seq_len = length $Target_seq_up350_down300;
			
			# Re-predict structure using RNAcofold	
			# To ensure RNAcofold predicts all pairs, I add a dangling end to sRNA and target binding site
			my $srna_binding_start_dangling = ($srna_binding_start > 1) ? ($srna_binding_start - 1) : $srna_binding_start;
			my $srna_binding_end_dangling = ($srna_binding_end < $srna_len) ? ($srna_binding_end + 1) : $srna_binding_end;
			my $target_binding_start_up350_dangling = ($target_binding_start_up350 > 1) ? ($target_binding_start_up350 - 1) : $target_binding_start_up350;
			my $target_binding_end_up350_dangling = ($target_binding_end_up350 < $target_up350_down300_seq_len) ? ($target_binding_end_up350 + 1) : $target_binding_end_up350;
			
			my $srna_binding_seq_dangling = substr $sub_srna_seq, ($srna_binding_start_dangling - 1), ($srna_binding_end_dangling - $srna_binding_start_dangling + 1);
			my $target_binding_seq_dangling = substr $Target_seq_up350_down300, ($target_binding_start_up350_dangling - 1), ($target_binding_end_up350_dangling - $target_binding_start_up350_dangling + 1);
			
			# In some cases, seed regions expand beyond binding region.
			my $srna_seed_start_real = ($srna_seed_start >= $srna_binding_start_dangling) ? $srna_seed_start : $srna_binding_start_dangling;
			my $srna_seed_end_real = ($srna_seed_end <= $srna_binding_end_dangling) ? $srna_seed_end : $srna_binding_end_dangling;
			my $seed_len_real = $srna_seed_end_real - $srna_seed_start_real + 1;
			my $srna_seed_restriction = "<" x ($srna_seed_start_real - $srna_binding_start_dangling) . "(" x $seed_len_real . "<" x ($srna_binding_end_dangling - $srna_seed_end_real); 
			
			# In some cases, seed regions expand beyond binding region.
			my $target_seed_start_up350_real = ($target_seed_start_up350 >= $target_binding_start_up350_dangling) ? $target_seed_start_up350 : $target_binding_start_up350_dangling;
			my $target_seed_end_up350_real = ($target_seed_end_up350 <= $target_binding_end_up350_dangling) ? $target_seed_end_up350 : $target_binding_end_up350_dangling;
			my $target_seed_restriction = ">" x ($target_seed_start_up350_real - $target_binding_start_up350_dangling) . ")" x $seed_len_real . ">" x ($target_binding_end_up350_dangling - $target_seed_end_up350); 
			
			# my $temp_file_binding = 'temp_for_RNAcofold_binding.txt'; #commented on April 9, 2010
			my $seq_binding = $srna_binding_seq_dangling . '&' . $target_binding_seq_dangling;
			my $restriction_binding = $srna_seed_restriction . '&' . $target_seed_restriction;
			#print "$seq_binding\n$restriction_binding\n";##
			my @rnacofold_lines = ` echo "$seq_binding\n$restriction_binding\n"|RNAcofold -C -noLP -noPS `;
			chomp $rnacofold_lines[1];
			my ($structure_dangling) = split /\s/, $rnacofold_lines[1];
			#print @rnacofold_lines;##
			#die "Finished prior to pairwise alignment\n";##
			# pairwise alignment
			my $srna_seed_start_in_binding = $srna_seed_start_real - $srna_binding_start_dangling + 1;
			my $target_seed_start_in_binding = $target_seed_start_up350_real - $target_binding_start_up350_dangling + 1;
			my $pairwise_array_ref = &pairwise_alignment($srna_binding_seq_dangling, $target_binding_seq_dangling, $structure_dangling, $srna_seed_start_in_binding, $target_seed_start_in_binding, $seed_len_real);
			my @pairwise_array = @$pairwise_array_ref;
			my $pairwise_srna_seq = $pairwise_array[0];
			my $pairwise_bar = $pairwise_array[1];
			my $pairwise_target_seq = $pairwise_array[2];
			
			# output
			# Convert the target position to start codon
			my $target_binding_start_AUG = ($target_binding_start_up350 > $target_name_upstream_up350{$Target_Name}) ? ($target_binding_start_up350 - $target_name_upstream_up350{$Target_Name}) : ($target_binding_start_up350 - $target_name_upstream_up350{$Target_Name} - 1);
			my $target_binding_end_AUG = ($target_binding_end_up350 > $target_name_upstream_up350{$Target_Name}) ? ($target_binding_end_up350 - $target_name_upstream_up350{$Target_Name}) : ($target_binding_end_up350 - $target_name_upstream_up350{$Target_Name} - 1);
			my $target_binding_site_AUG = $target_binding_start_AUG . ',' . $target_binding_end_AUG;
			my $sRNA_seed_position_new = $srna_seed_start_real . ',' . $srna_seed_end_real;
			my $target_seed_start_real_AUG = ($target_seed_start_up350_real > $target_name_upstream_up350{$Target_Name}) ? ($target_seed_start_up350_real - $target_name_upstream_up350{$Target_Name}) : ($target_seed_start_up350_real - $target_name_upstream_up350{$Target_Name} - 1);
			my $target_seed_end_real_AUG = ($target_seed_end_up350_real > $target_name_upstream_up350{$Target_Name}) ? ($target_seed_end_up350_real - $target_name_upstream_up350{$Target_Name}) : ($target_seed_end_up350_real - $target_name_upstream_up350{$Target_Name} - 1);
			my $target_seed_position_new = $target_seed_start_real_AUG . ',' . $target_seed_end_real_AUG;
			
			# Get gene product from hash table
			print OUT "\tProbability = $Probability\t$Energy\t$Hybrid_plus_open_both\t$Hybrid_seed\t$Hybrid_open_seed\n";  
			
			# align name in pairwise alignment
			my $alignment_srna_name = 'sRNA(' . $sub_srna_name . ') ';
			my @alignment_srna_name_array = split //, $alignment_srna_name;
			my $alignment_target_name = 'Target(' . $Target_Name . ') ';
			my @alignment_target_name_array = split //, $alignment_target_name;
			
			my $alignment_srna_name_array_len = scalar @alignment_srna_name_array;
			my $alignment_target_name_array_len = scalar @alignment_target_name_array;
			
			if ($alignment_srna_name_array_len <= $alignment_target_name_array_len) {
				for (my $i = 0; $i < ($alignment_target_name_array_len - $alignment_srna_name_array_len); $i++) {
					push @alignment_srna_name_array, ' ';
				}
			}
			else {
				for (my $i = 0; $i < ($alignment_srna_name_array_len - $alignment_target_name_array_len); $i++) {
					push @alignment_target_name_array, ' ';
				}
			}
			my $alignment_srna_name_aligned = join '', @alignment_srna_name_array;
			my $alignment_target_name_aligned = join '', @alignment_target_name_array;
			
			# align start in pairwise alignment
			my $alignment_srna_start =  $srna_binding_start_dangling;
			my @alignment_srna_start_array = split //, $alignment_srna_start;
			
			my $target_binding_end_AUG_dangling = ($target_binding_end_up350_dangling > $target_name_upstream_up350{$Target_Name}) ? ($target_binding_end_up350_dangling - $target_name_upstream_up350{$Target_Name}) : ($target_binding_end_up350_dangling - $target_name_upstream_up350{$Target_Name} - 1);
			my $alignment_target_end = $target_binding_end_AUG_dangling;
			my @alignment_target_end_array = split //, $alignment_target_end;
			
			my $alignment_srna_start_array_len = scalar @alignment_srna_start_array;
			my $alignment_target_end_array_len = scalar @alignment_target_end_array;
			
			if ($alignment_srna_start_array_len <= $alignment_target_end_array_len) {
				for (my $i = 0; $i < ($alignment_target_end_array_len - $alignment_srna_start_array_len); $i++) {
					splice @alignment_srna_start_array, 0, 0, ' ';
				}
			}
			else {
				for (my $i = 0; $i < ($alignment_srna_start_array_len - $alignment_target_end_array_len); $i++) {
					splice @alignment_target_end_array, 0, 0, ' ';
				}
			}
			my $alignment_srna_start_aligned = join '', @alignment_srna_start_array;
			my $alignment_target_end_aligned = join '', @alignment_target_end_array;
			
			# align end in pairwise alignment
			my $alignment_srna_end =  $srna_binding_end_dangling;
			my @alignment_srna_end_array = split //, $alignment_srna_end;
			
			my $target_binding_start_AUG_dangling = ($target_binding_start_up350_dangling > $target_name_upstream_up350{$Target_Name}) ? ($target_binding_start_up350_dangling - $target_name_upstream_up350{$Target_Name}) : ($target_binding_start_up350_dangling - $target_name_upstream_up350{$Target_Name} - 1);
			my $alignment_target_start = $target_binding_start_AUG_dangling;
			my @alignment_target_start_array = split //, $alignment_target_start;
			
			my $alignment_srna_end_array_len = scalar @alignment_srna_end_array;
			my $alignment_target_start_array_len = scalar @alignment_target_start_array;
			
			if ($alignment_srna_end_array_len <= $alignment_target_start_array_len) {
				for (my $i = 0; $i < ($alignment_target_start_array_len - $alignment_srna_end_array_len); $i++) {
					splice @alignment_srna_end_array, 0, 0, ' ';
				}
			}
			else {
				for (my $i = 0; $i < ($alignment_srna_end_array_len - $alignment_target_start_array_len); $i++) {
					splice @alignment_target_start_array, 0, 0, ' ';
				}
			}
			my $alignment_srna_end_aligned = join '', @alignment_srna_end_array;
			my $alignment_target_start_aligned = join '', @alignment_target_start_array;
				
			my $total_len = scalar @alignment_srna_name_array + scalar @alignment_srna_start_array + 1;
			my $blank_str = '';
			for (my $i = 0; $i < $total_len; $i++) {
				$blank_str .= ' ';
			}
			
			print OUT "\t$alignment_srna_name_aligned $alignment_srna_start_aligned $pairwise_srna_seq $alignment_srna_end_aligned\n";
			print OUT "\t$blank_str $pairwise_bar\n";
			print OUT "\t$alignment_target_name_aligned $alignment_target_end_aligned $pairwise_target_seq $alignment_target_start_aligned\n";
		}
	}
	unlink $temp_binding_site;
	unlink $temp_binding_site_sorted;
	close OUT;
	close OUTBIND;
	
	return 0;	
}

sub pairwise_alignment {
	my ($sub_srna_binding_seq, $sub_target_binding_seq, $sub_structure, $sub_srna_seed_start, $sub_target_seed_start, $sub_seed_len) = @_;

	my ($srna_notation, $target_notation) = split /&/, $sub_structure;
	$srna_notation .= '&';

	$sub_srna_binding_seq = lc $sub_srna_binding_seq;
	$sub_target_binding_seq = lc $sub_target_binding_seq;

	$sub_srna_binding_seq =~ s/t/u/g;
	$sub_target_binding_seq =~ s/t/u/g;

	my @srna_notation_array = split //, $srna_notation;
	my @srna_binding_seq_array = split //, $sub_srna_binding_seq;

	my @target_binding_seq_array = split //, $sub_target_binding_seq;

	#uppercase seed region
	for (my $i = ($sub_srna_seed_start - 1); $i < ($sub_srna_seed_start + $sub_seed_len - 1); $i++) {
		$srna_binding_seq_array[$i] = uc $srna_binding_seq_array[$i];
	}

	for (my $i = ($sub_target_seed_start - 1); $i < ($sub_target_seed_start + $sub_seed_len - 1); $i++) {
		$target_binding_seq_array[$i] = uc $target_binding_seq_array[$i];
	}

	my $target_notation_rev = reverse $target_notation;
	my @target_notation_rev_array = split //, $target_notation_rev;
	my @target_binding_seq_rev_array = reverse @target_binding_seq_array;

	my $i = 0;
	while (1) {
		if ($srna_notation_array[$i] eq '&') {
			last;
		}
		else {
			if ($srna_notation_array[$i] eq '(') {
				if ($target_notation_rev_array[$i] eq '.') {
					splice @srna_notation_array, $i, 0, '-';
					splice @srna_binding_seq_array, $i, 0, '-';
				}
			}
			elsif ($srna_notation_array[$i] eq '.') {
				if ($i == scalar @target_notation_rev_array) {
					push @target_notation_rev_array, '-';
					push @target_binding_seq_rev_array, '-';
				}
				else {
					if ($target_notation_rev_array[$i] eq ')') {
						splice @target_notation_rev_array, $i, 0, '-';
						splice @target_binding_seq_rev_array, $i, 0, '-';
					}
				}
			}
			$i++;
		}
	}

	my $pairwise = '';
	my $len_extend = scalar @srna_notation_array;
	for ($i = 0; $i < $len_extend; $i++) {
		if (($srna_notation_array[$i] eq '(') and ($target_notation_rev_array[$i] eq ')')) {
			if (($srna_binding_seq_array[$i] eq 'a') or ($srna_binding_seq_array[$i] eq 'A') or ($srna_binding_seq_array[$i] eq 'c') or ($srna_binding_seq_array[$i] eq 'C')) {
				$pairwise .= '|';
			}
			elsif (($srna_binding_seq_array[$i] eq 'g') or ($srna_binding_seq_array[$i] eq 'G')) {
				if (($target_binding_seq_rev_array[$i] eq 'u') or ($target_binding_seq_rev_array[$i] eq 'U')) {
					$pairwise .= '.';
				}
				else {
					$pairwise .= '|';
				}
			}
			else {  # $srna_binding_seq_array[$i] eq 'u'
				if (($target_binding_seq_rev_array[$i] eq 'g') or ($target_binding_seq_rev_array[$i] eq 'G')) {
					$pairwise .= '.';
				}
				else {
					$pairwise .= '|';
				}
			}
		}
		else {
			$pairwise .= ' ';
		}
	}

	my $srna_binding_seq_new = join '', @srna_binding_seq_array;
	my $target_binding_seq_rev_new = join '', @target_binding_seq_rev_array;
	
	my @array = ();
	push @array, $srna_binding_seq_new;
	push @array, $pairwise;
	push @array, $target_binding_seq_rev_new;
	
	return \@array;
}

sub energy_open_binding_site {
	my ($seq, $binding_start, $binding_end) = @_;
	my (@rnafold_results, $energy_start, $energy_end, $energy_structured, $energy_unstructured);
	my ($len, $restriction, $temp_file);
	# calculate the structured energy
	@rnafold_results = `echo $seq|RNAfold -p0 -noPS`;
	$energy_structured = substr($rnafold_results[2],index($rnafold_results[2],"-"),index($rnafold_results[2],"k") - index($rnafold_results[2],"-") - 1);
	# calculate the unstructured energy
	$len = length $seq;
	$restriction = "." x ($binding_start - 1) . "x" x ($binding_end - $binding_start + 1) . "." x ($len - $binding_end); 
	@rnafold_results = ` echo "$seq\n$restriction\n"|RNAfold -p0 -C -noPS `;
	$energy_unstructured = substr($rnafold_results[2],index($rnafold_results[2],"-"),index($rnafold_results[2],"k") - index($rnafold_results[2],"-") - 1);
	my $open_binding_site_energy = $energy_unstructured - $energy_structured;
	
	return $open_binding_site_energy;
}

sub unstructured_energy_binding_site {
	my ($seq, $binding_start, $binding_end) = @_;
	my (@rnafold_results, $energy_start, $energy_end, $energy_structured, $energy_unstructured);
	my ($len, $restriction, $temp_file);
	# calculate the unstructured energy
	$len = length $seq;
	$restriction = "." x ($binding_start - 1) . "x" x ($binding_end - $binding_start + 1) . "." x ($len - $binding_end); 
	@rnafold_results = ` echo "$seq\n$restriction\n"|RNAfold -p0 -C -noPS `;
	$energy_unstructured = substr($rnafold_results[2],index($rnafold_results[2],"-"),index($rnafold_results[2],"k") - index($rnafold_results[2],"-") - 1);
	
	return $energy_unstructured;
}	

sub energy_open_seed {
	my ($seq, $binding_start, $binding_end) = @_;
	my (@rnafold_results, $energy_start, $energy_end, $energy_structured, $energy_unstructured);
	my ($len, $restriction, $temp_file);
	# calculate the structured energy
	@rnafold_results = `echo $seq|RNAfold -noPS`;
	$energy_start = rindex $rnafold_results[1], '(';
	$energy_end = rindex $rnafold_results[1], ')';
	$energy_structured = substr $rnafold_results[1], ($energy_start + 1), ($energy_end - 1 - $energy_start);
	$energy_structured =~ s/\s//g;
	# calculate the unstructured energy
	$len = length $seq;
	$restriction = "." x ($binding_start - 1) . "x" x ($binding_end - $binding_start + 1) . "." x ($len - $binding_end); 
	@rnafold_results = ` echo "$seq\n$restriction\n"|RNAfold -C -noPS `;
	$energy_start = rindex $rnafold_results[1], '(';
	$energy_end = rindex $rnafold_results[1], ')';
	$energy_unstructured = substr $rnafold_results[1], ($energy_start + 1), ($energy_end - 1 - $energy_start);
	$energy_unstructured =~ s/\s//g;
	my $open_binding_site_energy = $energy_unstructured - $energy_structured;
	
	return $open_binding_site_energy;
}

sub unstructured_energy_seed {
	my ($seq, $binding_start, $binding_end) = @_;
	my (@rnafold_results, $energy_start, $energy_end, $energy_structured, $energy_unstructured);
	my ($len, $restriction, $temp_file);
	# calculate the unstructured energy
	$len = length $seq;
	$restriction = "." x ($binding_start - 1) . "x" x ($binding_end - $binding_start + 1) . "." x ($len - $binding_end); 
	@rnafold_results = ` echo "$seq\n$restriction\n"|RNAfold -C -noPS `; 
	$energy_start = rindex $rnafold_results[1], '(';
	$energy_end = rindex $rnafold_results[1], ')';
	$energy_unstructured = substr $rnafold_results[1], ($energy_start + 1), ($energy_end - 1 - $energy_start);
	$energy_unstructured =~ s/\s//g;
	
	return $energy_unstructured;
}	

sub in_array {
   my ($array_ref, $element) = @_;
   my $flag = 0;
   my @array = @$array_ref;
   print "in_array @array\n" if !defined $element;
   foreach (@array) {
      if ($element eq $_) {
         $flag = 1;
         last;
      }
   }
   return $flag;
}


sub basecontent {
	my ($seqstr) = @_;
	$seqstr =~ s/U/T/g;
  $seqstr =~ s/u/t/g;
  my @features = ();
  my $single = 0;
  
  my @a = $seqstr =~ /a/gi;
  my @c = $seqstr =~ /c/gi;
  my @g = $seqstr =~ /g/gi;
  my @t = $seqstr =~ /t/gi;
  my $a = @a;
  my $c = @c;
  my $g = @g;
  my $t = @t;
  $single = $a + $c + $g +$t;
  
	$a = ($a / $single) * 100;
	$t = ($t / $single) * 100;
	push @features, $a;
	push @features, $t;
	
	my $features_ref = \@features;
	return $features_ref;
}	
