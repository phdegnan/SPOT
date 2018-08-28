#!/usr/bin/perl

  #use strict;
  #use warnings;
  use SOAP::Lite;
  use HTTP::Cookies;
  use List::Util qw(shuffle);
  my $hintaRNAout = $ARGV[0];
  my $pvalcutoff = $ARGV[1];

  if(not $pvalcutoff) {
      $pvalcutoff = 0.01;
  }  

  my @genelist = ();
  my @background = ();
  my $GeneID = "";

  open(MYDATA, $hintaRNAout) or die("Error: cannot open file" . $hintaRNAout . "\n");
  my @hintaRNAoutlines = <MYDATA>;
  close MYDATA;
  my $pid =0; ##added PHD
  foreach(@hintaRNAoutlines) {
      my @split = split(/,/, $_);
      my $pval = $split[0];
      if($split[1]) { # check if entry is not blank for org. of interest
          if ($split[1] =~ m/GeneID:(\d+)\)/) {
              $GeneID = $1;
          #}elsif($split[1] =~ m/GI:(\d+)\)/) {##PD
          #   $GeneID = $1;
	  #    $pid++;
		if(length($GeneID) >= 8 ){$pid++;} ##added PHD
          }
          if($pval <= $pvalcutoff) {
               push(@genelist, $GeneID);
          }
          push(@background, $GeneID);
      }  
  }
  
  shift(@genelist);
  shift(@background); 
  #@RANDOM=shuffle @background;
  #@background=@RANDOM[0..399];
  my $lengthgenelist = scalar(@genelist);
  my $lengthback = scalar(@background);
 
  print "length gene list: $lengthgenelist\n";
  
  print "length background: $lengthback\n";
  
  my $inputIds = join(',', @genelist);
  my $inputIdsBack = join(',', @background);
#print "\n$inputIds\n";
#print "\n$inputIdsBack\n";

  my $soap = SOAP::Lite                             
     -> uri('http://service.session.sample')                
     -> proxy('http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService',
                cookie_jar => HTTP::Cookies->new(ignore_discard=>1));

 #user authentication by email address
 #For new user registration, go to http://david.abcc.ncifcrf.gov/webservice/register.htm 
 my $check = $soap->authenticate('pdegnan@illinois.edu')->result; # set your email here
  	print "\nUser authentication: $check\n";

 if (lc($check) eq "true") { 	
 #addList
 #WARNING: user should limit the number of genes in the list within 3000, or DAVID will turn down the service due to resource limitation.
 my $idType = "";
 if(($pid/$lengthgenelist ) >= 0.90){	##PD
 	$idType = 'PROTEIN_GI_ACCESSION';
 }else{
 	$idType = 'ENTREZ_GENE_ID'; 	##PD
 }
 my $listName = 'make_up';
 my $listType=0;
 my $list = $soap ->addList($inputIds, $idType, $listName, $listType)->result;
  	#print "\n[$idType][$inputIds]\n";
	print "\n$list of list was mapped\n"; 


 my $listNameBack = 'back';
 my $listTypeBack=1;
 my $listBack = $soap ->addList($inputIdsBack, $idType, $listNameBack, $listTypeBack)->result;
        #print "\n[$idType][$inputIdsBack]\n";
	print "\n$listBack of list was mapped\n";
 
 #set user defined categories 
 #my $categories = $soap ->setCategories("abcd,BBID,BIOCARTA,COG_ONTOLOGY,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE")->result;
 #to user DAVID default categories, send empty string to setCategories():
 my $categories = $soap ->setCategories("BBID,BIOCARTA,GOTERM_BP_FAT,COG_ONTOLOGY,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE")->result;
 print "\nValid categories: \n$categories\n\n";  
 
open (termClusterReport, ">", "termClusterReport.txt");

my $overlap=3;
my $initialSeed = 3;
my $finalSeed = 3;
my $linkage = 0.5;
my $kappa = 50;
my $termClusterReport = $soap->getTermClusterReport($overlap,$initialSeed,$finalSeed,$linkage,$kappa);

my %hash =  %{$termClusterReport->result}; # pull first result
#foreach $keys (keys(%hash)){
#	print "[$keys][$hash{$keys}]\n";
#}
#}
#__END__

my $switch = 0;

    $swtich = 1;
	my @simpleTermClusterRecords = $termClusterReport->paramsout; 	## array elements = simpleTermClusterRecord=HASH(0x2238258) only 5 things...
	print "Total TermClusterRecords: ".(@simpleTermClusterRecords+1)."\n\n"; 
	my @simpleTermClusterRecordKeys = keys %{$termClusterReport->result};	## array elements = [score][name][simpleChartRecords]	
	my @simpleTermClusterRecordValues = values %{$termClusterReport->result}; ## [0.6269632674571864[GO:0005976~polysaccharide metabolic process][ARRAY(0x2243ec8)]	
	@chartRecords = @{$hash{"simpleChartRecords"}};
	
	print termClusterReport "Annotation Cluster 1\tEnrichment Score:  $hash{'score'}\n";
	print termClusterReport "Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n";
	for $j (0 .. (@chartRecords-1))
	{			
		%chartRecord = %{$chartRecords[$j]};	
			
		my $categoryName = $chartRecord{"categoryName"};
		my $termName = $chartRecord{"termName"};
		my $listHits = $chartRecord{"listHits"};
		my $percent = $chartRecord{"percent"};
		my $ease = $chartRecord{"ease"};
		my $Genes = $chartRecord{"geneIds"};
		my $listTotals = $chartRecord{"listTotals"};
		my $popHits = $chartRecord{"popHits"};
		my $popTotals = $chartRecord{"popTotals"};
		my $foldEnrichment = $chartRecord{"foldEnrichment"};
		my $bonferroni = $chartRecord{"bonferroni"};
		my $benjamini = $chartRecord{"benjamini"};
		my $FDR = $chartRecord{"afdr"};
		
		print termClusterReport "$categoryName\t$termName\t$listHits\t$percent\t$ease\t$Genes\t$listTotals\t$popHits\t$popTotals\t$foldEnrichment\t$bonferroni\t$benjamini\t$FDR\n";
	
	}	
	for $k (0 .. (@simpleTermClusterRecords-1))
	{	
		my $itr=$k+2;
		@simpleTermClusterRecordValues = values %{$simpleTermClusterRecords[$k]};
		%current=%{$simpleTermClusterRecords[$k]};
		@chartRecords = @{$current{'simpleChartRecords'}};
		
		print termClusterReport "\nAnnotation Cluster $itr\tEnrichment Score:  $current{'score'}\n";
		print termClusterReport "Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n";
		for $j (0 .. (@chartRecords-1))
		{			
			%chartRecord = %{$chartRecords[$j]};
			
			my $categoryName = $chartRecord{"categoryName"};
			my $termName = $chartRecord{"termName"};
			my $listHits = $chartRecord{"listHits"};
			my $percent = $chartRecord{"percent"};
			my $ease = $chartRecord{"ease"};
			my $Genes = $chartRecord{"geneIds"};
			my $listTotals = $chartRecord{"listTotals"};
			my $popHits = $chartRecord{"popHits"};
			my $popTotals = $chartRecord{"popTotals"};
			my $foldEnrichment = $chartRecord{"foldEnrichment"};
			my $bonferroni = $chartRecord{"bonferroni"};
			my $benjamini = $chartRecord{"benjamini"};
			my $FDR = $chartRecord{"afdr"};	
			
			print termClusterReport "$categoryName\t$termName\t$listHits\t$percent\t$ease\t$Genes\t$listTotals\t$popHits\t$popTotals\t$foldEnrichment\t$bonferroni\t$benjamini\t$FDR\n";		
		}
	}
	close termClusterReport;
	print "termClusterReport.txt generated\n";
}
