#!/usr/bin/perl

  #use strict;
  #use warnings;
  use SOAP::Lite;
  use HTTP::Cookies;

  #call: /usr/local/perl/bin/perl termClusterReport.pl intarna.csv > david.out 

  my $intaRNAout = $ARGV[0];

  my @genelist = ();
  my @background = ();
  my $GeneID = "";
  my $top50 = 1; # enrichment done for top50

  open(MYDATA, $intaRNAout) or die("Error: cannot open file" . $intaRNAout . "\n");
      my @intaRNAoutlines = <MYDATA>;
  close MYDATA;
  my $pid =0; ##added PHD
  foreach(@intaRNAoutlines) {
      my @split = split(/;/, $_);
          #print "$split[1]\n";
          if ($split[19] =~ m/(^\d+)/ || $split[21] =~ m/(^\d+)/ { ##21??
              $GeneID = $1;
	      #print "[" . length($GeneID). "]\n";
	      if(length($GeneID) == 8 ){$pid++;} ##added PHD
          }
	  #$GeneID = $split[21];
          if($top50 <= 50) {
               #print "$GeneID\n";
               push(@genelist, $GeneID);
               $top50++;
          }
          push(@background, $GeneID);
  }
  
  my $lengthgenelist = scalar(@genelist);
  my $lengthback = scalar(@background);
 
  print "length gene list: $lengthgenelist\n";
  print "length background: $lengthback\n";
  
  my $inputIds = join(',', @genelist);
  my $inputIdsBack = join(',', @background);

#  print ">>>>$inputIds\n\n\n";
#  print "<<<<$inputIdsBack\n\n\n";
  
##############
############## at this point we have out gene list and background entrez gene IDs parsed from the hintaRNA.csv
##############

  my $soap = SOAP::Lite                             
     -> uri('http://service.session.sample')                
     -> proxy('http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService',
                cookie_jar => HTTP::Cookies->new(ignore_discard=>1));

 #user authentication by email address
 #For new user registration, go to http://david.abcc.ncifcrf.gov/webservice/register.htm 
 my $check = $soap->authenticate('patrick.degnan@ucr.edu')->result;
  	print "\nUser authentication: $check\n";

 if (lc($check) eq "true") { 	
 #addList
 #WARNING: user should limit the number of genes in the list within 3000, or DAVID will turn down the service due to resource limitation.
 my $idType = "";
 if(($pid/$lengthgenelist) > 0.90){
 	print "Using GI numbers\n"; ##PD
 	$idType = 'PROTEIN_GI_ACCESSION'; ##PD
 }else{
 	$idType = 'ENTREZ_GENE_ID';
 }
 #print "[$idType]\n";
 my $listName = 'make_up';
 my $listType=0;
 my $list = $soap ->addList($inputIds, $idType, $listName, $listType)->result;
  	print "\n$list of list was mapped\n"; 


 my $listNameBack = 'back';
 my $listTypeBack=1;
 my $listBack = $soap ->addList($inputIdsBack, $idType, $listNameBack, $listTypeBack)->result;
        print "\n$listBack of list was mapped\n";

 #set user defined categories 
 #my $categories = $soap ->setCategories("abcd,BBID,BIOCARTA,COG_ONTOLOGY,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE")->result;
 #to user DAVID default categories, send empty string to setCategories():
 my $categories = $soap ->setCategories("BBID,BIOCARTA,GOTERM_BP_FAT,COG_ONTOLOGY,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE")->result;
 print "\nValid categories: \n$categories\n\n";  

open (termClusterReport, ">", "termClusterReport.txt");

#test getTermClusterReport(int overlap,int initialSeed, int finalSeed, double linkage, int kappa)

my $overlap=3;
my $initialSeed = 2;
my $finalSeed = 2;
my $linkage = 0.5;
my $kappa = 85;
my $termClusterReport = $soap->getTermClusterReport($overlap,$initialSeed,$finalSeed,$linkage,$kappa);
#my $termClusterReport = $soap->getTermClusterReport();

my %hash =  %{$termClusterReport->result};

my $switch = 0;

#foreach my $key (keys %hash) {
#	print "[$key][$hash{$key}]\n";
#}
    $swtich = 1;


#if($switch) {
	my @simpleTermClusterRecords = $termClusterReport->paramsout; 	
	print "Total TermClusterRecords: ".(@simpleTermClusterRecords+1)."\n\n"; 
	my @simpleTermClusterRecordKeys = keys %{$termClusterReport->result};		
	my @simpleTermClusterRecordValues = values %{$termClusterReport->result};
		
	@chartRecords = @{$hash{'simpleChartRecords'}}; ##CHANGED PHD
	
	print termClusterReport "Annotation Cluster 1\tEnrichment Score:  $hash{'score'}\n"; ##CHANGED phd
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
		%current=%{$simpleTermClusterRecords[$k]};  ##added PHD
		@chartRecords = @{$current{'simpleChartRecords'}}; ##CHANGED phd
		
		print termClusterReport "\nAnnotation Cluster $itr\tEnrichment Score:  $current{'score'}\n"; ##CHANGED phd
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
#} else {

#print "termClusterReport.txt failed, no significant enrichment detected!\n";

#}
}
__END__

#push(@simpleTermClusterRecords,$termClusterReport->result);

