=head1 LICENCE

Copyright 2015 EMBL - European Bioinformatics Institute (EMBL-EBI)

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

=head1 NAME

import_iMits_RESTFULL.pl - iMITS/EMMA Database syncrhonisation

=head1 SYNOPSIS

=cut

use strict;
use LWP::UserAgent;
use DBI;
use JSON qw( decode_json );
use Try::Tiny;
use Getopt::Long;
use Log::Log4perl qw(:easy);
Log::Log4perl->init("log.conf");
my $log = Log::Log4perl->get_logger("iMitsImport");

## load configuration
eval { require "CFG.pl" } ||
  die "Configuration file 'CFG.pl' was not found.\n";

my $countStock = 0;
my %stockStrains = ();
my @stocks = ();
my %colony_names = ();
my %unknown_alleles = ();
my %endonuclease_mediated_mutations = ();

my %recessive_alleles_chromosome_location = (
					     'Tyrc-Brd' =>  7,
					     'Atm1Brd' => 2
					    );

my @production_centers = 
    ('GSF',
     'HMGU',
     'ICS',
     'Monterotondo',
     'MRC - Harwell',
     'Harwell',
     'WTSI',
     'CIPHE',
     'SEAT',
     'IMG',
     'INFRAFRONTIER-Oulu',
     'INFRAFRONTIER-VETMEDUNI');

my %tm_attempts = ( 1 => 'first', 2 => 'second', 3 => 'third', 4 => 'fourth', 5 => 'fifth', 6 => 'sixth', 7 => 'seventh');
my %ilar_codes = {};

my ($phenotyping_arg, $colony_arg, $start_page, $last_page, $per_page, $filter, $debug, $help) = (0, undef, 1, undef, 100, '', 0, 0);
my $opt = GetOptions('phenotyping', \$phenotyping_arg, 
		     'debug', \$debug, 
		     'colony:s', \$colony_arg,
		     'start_page:s', \$start_page,
		     'last_page:s', \$last_page,
		     'per_page:s', \$per_page,
		     'filter:s', \$filter,
		     'help', \$help

		    );

if ($help) {

  print STDERR  "Usage: $0 --phenotyping=<Phenotyping attempts> --colony=<iMits colony name>\n";
  print STDERR "--phenotyping       Flag: parse phenotyping attempts and not the micro injection attempts\n";
  print STDERR "--colony            Only retrieve and record the specific colony\n";
  #print STDERR "--production_center ";
  print STDERR "--first_page        What page to start downloading the data from the RESTful endpoint\n";
  print STDERR "--last_page         What page to stop downloading the data from the RESTful endpoint\n";
  print STDERR "--filter            Filter to add to the RESTful query\n";
  exit 1;

}

my $emma_dbh = DBI->connect($CFG::DSN, $CFG::USER, $CFG::PASSWD,
			    {
			     InactiveDestroy => 1, RaiseError => 1, PrintError => 1}
			   ) or die "Unable to connect: $DBI::errstr\n";

my $komp2_dbh;

build_ilar_codes();
processAttempts();

$emma_dbh->disconnect();
$komp2_dbh->disconnect();

exit 0;

sub callRESTFUL {
  my ($url) = @_;
  my $success = 0;
  my $json_result = undef;
  my $nb_tries = 0;

  while ($success == 0 && $nb_tries < 5) {
	
    $nb_tries++;
	
    print $url . "\n";
        
    my $curl_command = "curl --globoff --basic --user $CFG::IMITSUSER:$CFG::IMITSPASSWD -X GET '$url' --silent";
    $json_result =`$curl_command`;
	
    #print $json_result . "\n";
	
    if ($json_result =~ m/Please try again later/) {
      print STDERR "Error Fetching data from iMITS.\n";
    } else {
      $success = 1;
    }
  }
    
  my $decode_properly = 1;
  my $array_decoded = undef;
    
  try {
	
    $array_decoded = decode_json($json_result);
	
  } catch {
    $decode_properly == 0;
    print "ERROR: Fetching Data from iMITS end point for URL $url: $_\n";
    exit 1;
  };
  return $array_decoded;
    
}

sub getParentStrainData {
  my ($colony_arg) = @_;

  $colony_arg =~ s/\s/%20/g;
  my $url = "https://www.mousephenotype.org/imits/mi_attempts.json?colony_name_eq=$colony_arg";

  my $array_decoded = &callRESTFUL($url);
  my $decoded = $array_decoded->[0];
  my $pipeline = $decoded->{'pipeline_name'};
  my $consortium = $decoded->{'consortium_name'};
  my $epd_id = $decoded->{'es_cell_name'}; 
  my $es_cell_info = &getESCellData($epd_id);
  my @copy = @{$es_cell_info};
  push(@copy, $epd_id);
  push(@copy, $pipeline);
  push(@copy, $consortium);
    
  return \@copy;
} 

sub getESCellData {
  my ($es_cell_name) = @_;

  print "Get ES Cell data for $es_cell_name\n";
    
  my $array_decoded = &callRESTFUL("http://www.i-dcc.org/imits/targ_rep/es_cells.json?name_eq=$es_cell_name");
    
  my $decoded = $array_decoded->[0];

  my @result = ();
  push (@result, $decoded->{'ikmc_project_id'});
  push (@result, $decoded->{'strain'});
  push (@result, $decoded->{'parental_cell_line'});

  return \@result;

}

sub processAttempts {


  my $total_entries = 0;

  my $current_page = $start_page;

  while (1) {

    my $url;
    if ($phenotyping_arg) {

      if (defined($colony_arg)) {
		
	my $colony_id = $colony_arg;
	$colony_id =~ s/\s/%20/g;
	$url = "https://www.mousephenotype.org/imits/phenotype_attempts.json?colony_name_eq=$colony_id";
			
      } else {
	    
	$url = "https://www.mousephenotype.org/imits/phenotype_attempts.json?page=$current_page&per_page=$per_page";
	if (length($filter) > 0) {
	  $filter =~ s/\s/%20/g;
	  $url .= "&" . $filter;
	}
      }
	    
    } else {
	
      if (defined($colony_arg)) {
		
	my $colony_id = $colony_arg;
	$colony_id =~ s/\s/%20/g;
	$url = "https://www.mousephenotype.org/imits/mi_attempts.json?colony_name_eq=$colony_id";

      } else {

	$url = "https://www.mousephenotype.org/imits/mi_attempts.json?page=$current_page&per_page=$per_page";
	if (length($filter) > 0) {
	  $filter =~ s/\s/%20/g;
	  $url .= "&" . $filter;
	}
      }
    }
	
    $log->info($url);
    my $array_decoded = &callRESTFUL($url);
	
    my $nb_entries = 0;

    if (defined($array_decoded)) {
	    
      my $count = 0;
      $nb_entries = scalar @{$array_decoded};
      foreach my $decoded (@{$array_decoded}) {
		
	my ($strain_name,$mutation_type,$marker_symbol,$allele,$sup_allele,$phenotype_allele,$cell_line_bg,$distribution_center,$pre_genetic_description,$genetic_description,$epd_id,$cell_line_bgsymbol,$mta_file,$str_status,$epd_id,$mi_date,$maintained_bg,$maintained_bgsymbol,$str_access,$phenotype_description,$full_allele_name,$allele_project,$repository,$source,$project,$theTime,$material_deposited,$colony_prefix,$mi_attempt_colony_prefix,$imits_status,$genotype_confirmed,$consortium,$freezing_started,$archived,$notes,$breeding,$archiving_method_id,$males,$females,$male_bg_id,$female_bg_id,$embryo_state,$maintenance,$mutant_viable,$mutant_fertile,$require_homozygous,$immunocompromised,$current_sanitary_status,$animal_husbandry,$project_id,$emma,$pipeline,$dist_ilar_code,$ilar_code,$escell_strain,$test_cross,$back_cross,$genotype_file,$tm_es_line,$tat_cre,$cre_deleter_strain,$is_cre_excision_required,$rederivation_complete,$mouse_allele_type,$is_active,$is_crispr_cas9,$impc_data_available,$mgi_accession_id);

	$count++;
	my $phenotype_attempts=$decoded->{'phenotype_attempts_count'};
	my $urlallele;
	$mgi_accession_id=$decoded->{"mgi_accession_id"};
	$log->debug("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n$phenotype_attempts\n$mgi_accession_id");

	$colony_prefix = $decoded->{'colony_name'};
	$log->info("Colony name:        $colony_prefix");

	# get the parent colony prefix

	if ($phenotyping_arg) {

	  $mi_attempt_colony_prefix = $decoded->{"mi_attempt_colony_name"};
	  $log->info("MI Attempt Colony:  $mi_attempt_colony_prefix");
	}

	$is_active = $decoded->{'is_active'};
	$repository = $decoded->{'production_centre_name'};

	if ($is_active eq 'false') {

	  $log->warn("Project not active anymore - Check EMMA distribution then skip");

	  if ($phenotyping_arg && $colony_prefix ne $mi_attempt_colony_prefix) {
	    my $rejected = &reject_strain_access_in_emma($colony_prefix);
	    $log->warn("REJECTED: $mi_attempt_colony_prefix => $colony_prefix") if ($rejected);
	  } else {
	    my $rejected = &reject_strain_access_in_emma($colony_prefix);
	    $log->warn("REJECTED: $colony_prefix") if ($rejected);
	  }

	  next;

	}

	if (!is_emma_repository($repository)) {
	  $log->warn("$repository is not an EMMA production node");
	  next;
	}

	$mi_date = $decoded->{'mi_date'};
	$imits_status = $decoded->{'status_name'};
	$genotype_confirmed = ($imits_status eq 'Genotype confirmed');
	$pipeline = $decoded->{'pipeline_name'};
	$consortium = $decoded->{'consortium_name'};

	$is_crispr_cas9 = ($decoded->{'mi_plan_mutagenesis_via_crispr_cas9'} eq 'true') ? 1 : 0;

	$marker_symbol = $decoded->{'marker_symbol'};
	$epd_id = $decoded->{'es_cell_name'};

	$log->info("Repository:         $repository");

	$test_cross = $decoded->{'test_cross_strain_name'};
	$back_cross = $decoded->{'colony_background_strain_name'};

	# For the time being we only report CRISPR CAS9 strains assuming the nomenclature will be correct

	if ($is_crispr_cas9 == 1 && $genotype_confirmed) {
		    
	  my $ilar = &get_ilar_code($repository);
	  $endonuclease_mediated_mutations{$colony_prefix} = $marker_symbol . "<em#" . $ilar . ">";
	  #call a subroutine to  grab colonies_attributes field from json response to get f1 colonies produced
	  #by this mi then iterate over ignoring genotype confirmed false results.
	  #&crispr_lines;
	  $pipeline = "CRISPR/CAS9"; #pipeline_name not suppied from json response for crisp/cas9 lines
	  my($colonies_attributes,$attributes_genotype_confirmed,$attributes_allele_name,
	     $attributes_mgi_allele_id,$attributes_name);
	  $colonies_attributes = $decoded->{"colonies_attributes"};
	  #Iterate over colonies ignoring genotype_confirmed=false
	  for (my $i = 0; $i < scalar(@{$colonies_attributes}); $i++) { 
	    $attributes_genotype_confirmed = $colonies_attributes->[$i]->{"genotype_confirmed"};
	    next if $attributes_genotype_confirmed ne 'true';
	    $attributes_allele_name = $colonies_attributes->[$i]->{"allele_name"};
	    $attributes_mgi_allele_id = $colonies_attributes->[$i]->{"mgi_allele_id"};
	    $attributes_name = $colonies_attributes->[$i]->{"name"}; 
	    print "Genotype Confirmed: " .$attributes_genotype_confirmed . "\n";
	    print "Allele name: " .$attributes_allele_name . "\n";
	    print "MGI Allele ID: " .$attributes_mgi_allele_id . "\n";
	    print "Name: " .$attributes_name  . "\n";
	  }
		    
	}


	if ($phenotyping_arg) {

	  $mouse_allele_type = $decoded->{"mouse_allele_type"};
	  $is_cre_excision_required = ($decoded->{"cre_excision_required"}) ? 1 : 0;
	  $impc_data_available = ($imits_status eq "Phenotyping Started" || $imits_status eq "Phenotyping Complete");
	  $consortium .= "/IMPC";

	  # For "a" allele type, we update the database to indicate that the strain is part of IMPC
	  # This is achieved by modifying the consortium information for this strain
	  # For other allele types, this is part of the emma_update procedure
	  print "PHENOTYPING IS SET TO $phenotyping_arg + IMPC DATA AVAILABILITY = $impc_data_available + MOUSE ALLELE TYPE = $mouse_allele_type + IMITS STATUS = $imits_status + GENOTYPE CONFIRMED = + $genotype_confirmed\n\n ";
	  if ($mouse_allele_type eq 'a') {
	    $log->info("ABOUT TO UPDATE CONSORTIUM INFO");
	    &update_strain_consortium_information($colony_prefix, $consortium);

				# If IMPC data are available, mark also this information in the database for this colony

	    if ($impc_data_available) {
	      $log->info("IMPC DATA IS AVAILABLE");
	      &set_impc_phenotype_data_exists($colony_prefix, $impc_data_available,$phenotype_attempts);
	    }
	  }

	}
		
	# Get the distribution center information
	# distribution_centres_attributes

	my $distribution_info = $decoded->{"distribution_centres_attributes"};

	$emma = 0;

	# Phenotyping occurs for tm1b, tm1c, tm1d, tm1e
	# Other use case: mouse_allele_type": ".1"
	if ($genotype_confirmed || ($phenotyping_arg && ($mouse_allele_type =~ /^[bcde]/ || $mouse_allele_type =~ /^\.\d{1}/))) {

	  # if there are distribution centers, then loop through and look for the EMMA network

	  if ($distribution_info && scalar(@{$distribution_info}) > 0) {

	    # this field is always empty
	    #print "Mouse allele symbol = " . $decoded->{"mouse_allele_symbol_superscript"} . "\n";
	    # Example "distribution_centres_attributes":                                                                                                                                                                                                      
	    # [{"distribution_network":"EMMA","end_date":null,"id":430,"is_distributed_by_emma":true,"start_date":null,                                                                                                                                       
	    # "deposited_material_name":"Live mice","centre_name":"WTSI","_destroy":false}]                                                                                                                                                                   
	    # A center can choose to distribute different type of material. So we keep all of them and concatenate them.

	    my $material_index = 0;
			
	    for (my $i = 0; $i < scalar(@{$distribution_info}); $i++) {
			    
	      if ($distribution_info->[$i]->{"distribution_network"} eq "EMMA") {
			$log->info("Deposited material: " . $distribution_info->[$i]->{"deposited_material_name"});
			#$allele = $decoded->{"mouse_allele_symbol_superscript"};
			$material_deposited = (($material_index == 0) ? "" : "/") .  $distribution_info->[$i]->{"deposited_material_name"};
			$distribution_center = $distribution_info->[$i]->{"centre_name"};
			$material_index++;
			$emma = 1;
	      }
	    }
	  }			# distribution info available
		    
	  #
	  # First check the strain exists in EMMA
	  #
    
	  if ($emma == 0) {
			
	    print "EMMA: not distributed by EMMA\n";
	    if ($phenotyping_arg && $colony_prefix ne $mi_attempt_colony_prefix) {
	      my $rejected = &reject_strain_access_in_emma($colony_prefix);
	      print "REJECTED: $mi_attempt_colony_prefix => $colony_prefix\n" if ($rejected);
	    } else {
	      my $rejected = &reject_strain_access_in_emma($colony_prefix);
	      print "REJECTED: $colony_prefix\n" if ($rejected);
	    }
	    next;
	  } else {
	    $log->error("PASSED: $colony_prefix");
	  }


	  # it could be that a line is not distributed anymore by a center
	  # in this case we need to check all distribution center
	  # example Biat                                                                                                                                                                                                                                                
	  # We only proceed if the genotype confirmed information is present
		    
	  if ($phenotyping_arg) {

	    # Update the information of the parent strain as IMPC related since
	    # the colony.
			
	    # get parent pipeline and epd_id
	    # this works for cre excision ONLY
	    my $parent_info = &getParentStrainData($mi_attempt_colony_prefix);

	    $project_id = $parent_info->[0];
	    $escell_strain = $parent_info->[1];
	    $tm_es_line = $parent_info->[2];		    
	    $epd_id = $parent_info->[3];
	    $pipeline = $parent_info->[4];
	    my $mi_attempt_consortium = $parent_info->[5];

	    if ($mi_attempt_colony_prefix ne $colony_prefix) {
			    
	      # indicate that it's an IMPC mice
	      &update_strain_consortium_information($mi_attempt_colony_prefix, $mi_attempt_consortium . "/IMPC");
			    
	    }

	    if ($imits_status ne 'Cre Excision Started') {
			    
	      # "status_dates":{"Phenotype Attempt Registered":"2012-06-15","Cre Excision Started":"2012-10-01","Cre Excision Complete":"2013-02-08"}                                                                                         
	      # store the cre excision complete date instead of the micro-injection date                                                                                                                                                      
			    
	      $mi_date = $decoded->{"status_dates"}->{"Cre Excision Complete"}; # e.g. 2013-02-08                                                                                                                                             
	      $genotype_confirmed = (defined($mi_date) && $mi_date =~ /^\d{4}\-\d{2}\-\d{2}$/) ? 1 : 0;
	    }
			
			
	    $tat_cre = $decoded->{"tat_cre"} eq 'true';
			
	    if (!$tat_cre) {
	      $cre_deleter_strain = $decoded->{"deleter_strain_name"}; #"MGI:3046308: Hprt<tm1(CMV-cre)Brd>"                                                                                                                                  
	      print "Deleter strain name: " . $cre_deleter_strain . "\n";
	    }
	  } else {

	    # ES Cell information
	    my $escell_info = &getESCellData($epd_id);
			
	    $project_id = $escell_info->[0];
	    $escell_strain = $escell_info->[1];
	    $tm_es_line = $escell_info->[2];
			
	  }

	  $mi_date = "$mi_date 00:00:00"; #hack as data format in new iMits BioMart changed                                                                                                                                                                              
		    
	  #######
	  # MICE SHENANIGANS
	  # BLACK LISTED STRAINS/
	  # MICRO-INJECTION FLIPPING DATES
	  # ALLELE NAME CORRECTION
	  #######

	  $mi_date = '2013-07-25 00:00:00' if ($colony_prefix eq 'MGJK');
	  $mi_date = '2012-12-18 00:00:00' if ($colony_prefix eq 'MFXF');
	  $mi_date = '2012-12-18 00:00:00' if ($colony_prefix eq 'MFXD');
	  $mi_date = '2010-09-15 00:00:00' if ($colony_prefix eq 'MDCK');
	  print "Date:               $mi_date\n";
		    
	  ### SPECIAL RULES - FOR SANGER (error in IMITS)
	  if ($marker_symbol eq 'Spata25' && $colony_prefix eq 'MGJR') {
	    $allele = 'Spata25<tm1(KOMP)Wtsi>';
	  }

	  print "Consortium:         $consortium\n";

	  # Choose correct allele information
	  if ($phenotyping_arg) {

	    # we have to use a different field in this case:
	    $allele = $decoded->{"mouse_allele_symbol"};
	    # extract the pipeline from the parent strain
	    #mi_attempt_colony_name": "GSF-HEPD0538_7_A07-1"

	  } else {

	    if (length($decoded->{'mouse_allele_symbol'}) > 0) {
	      $allele = $decoded->{'mouse_allele_symbol'};
	    } elsif (length($decoded->{'es_cell_allele_symbol'}) > 0) {
	      $allele = $decoded->{'es_cell_allele_symbol'};
	    }
	  }
		    
	  if ($allele) {
	    if ($allele =~ /EGFP_CreERT2/) {
			    
	      $allele =~ s/EGFP_CreERT2/EGFP\/Cre\/ERT2/; 
			    
	    }

	    $sup_allele = $allele;
	    $urlallele = $allele;			
	    $allele =~ s/<sup>/</;
	    $urlallele =~ s/<sup>//;
	    $urlallele =~ s/<\/sup>//;
	    $urlallele =~ s/^$marker_symbol//;
	    $log->debug("URL Allele3:      $urlallele");
	    $allele =~ s/<\/sup>/>/;

	    $log->debug("Allele:             $allele");
	    my $t_allele = 'targeted mutation';
	    my $p_allele = "Wellcome Trust Sanger Institute";
			
	    # Check the allele type                                                                                                                                                                                                            
	    # regular expression to capture the TM allele structure.                                                                                                                                                                               
	    if ($allele =~ /tm(\d+[bcde]{0,1}\.{0,1}\d*)/) {
	      $t_allele = 'targeted mutation ' . $1;
	    }
			
	    # check who produced                                                                                                                                                                                                                   
	    if ($allele =~ /Wtsi/) {
	      $p_allele = "Wellcome Trust Sanger Institute";
	    } elsif ($allele =~ /Hmgu/) {
	      $p_allele = "Helmholtz Zentrum Muenchen GmbH";
	    } elsif ($allele =~ /Vlcg/) {
	      $p_allele = "Velocigene";
	    } elsif ($allele =~ /Mbp/) {
	      $p_allele = "Mouse Biology Program, UCDavis";
	    } else {
	      $log->error("ERROR: Production Center - The allele nomenclature $allele is not valid for this script - PLEASE FIX THIS - and restart the script");
	      # go to next entry
	      next;
	    }
			
	    # e.g. targeted mutation 1b, Wellcome Trust Sanger Institute                                                                                                                                                                           
	    $allele_project = &check_pipeline_project($pipeline, $allele);
	    $full_allele_name = $allele_project . ' ' . $t_allele . ', ' . $p_allele;
	    $log->info("Full allele name:   " . $full_allele_name);		    
	    $log->info("Pipeline:           $pipeline");

	  } else {		# for GT not yet assigned
	    $allele = "$marker_symbol<not yet available>";
	    $full_allele_name = $allele_project . ' targeted mutation, Wellcome Trust Sanger Institute';
	  }




	  # BUSINESS RULE 1: only take EUCOMM lines or KOMP lines produced by HMGU, ICS, Monterontondo, Harwell, CIPHE, SEAT (CNRS Villejuif) or distributed by CNRS
	  # BUSINESS RULE 2: gene Wls with colony prefix gpr177 albino must be skipped
		    
	  # BUSINESS RULE 3: for Cre excised mice, only import Sanger lines for the time being
		    
	  my $is_emma_compliant = 
	    &emma_compliant($marker_symbol, 
			    $colony_prefix, 
			    $pipeline, 
			    $repository, 
			    $distribution_center, 
			    $phenotyping_arg, 
			    $genotype_confirmed, 
			    $emma);

	  next if (!$is_emma_compliant);

	  # convert iMits centre codes to an EMMA standard one
	  $repository = &convert_center_information($repository);
	  # print "ALTERED DISTRIBUTION CENTRE IS NOW " .$distribution_center . "\n";
	  $distribution_center = &convert_center_information($distribution_center);
	  #print "ALTERED DISTRIBUTION CENTRE IS NOW " .$distribution_center . "\n";
	  $dist_ilar_code = &get_ilar_code($distribution_center);
	  $ilar_code = "/" . &get_ilar_code($repository);

	  ## CORRECT ES CELL, TESTCROSS and BACKCROSS INFORMATION
	  $back_cross = &correct_background_information($back_cross);

	  #
	  # STRAIN NOMENCLATURE RULES
	  #

	  if ($phenotyping_arg) {

	    if ($back_cross eq 'C57BL/6NTac') {
	      $strain_name = "C57BL/6N-$allele$ilar_code"; 
	    } elsif (
		     $back_cross eq 'C57BL/6Brd-Tyr;C57BL/6N' || 
		     $back_cross eq 'C57BL/6Brd-Tyr<c-Brd>;C57BL/6N;C57BL/6NTac' || 
		     $back_cross eq 'C57BL/6J-Tyr<c-Brd>;C57BL/6N' ||
		     $back_cross eq 'C57BL/6Brd-Tyr<c-Brd>;C57BL/6N'
		    ) {
	      #  we can consider "C57BL/6NTac" a subset of "C57BL/6N"
	      $strain_name = "B6Brd;B6N-Tyr<c-Brd> $allele$ilar_code";
	    } else {
	      $strain_name = "$back_cross-$allele$ilar_code";
	    }
			
	    print "STRAIN NAME:\t" . $strain_name . " ($back_cross)\n";
			
	    # make sure we don't write to the database
	    #next;

	  } else {

	    $escell_strain = &correct_background_information($escell_strain);
	    $test_cross = &correct_background_information($test_cross);
			
	    if ($debug) {
	      print "ES CELL:             $escell_strain\n";
	      print "TEST CROSS:          $test_cross\n";
	      print "BACK CROSS:          $back_cross\n";
	    }
			
			

	    if (!$escell_strain) {
	      $cell_line_bgsymbol = 'B6N';
	      $cell_line_bg = 'C57BL/6N'; 
	    } else {
	      $cell_line_bg = $back_cross; # ? why added this default value - leave for now
	      if ($escell_strain eq '129S7') {
		$cell_line_bgsymbol = '129S7';
		$cell_line_bg = '129S7/SvEvBrd-Hprt<b-m2>';
	      } elsif ($escell_strain eq '129P2/OlaHsd') {
		$cell_line_bgsymbol = '129P2';
	      } elsif ($escell_strain eq 'C57BL/6N') {
		$cell_line_bgsymbol = 'B6N';
	      } elsif ($escell_strain =~ 'A\<tm1Brd\>') {
		$cell_line_bgsymbol = 'B6N-A<tm1Brd>/a';
	      }
	    }

	    print "\n----------------\n";

	    print "$escell_strain\t$test_cross\t$back_cross\n";
	    #exit 1;

	    # decision July 1st 2010 - not to display maintained background for EUCOMM lines
	    $maintained_bg = '';
	    $maintained_bgsymbol = '';

	    # COMPUTE STRAIN NOMENCLATURE
	    # Calculate Strain Name from combination of ES cell, Test cross and Back cross backgrounds. 
	    # For all non-WTSI centres there are definite rules involving all 3 that are agreed with IMSR and the EMMA curators
	    # For WTSI there are agreed rules with them using the ES cell and test cross backgrounds only - in the near future Dave Melvin and Vivek will give the go ahead to switch to same rules as the other centres
	    # and their backcross data will be full populated
	    # Lines that do not fall into one of the agreed patters get assigned STOCK
			
	    $strain_name = 
	      &compute_strain_nomenclature($escell_strain, 
					   $test_cross, 
					   $back_cross, 
					   $marker_symbol, 
					   $allele, 
					   $repository, 
					   $ilar_code);
	  }

	  ### GENETIC DESCRIPTION

	  $genetic_description = "This mouse line originates from $allele_project ES clone $epd_id. For further details on the construction of this clone see the page at <a href=\"http://www.knockoutmouse.org/martsearch/ikmc_project/$project_id\">the IKMC portal</a>.";
	  # Removal of the targeting cassette using Flp recombinase is required to convert the targeted into a conditional allele - <a href=\"http://www.knockoutmouse.org/kb/entry/105/\">more information on conversion to the b,c and d allele forms</a>. Click <a href=http://www.knockoutmouse.org/about/targeting-strategies>here</a> for more information on $allele_project final vectors";

	  ## PHENOTYPIC DESCRIPTION
	  ## We remove all the old rules on the phenotype description linked to EUMODIC data/Sanger MGP
	  $phenotype_description = "Potential phenotyping data in the <a href=\"https://www.mousephenotype.org/data/search?q=$marker_symbol\">IMPC portal</a>";
	  if ($repository eq 'SANG') {
	    $maintenance = "Current information may be viewed at the <a href=\"http://www.sanger.ac.uk/mouseportal/search?query=$marker_symbol\">Sanger mouse portal</a>: Viability at weaning, Fertility, General Observations";
	  }
		    
	  # MTAs specification: 
	  # The MI date can influence what MTA to associate to the strain
	  $mta_file = 
	    &compute_mta_file($pipeline, 
			      $repository, 
			      $distribution_center, 
			      $consortium, 
			      $phenotyping_arg,
			      $allele,
			      $mi_date);

	  print STDERR "MTA: $mta_file\n";

	  # Default mutation type / mutation subtype: Targeted / Targeted Conditional
	  $mutation_type = 'TM/TC';

	  # Check allele nomenclature
	  if ($allele =~ /tm\d+a/) {
	    $mutation_type = 'TM/KOF';
	    print STDERR "$allele => " . $mutation_type . " ($marker_symbol)\n";
	  } elsif ($allele =~ /tm\d+e/) {
	    $mutation_type = 'TM/TNC';
	    $genetic_description = "This mouse line originates from $allele_project ES clone $epd_id. For further details on the construction of this clone see the page at <a href=\"http://www.mousephenotype.org/martsearch_ikmc_project/martsearch/ikmc_project/$project_id\">the IMPC portal</a>. The targeted allele has lost the 3' loxP site. These mutations cannot be converted into conditional alleles. Click <a href=http://www.knockoutmouse.org/about/targeting-strategies>here</a> for more information on $allele_project final vectors.";
	  } elsif ($allele =~/tm\d+b/) {
	    # tm1b: Reporter-tagged deletion allele (post-Cre)
	    $mutation_type = 'TM/RTD';
	    if ($repository eq 'SANG') {
	      if ($tat_cre) {
		$genetic_description = "This mouse line originates from $allele_project ES clone $epd_id. For further details on the construction of this clone see the page at <a href=\"http://www.mousephenotype.org/martsearch_ikmc_project/martsearch/ikmc_project/$project_id\">the IKMC page on the IMPC portal</a>. The critical exon(s) were flanked by loxP sites, and subsequent cre expression excised this critical sequence resulting in a knockout reporter allele using a cell permeable HTN-Cre as described in <a href=http://www.ncbi.nlm.nih.gov/pubmed/24197666>doi:10.1007/s11248-013-9764-x</a>. Click <a href=http://www.mousephenotype.org/martsearch_ikmc_project/about/targeting-strategies>here</a> for more information on $allele_project final vectors.";
	      } else {
				# uses an X-linked CMV-Cre driver for this purpose 
		$genetic_description = "This mouse line originates from $allele_project ES clone $epd_id. For further details on the construction of this clone see the page at <a href=\"http://www.mousephenotype.org/martsearch_ikmc_project/martsearch/ikmc_project/$project_id\">the IKMC page on the IMPC portal</a>. The critical exon(s) were flanked by loxP sites, and subsequent cre expression excised this critical sequence resulting in a knockout reporter allele using an X-linked CMV-Cre driver (EM:05414; C57BL/6N-Hprt<sup>tm1(CMV-cre)Wtsi</sup>/WtsiOulu) for this purpose. Click <a href=http://www.mousephenotype.org/martsearch_ikmc_project/about/targeting-strategies>here</a> for more information on $allele_project final vectors.";
	      }
	    } else {
	      if (!$tat_cre && defined($cre_deleter_strain)) {
				# cre deleter strain in HTML
		my $html_deleter = $cre_deleter_strain;
		$html_deleter =~ s/</<sup>/g;
		$html_deleter =~ s/>/<\/sup>/g;

		$genetic_description = "This mouse line originates from $allele_project ES clone $epd_id. For further details on the construction of this clone see the page at <a href=\"http://www.mousephenotype.org/martsearch_ikmc_project/martsearch/ikmc_project/$project_id\">the IKMC page on the IMPC portal</a>. The critical exon(s) were flanked by loxP sites, and subsequent cre expression excised this critical sequence resulting in a knockout reporter allele ($cre_deleter_strain). Click <a href=http://www.mousephenotype.org/martsearch_ikmc_project/about/targeting-strategies>here</a> for more information on $allele_project final vectors.";
	      } else {
		$genetic_description = "This mouse line originates from $allele_project ES clone $epd_id. For further details on the construction of this clone see the page at <a href=\"http://www.mousephenotype.org/martsearch_ikmc_project/martsearch/ikmc_project/$project_id\">the IKMC page on the IMPC portal</a>. The critical exon(s) were flanked by loxP sites, and subsequent cre expression excised this critical sequence resulting in a knockout reporter allele. Click <a href=http://www.mousephenotype.org/martsearch_ikmc_project/about/targeting-strategies>here</a> for more information on $allele_project final vectors.";
	      }
	    }
	    print "Based on allele infotmation: " . $mutation_type . "\t" . $genetic_description . "\n";
	  } elsif ($allele =~/tm(\d+)(c)/) {
	    # tm1c: Conditional allele (post-Flp)
	    $mutation_type = 'TM/TC';
	    my $tm_att = $1;
	    my $tm_c = "tm$1$2"; 
	    $genetic_description = "This line originates from $allele_project ES clone $epd_id, after breeding with a Flp recombinase deleter line to convert the original targeted allele tm".$tm_att."a (knock-out first allele) into a conditional allele " .$tm_c . ". For further details on the construction of this clone see the page at <a href=\"http://www.mousephenotype.org/martsearch_ikmc_project/martsearch/ikmc_project/$project_id\">the IKMC page on the IMPC portal</a>.";
	    # Removal of the targeting cassette using Flp recombinase is required to convert the targeted into a conditional allele - more information on conversion to the b,c and d allele forms. Click <a href=http://www.mousephenotype.org/martsearch_ikmc_project/about/targeting-strategies>here</a> for more information on $allele_project final vectors.";
	    print "Based on allele information: " . $mutation_type . "\t" . $genetic_description . "\n";

	  } elsif ($allele =~/tm(\d+)(d)/) {

	    # tm1d: Deletion allele (post-Flp and Cre with no reporter)
	    $mutation_type = 'TM/TC';
	    my $tm_att = $1;
	    my $tm_d = "tm$1$2";
	    $genetic_description = "This line originates from $allele_project ES clone $epd_id, after breeding with a Flp recombinase deleter line to convert the original targeted allele tm".$tm_att."a (knock-out first allele) into a deletion allele " .$tm_d . " (post-Flp and Cre with no reporter). For further details on the construction of this clone see the page at <a href=\"http://www.mousephenotype.org/martsearch_ikmc_project/martsearch/ikmc_project/$project_id\">the IKMC page on the IMPC portal</a>.";
	    #Removal of the targeting cassette using Flp recombinase is required to convert the targeted into a conditional allele - more information on conversion to the b,c and d allele forms. Click <a href=http://www.mousephenotype.org/martsearch_ikmc_project/about/targeting-strategies>here</a> for more information on $allele_project final vectors.";
	    print "Based on allele information: " . $mutation_type . "\t" . $genetic_description . "\n";

	  } elsif ($allele =~/EGFP\/Cre\/ERT2/) {
			
	    print "EUCOMMTOOLS Cre allele: $allele $repository\n";
	    $mutation_type = 'TM/KI';
	    if ($repository eq 'SANG') {
	      #print "Genetic description PRE: $genetic_description\n";
	      $phenotype_description="";
	      $pre_genetic_description = "Quality Control/Disclaimer: The structure of the targeted mutation is not fully verified. It is recommended that the recipient confirms the mutation structure by Southern blotting. The Cre expression pattern is expected to mirror that of the host gene, but this may not be fully described and/or may be modified in the targeted allele. Although characterisation of the Cre expression pattern is an objective of EUCOMMTools this is an ongoing process which may be incomplete. It is therefore recommended that the recipient verifies the Cre expression pattern. For the $allele allele, the promoter-driven selection cassette is flanked with Rox sites and may be removed using Dre recombinase.";
	      $genetic_description = $genetic_description . "<br/>\n" . $pre_genetic_description;
	      print "Genetic description: $genetic_description\n";
	    } elsif ($repository eq 'CNR') {
	      $pre_genetic_description = "Quality Control/Disclaimer: The structure of the targeted mutation is not fully verified. It is recommended that the recipient confirms the mutation structure by Southern blotting. The Cre expression pattern is expected to mirror that of the host gene, but this may not be fully described and/or may be modified in the targeted allele. Although characterisation of the Cre expression pattern is an objective of EUCOMMTools this is an ongoing process which may be incomplete. It is therefore recommended that the recipient verifies the Cre expression pattern. For the $allele allele, the promoter-driven selection cassette is flanked with Rox sites and may be removed using Dre recombinase.";
	      $genetic_description = $genetic_description . "<br/>\n" . $pre_genetic_description;
	      print "Genetic description: $genetic_description\n";
	      print "MUTATION TYPE: $mutation_type\n";

	    } elsif ($repository eq 'ICS') {
	      $genetic_description = '';
	    }
			
	    # Alternative text when characterisation done
	    # Cre-ER<sup>T2</sup> fusion gene activity is inducible and observed only following tamoxifen administration. When $sup_allele mice are bred with mice containing loxP-flanked sequence, tamoxifen-inducible, Cre-mediated recombination will result in deletion of the floxed sequences in the cells of the offspring";

	  } elsif ($allele =~/tm(\d{1})\(/) {
	    print "Found Allele!!!! $allele\n";
	    # Example targeted mutation 1, Wellcome Trust Sanger Institute
	    my $tm_att = $1;
	    my $serial_number = $tm_attempts{$tm_att};
	    $mutation_type = 'TN/TC';
	  } elsif ($allele =~/tm(1\.\d{1})\(/) {

	    my $tm_att = $1;
	    $tm_att =~ s/\./_/g;
	    #my $serial_number = $tm_attempts{$tm_att};
	    $mutation_type = 'TM/$tm_att';

	  } else {
	    print "ERROR: Allele type is unknown " . $allele . "\n";
	    print $_ . "\n";
	    $unknown_alleles{$colony_prefix} = $allele;
	    next;
	  }

	  #next;

	  # In past lines went GLT and were flagged as to_emma=1 and hence were imported into the EMMA database. Subsequently if these were shown to have failed QC the status was left the same but the to_emma flag set to 0 and these line records should then be set to Retracted (R) in the EMMA database
	  #
	  if ($consortium == "/IMPC") {
	    print "\tAlllele: $urlallele\n";
	    my $NEWgenetic_description  = $genetic_description;
	    $NEWgenetic_description =~ s/martsearch_ikmc_project\/martsearch\/ikmc_project\/$project_id/data\/alleles\/$mgi_accession_id/;
	    print "OLD GENETIC DESC:         $genetic_description\n\n";
	    print "NEW GENETIC DESC:         $NEWgenetic_description\n\n";
	  } 
	  if ($emma == 1) {
	    $str_access = 'P';
	  } else {
	    $str_access = 'R';
	  }

	  print "STR ACCESS:         $str_access\n";
	  print "\tPROJECT ID: $project_id\n";
	  # calculate funding source based on date 17 after 1st July 2010 to 31st December 2011, 18 in 2012
	  # 45 | I3-p1    | 1st reporting period of Infrafrontier-I3 WP4
	  # 46 | I3-p2    | 2nd reporting period of Infrafrontier-I3 WP4    
	  # 47 | I3-p3    | 3rd reporting period of Infrafrontier-I3 WP4  

	  # http://perldoc.perl.org/functions/localtime.html
	  # ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)
	  # year = 1900 + timeData[5]

	  my @timeData = localtime(time);

	  #
	  # check the distribution center exists!
	  # added on the 26 March 2013 and discussed with Sabine
	  #
	  if ($distribution_center eq 'SANG' || (length($distribution_center) == 0 && $repository eq 'SANG')) {
	    # Default funding source for SANG: 'MGP'
	    $source = 34;
	  } elsif ($timeData[5] == 111 || ($timeData[5] == 110 && $timeData[4] >= 6)) {
	    $source = 17;
	  } elsif ($timeData[5] == 112) {
	    $source = 18;
	  } elsif ($timeData[5] == 113 || ($timeData[5] == 114 && $timeData[4] < 6)) {
	    # 01 Jan 2013 until 30 June 2014: 'I3-p1'
	    $source = 45;
	  } elsif ($timeData[5] == 115  || ($timeData[5] == 114 && $timeData[4] >= 6)) {
	    # 01 July 2014 until 31 Dec 2015: 'I3-p2'
	    $source = 46;
	  } elsif ($timeData[5] == 116) {
	    # 01 Jan 2016 until 31 Dec 2016: 'I3-p3'
	    $source = 47;
	  } else {
	    $source = 5;
	  }

	  if ($distribution_center eq 'SANG' || ($repository eq 'SANG' && $distribution_center eq 'IMG' && ($pipeline eq 'Sanger MGP' || $pipeline eq 'MGP'))) {
	    # MGP
	    $project = 5;
	  } else {
	    # used to be EUCOM (3), become LS from now on
	    $project = 7;
	  }
		    
	  # definine the current time as this is used to record when a record is first inserted from iMits into the EMMA database
	  my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	  my $year = 1900 + $yearOffset;
	  $month = $month + 1;
	  $theTime = "$year-$month-$dayOfMonth";
		    
	  # CHECK REPOSITORY
	  # For lines reassigned from the original repository to a new distribution centre 
	  # need to create a new record in the EMMA database representing this line at the 
	  #new centre and set the status of the original one to N(Not for distribution

	  # flag to indicate whether distribution_center different from repository
	  my $dist_rep_diff = ($distribution_center ne $repository);

	  if (defined($distribution_center) && length($distribution_center) > 0 && $distribution_center ne $repository) {

	    # we will rely on the MI DATE and the COLONY NAME once all the colony names will be loaded.
	    my $sql = undef;

	    if (defined($colony_prefix) && length($colony_prefix) > 0) {
	      $sql = "select freezing_started, archived, notes, breeding, archiving_method_id, males, females, male_bg_id, female_bg_id, embryo_state, maintenance, mutant_viable, mutant_fertile, require_homozygous, immunocompromised, current_sanitary_status, animal_husbandry, genotype_file from laboratories, archive, strains left join residues on res_id=residues.id where archive_id=archive.id and code_internal='$epd_id' and id_labo=lab_id_labo and submitted='$mi_date' and code = '$repository' and colony_prefix = '$colony_prefix'";
	    } else {
	      $sql = "select freezing_started, archived, notes, breeding, archiving_method_id, males, females, male_bg_id, female_bg_id, embryo_state, maintenance, mutant_viable, mutant_fertile, require_homozygous, immunocompromised, current_sanitary_status, animal_husbandry, genotype_file from laboratories, archive, strains left join residues on res_id=residues.id where archive_id=archive.id and code_internal='$epd_id' and id_labo=lab_id_labo and submitted='$mi_date' and code = '$repository'";
	    }
			
	    my ($freezing_started,$archived,$notes,$breeding,$archiving_method_id,$males,$females,$male_bg_id,$female_bg_id,$embryo_state,$maintenance,$mutant_viable,$mutant_fertile,$require_homozygous,$immunocompromised,$current_sanitary_status,$animal_husbandry, $genotype_file) = &execute_query($sql, $emma_dbh);

	    $sql = "select id_str from laboratories, archive, strains where archive_id=archive.id and code_internal='$epd_id' and id_labo=lab_id_labo and submitted='$mi_date' and code = '$distribution_center' and str_access = 'P' and genotype_file is null"; # and colony_prefix = '$colony_prefix'";
	    my ($new_id_str) = &execute_query($sql, $emma_dbh);
	    #my $new_geno_file = "EM0".$new_id_str."_geno.pdf";
	    my $genoFileEmmaIdPrefix;
	if ($new_id_str < 10000) {
        $genoFileEmmaIdPrefix = 'EM:0';
      } else {
              $genoFileEmmaIdPrefix = 'EM:';
              }
		my $new_geno_file = $genoFileEmmaIdPrefix.$new_id_str."_geno.pdf";
              if ($new_id_str && $genotype_file) {
	      system("cp /nfs/panda/emma/genotype_protocols/$genotype_file /nfs/panda/emma/genotype_protocols/$new_geno_file");
	      system("svn add /nfs/panda/emma/genotype_protocols/$new_geno_file --username philw --password 9IDP26TI");
	      system("svn commit -m \"\" /nfs/panda/emma/genotype_protocols/$new_geno_file --username philw --password 9IDP26TI");
	    }

	    if ($repository eq 'SANG') {
	      if ($maintenance) {
				#next if ($maintenance =~ /^Current.+/);
		my @maintenance = split(/Current/,$maintenance);
		$maintenance = $maintenance[0];
		$maintenance = $maintenance." Current information may be viewed at the <a href=\"http://www.sanger.ac.uk/mouseportal/search?query=$marker_symbol\">Sanger mouse portal</a>: Viability at weaning, Fertility, General Observations";
	      } else {
		$maintenance = "Current information may be viewed at the <a href=\"http://www.sanger.ac.uk/mouseportal/search?query=$marker_symbol\">Sanger mouse portal</a>: Viability at weaning, Fertility, General Observations";
	      }
	    }
	    if ($material_deposited ne 'Live mice') { # frozen material sent
	      $source = 35;
	      $str_status = 'ARRD';
	      #$str_status = 'ARCHD';# change decided at IT meeting Sept 2011 but realised it causes problems in Dec 2011 with the automatic emails
	      #print STDERR "Execute $material_deposited emma update for Marker:$marker_symbol\n";
	      #print STDERR "Embryo state: $embryo_state\n";

	      &update_emma_database("$strain_name$dist_ilar_code",
				    $mutation_type,
				    $marker_symbol,
				    $allele,
				    $cell_line_bg,
				    $distribution_center,
				    $genetic_description,
				    $tm_es_line,
				    $cell_line_bgsymbol,
				    $mta_file,
				    $str_status,
				    $epd_id,
				    $mi_date,
				    $genotype_confirmed,
				    $maintained_bg,
				    $maintained_bgsymbol,
				    $str_access,
				    $phenotype_description,
				    $full_allele_name,
				    $repository,
				    $source,
				    $project,
				    $consortium,
				    undef, # theTime
				    $material_deposited,
				    $freezing_started,
				    $archived,
				    $notes,
				    $breeding,
				    $archiving_method_id,
				    $males,
				    $females,
				    $male_bg_id,
				    $female_bg_id,
				    $embryo_state,
				    $maintenance,
				    $mutant_viable,
				    $mutant_fertile,
				    $require_homozygous,
				    $immunocompromised,
				    $current_sanitary_status,
				    $animal_husbandry,
				    $colony_prefix,
				    $dist_rep_diff,
				    $impc_data_available,
				    $phenotype_attempts
				   );
	    } else {
	      $str_status = 'ARRD'; # changed back from ARCHD to ARRD in Feb 2012 on Sabine's advice
	      #print STDERR "Execute $material_deposited emma update for Marker:$marker_symbol\n";
	      &update_emma_database("$strain_name$dist_ilar_code",
				    $mutation_type,
				    $marker_symbol,
				    $allele,
				    $cell_line_bg,
				    $distribution_center,
				    $genetic_description,
				    $tm_es_line,
				    $cell_line_bgsymbol,
				    $mta_file,
				    $str_status,
				    $epd_id,
				    $mi_date,
				    $genotype_confirmed,
				    $maintained_bg,
				    $maintained_bgsymbol,
				    $str_access,
				    $phenotype_description,
				    $full_allele_name,
				    $repository,
				    $source,
				    $project,
				    $consortium,
				    undef,
				    $material_deposited,
				    undef, # freezing started
				    undef, # archived
				    undef, # notes
				    undef, # breeding
				    undef, # archiving method id
				    undef, # males
				    undef, # females
				    undef, # male_bg_id
				    undef, # female_bg_id
				    undef, # embryo_state
				    $maintenance,
				    $mutant_viable,
				    $mutant_fertile,
				    $require_homozygous,
				    $immunocompromised,
				    $current_sanitary_status,
				    $animal_husbandry,
				    $colony_prefix,
				    $dist_rep_diff,
				    $impc_data_available,
				    $phenotype_attempts
				   );
	    }
	    $str_access = 'N';
	  }
		    
	  if ($str_access eq 'N') {
	    $str_status = 'ARCHD';
	  } else {
	    $str_status = 'ARRD';
	  }

	  print "Execute global emma update for Marker:$marker_symbol with access $str_access and status $str_status\n";

	  &update_emma_database($strain_name,
				$mutation_type,
				$marker_symbol,
				$allele,
				$cell_line_bg,
				$repository,
				$genetic_description,
				$tm_es_line,
				$cell_line_bgsymbol,
				$mta_file,
				$str_status,
				$epd_id,
				$mi_date,
				$genotype_confirmed,
				$maintained_bg,
				$maintained_bgsymbol,
				$str_access,
				$phenotype_description,
				$full_allele_name,
				$repository,
				$source,
				$project,
				$consortium,
				$theTime, 
				undef ,		 # material deposited
				undef,		 # freezing started
				undef,		 # archived
				undef,		 # notes
				undef,		 # breeding
				undef,		 # archiving method id
				undef,		 # males
				undef,		 # females
				undef,		 # male_bg_id
				undef,		 # female_bg_id
				undef,		 # embryo_state
				$maintenance,	 # maintenance
				undef ,		 # mutant_viable
				undef ,		 # mutant_fertile
				undef ,		 # require_homozigous
				undef ,		 # immunocompromised
				undef , # current_sanitary_status
				undef , # animal husbandry 
				$colony_prefix,
				$dist_rep_diff,
				$impc_data_available,
				$phenotype_attempts
			       ); #$maintenance);

		    
	  print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

	}			# genotype confirmed

      }				# foreach entry
	    
    }				# decoding worked
	
    $total_entries += $nb_entries;
	
    last if (defined($colony_arg) || $nb_entries < $per_page || $current_page == $last_page); #$current_page > 3);

    $current_page++;

  }				# while 1   

  print STDERR "\n======================\nNumber of CRISPR/CAS9 alleles: " . scalar(keys(%endonuclease_mediated_mutations)) . "\n";
  foreach my $emm (keys %endonuclease_mediated_mutations) {
    print STDERR $emm . " " . $endonuclease_mediated_mutations{$emm} . "\n";
  }
  print STDERR "\n======================\nNumber of unspecified alleles: " . scalar(keys(%unknown_alleles)) . "\n";
  foreach my $ua (keys %unknown_alleles) {
    print STDERR $ua . " " . $unknown_alleles{$ua} . "\n";
  }    
  print STDERR "\n======================\nNumber of STOCK strains: " . $countStock . "\n";
  foreach my $k (keys %stockStrains) {
    print STDERR $k . " " . $stockStrains{$k} ."\n";
  }
  foreach my $l (@stocks) {
    print STDERR $l . "\n";
  }
    
  print "Nb of page parsed: " . ($current_page-$start_page+1) . "\tNb of entries parsed $total_entries\n";
    
  exit 1;
    
}

sub update_emma_database{
  my ($strain_name,
      $mutation_type,
      $marker_symbol,
      $allele,
      $cell_line_bg,
      $current_repository,
      $genetic_description,
      $tm_es_line,
      $cell_line_bgsymbol,
      $mta_file,
      $str_status,
      $epd_id,
      $mi_date,
      $genotype_confirmed,
      $maintained_bg,
      $maintained_bgsymbol,
      $str_access,
      $phenotype_description,
      $full_allele_name,
      $original_repository,
      $source,
      $project,
      $consortium,
      $theTime,
      $material_deposited,
      $freezing_started,
      $archived,
      $notes,
      $breeding,
      $archiving_method_id,
      $males,
      $females,
      $male_bg_id,
      $female_bg_id,
      $embryo_state,
      $maintenance,
      $mutant_viable,
      $mutant_fertile,
      $require_homozygous,
      $immunocompromised,
      $current_sanitary_status,
      $animal_husbandry,
      $colony_name,
      $b_dist_center_differs,
      $impc_data_available,
      $phenotype_attempts
     ) = @_;
  $genetic_description =~ s/\"//g;
  $phenotype_description =~ s/\"//g;
  my $sql;
  my $rtool_id;

  print "=== EMMA_UPDATE PROCEDURE START Original rep: $original_repository Current rep: $current_repository Genotype? $genotype_confirmed ===\n";
  # SELECT PERSON/LAB from people and laboratories tables
  my ($id_per_contact,$id_labo) = &execute_query("SELECT id_per, id_labo from people, laboratories where people.lab_id_labo = laboratories.id_labo and people.surname= laboratories.name and laboratories.code=" . $emma_dbh->quote($current_repository), $emma_dbh);
  my ($id_per) = &execute_query("SELECT id_per from people, laboratories where people.lab_id_labo = laboratories.id_labo and people.surname= laboratories.name and laboratories.code=" . $emma_dbh->quote($original_repository), $emma_dbh);

  # SELECT THE STRAIN WITH SAME EPD/COLONY PREFIX/MI DATE/LABO
  my $id_str = &check_existing_strain($strain_name, $epd_id, $colony_name, $mi_date, $id_labo, 1, 1);
  print "\n\n At line 1148 \n\n";
  my $iCountRtool = &execute_query("SELECT COUNT(str_id_str) FROM rtools_strains WHERE str_id_str='$id_str'",$emma_dbh);
  if ($phenotype_attempts > 0) {
    print "PHENOTYPING ATTEMPTS IS GREATER THAN 0 = $phenotype_attempts\n\n";
    #set id_str to null as do not want to update tma1 strain
    #$id_str = undef;
  }
    
  # STOP HERE if there no strain and no genotype confirmed mice 
  if (!$id_str && !$genotype_confirmed) {
    print "No existing strain identifier found. Exiting procedure.\n";
    return $id_str;
  }

  if ($debug) {
    print "EMMA_UPDATE PROCEDURE\n";
    if ($id_str) {
      print "emma_update:\tUpdate $colony_name $id_str\tDistribution access:$str_access\n";
    } else {
      print "emma_update:\tCreate $colony_name\tDistribution access:$str_access\n";
    }
  }

  # OTHERWISE CREATE NEW BACKGROUND RECORDS WHERE NEEDED
  my ($id_bg) = &execute_query("SELECT id_bg FROM backgrounds WHERE name = '$cell_line_bg'", $emma_dbh);
  if (!$id_bg) {
    my $sql = "INSERT INTO backgrounds (name,symbol) VALUES ('$cell_line_bg','$cell_line_bgsymbol')";
    &execute_query($sql, $emma_dbh);
  }
  ($id_bg) = &execute_query("SELECT id_bg FROM backgrounds WHERE name = '$maintained_bg'", $emma_dbh);
  if (!$id_bg) {
    $sql = "INSERT INTO backgrounds (name,symbol) VALUES ('$maintained_bg','$maintained_bgsymbol')";
    &execute_query($sql, $emma_dbh);
  }

  # CREATE/UPDATE NEW STRAINS AND LINK TABLE RECORDS IF APPROPRIATE FOR GENOTYPE CONFIRMED MICE
  # REVERT TO RETRACTED for not Genotype Confirmed one
    
  if ($id_per_contact && $id_per) { # excludes UCD and BCM lines that slip through
	
    # filtering empty strings / patch

    $males = ($males eq '' || !defined($males)) ? 'NULL' : $emma_dbh->quote($males);
    $females = ($females eq '' || !defined($females)) ? 'NULL' : $emma_dbh->quote($females);
    $embryo_state = ($embryo_state eq '' || !defined($embryo_state)) ? 'NULL' : $emma_dbh->quote($embryo_state);
    $male_bg_id = ($male_bg_id eq '') ? 'NULL' : $emma_dbh->quote($male_bg_id);
    $female_bg_id = ($female_bg_id eq '') ? 'NULL' : $emma_dbh->quote($female_bg_id);

    # NEW STRAIN
	
    if (!$id_str && $genotype_confirmed) {
      #print STDERR "NO strain ID!: $current_repository = $id_per_contact ;  $original_repository = $id_per\n";
      if ($theTime eq '') {
	#print STDERR "NO time!\n";
	$sql = "INSERT INTO archive (submitted,evaluated,glt,lab_id_labo,freezing_started, archived, notes, breeding, archiving_method_id, males, females, male_bg_id, female_bg_id, embryo_state) VALUES ('$mi_date','$mi_date',NULL,$id_labo,'$freezing_started','$archived','$notes','$breeding','$archiving_method_id',$males,$females,$male_bg_id,$female_bg_id,$embryo_state)";
      } else {
	#print STDERR "time = $theTime\n";
	$sql = "INSERT INTO archive (submitted,evaluated,glt,lab_id_labo,freezing_started, archived, notes, breeding, archiving_method_id, males, females, male_bg_id, female_bg_id, embryo_state) VALUES ('$mi_date','$mi_date','$theTime',$id_labo,'$freezing_started','$archived','$notes','$breeding','$archiving_method_id',$males,$females,$male_bg_id,$female_bg_id,$embryo_state)";

      }
      &execute_query($sql, $emma_dbh);
      $sql = "INSERT INTO residues (current_sanitary_status,animal_husbandry) VALUES ('$current_sanitary_status','$animal_husbandry')";
      &execute_query($sql, $emma_dbh);

      my ($archive_id) = &execute_query("select max(id) from archive", $emma_dbh);
      my ($res_id) = &execute_query("SELECT max(id) FROM residues", $emma_dbh);

      my $available_to_order = 'no';
      if ($current_repository eq 'MRC') {
	$available_to_order = 'yes';
      }

      $sql = "INSERT INTO strains (available_to_order, res_id, code_internal, name, bg_id_bg, per_id_per, per_id_per_contact, archive_id, charact_gen, pheno_text, str_status,str_access,str_type,mta_file,maintenance, mutant_viable, mutant_fertile, require_homozygous, immunocompromised, colony_prefix, ls_consortium) VALUES ('$available_to_order','$res_id','$epd_id','$strain_name','3270','$id_per','$id_per_contact','$archive_id',\"$genetic_description\",\"$phenotype_description\",'$str_status','$str_access','MSR','$mta_file','$maintenance','$mutant_viable','$mutant_fertile','$require_homozygous','$immunocompromised', '$colony_name', '$consortium')";
      &execute_query($sql, $emma_dbh);
      my ($idStr) = &execute_query("SELECT id_str FROM strains WHERE name ='$strain_name'",$emma_dbh);
my $emmaIdPrefix = '';
if ($idStr < 10000) {
	$emmaIdPrefix = 'EM:0';
#$sql = "UPDATE strains set emma_id = concat('EM:0',cast(id_str as char)) WHERE name ='$strain_name'";
      } else {
	$emmaIdPrefix = 'EM:';
}
$sql = "UPDATE strains set emma_id = concat('$emmaIdPrefix','$idStr') WHERE name ='$strain_name'";

&execute_query($sql, $emma_dbh);
	    
      my ($new_id_str) = &execute_query("SELECT max(id_str) FROM strains", $emma_dbh);
      $rtool_id = '9';
      if ( index($consortium, 'EUCOMMToolsCre' ) != -1) {
	$rtool_id = '10';
      } 
      $sql = "INSERT INTO rtools_strains (str_id_str, rtls_id) VALUES ($new_id_str,$rtool_id)";
      &execute_query($sql, $emma_dbh);
	    
      $sql = "INSERT INTO availabilities_strains (str_id_str, avail_id) VALUES ($new_id_str,'3')";
      &execute_query($sql, $emma_dbh);
	    
      $sql = "INSERT INTO projects_strains (str_id_str, project_id) VALUES ($new_id_str,'$project')";
      &execute_query($sql, $emma_dbh);

      print STDERR "$strain_name reporting period: $source\n";
      $sql = "INSERT INTO sources_strains (str_id_str, sour_id, lab_id_labo) VALUES ($new_id_str,'$source',$id_labo)";
      &execute_query($sql, $emma_dbh);

      # SELECT id_str (in this case, we retrieve the strain with a different lab id to reassign )
      # &execute_query("SELECT id_str FROM strains, archive WHERE archive_id = id and code_internal = '$epd_id' AND submitted = '$mi_date' and lab_id_labo != $id_labo");
      my ($original_id_str) = &check_existing_strain($strain_name, $epd_id, $colony_name, $mi_date, $id_labo, 0, 1);

      if ($original_id_str) {
if ($original_id_str < 10000) { 
	$emmaIdPrefix = 'EM:0';
} else {
	$emmaIdPrefix = 'EM:';
}
	$sql = "UPDATE web_requests SET str_id_str = '$new_id_str', strain_id = CONCAT('$emmaIdPrefix','$new_id_str'), lab_id_labo = $id_labo WHERE str_id_str='$original_id_str' AND req_status ='TO_PR'";
	&execute_query($sql, $emma_dbh);
      }

      # this will be used later for the mutations.
      $id_str = $new_id_str;

    }				# NEW STRAIN 
	
    # EXISTING STRAIN
    elsif ($id_str) {

      # Update strain records
      # There are 2 cases: 1. genotype confirmed: update the record normally. 2. Micro-injection aborted: the strain access should be R (for RETRACTED)

      # since the micro injection date is the same we also update the colony_prefix otherwise we have to fix things manually!                                                
      # We should update the database only if something has changed! 
      # below line added by phiw to resolve reccurring issue with 2 lines being available for MDCK/Cyflp1                          
      $str_access = 'R' if ($colony_name eq 'MDCK' && $id_str eq '6774');
      print STDERR "Update strain record of $id_str\n";
      my ($current_strain_name, $current_strain_access, $current_strain_status, $archive_id, $current_mta_file) = &execute_query("SELECT name, str_access, str_status, archive_id, mta_file FROM strains WHERE id_str=$id_str", $emma_dbh);

      # first time we run the script we update the colony_name                                                                                                                
      # however this will be removed since we'll check the colony prefix later                                                                                                
      print STDERR "Update strain colony_prefix to '$colony_name'\n";
      $sql = "UPDATE strains SET colony_prefix='$colony_name' WHERE id_str='$id_str'";
      &execute_query($sql, $emma_dbh);

      # Update consortium information associated to the strains
      print STDERR "Update strain consortium information to '$consortium'\n";
      $sql = "UPDATE strains SET ls_consortium='$consortium', charact_gen=\"$genetic_description\" WHERE id_str='$id_str'";
      &execute_query($sql, $emma_dbh);
      #Update rtools
      if ( index($consortium, 'EUCOMMToolsCre' ) != -1 && $iCountRtool >= 1) {
	print STDERR "Inserting rtools value of 10 for EUCOMMTOOLSCre lines '$consortium'\n";
	$rtool_id = '10';

	$sql = "UPDATE rtools_strains SET  rtls_id=$rtool_id WHERE str_id_str=$id_str AND rtls_id=9";
      }
      &execute_query($sql, $emma_dbh);
      # Check MTA information
      if (defined($mta_file) && defined($current_mta_file) && $mta_file ne $current_mta_file) {
	$sql = "UPDATE strains SET mta_file='$mta_file' WHERE id_str='$id_str'";
	&execute_query($sql, $emma_dbh);
      }

      # Update IMPC information
      &set_impc_phenotype_data_exists($colony_name, $impc_data_available,$phenotype_attempts);
      print "IMPC DATA HAS BEEN UPDATED AT LINE 1268 colony prefix is $colony_name\n\nIMPC data availability is $impc_data_available\n\npheno attempts is $phenotype_attempts\n\n";
      if ($strain_name ne $current_strain_name) {
	print STDERR "Set new strain name for $id_str\t'$current_strain_name' to '$strain_name' (colony_prefix = '$colony_name')\n";
	$sql = "UPDATE strains SET name='$strain_name' WHERE id_str='$id_str'";	# need to update name even if str_access = N or C to keep below code working                                                                                                     

	&execute_query($sql, $emma_dbh);
      }

      ## Check whether the distribution center was changed
      ## This happens sometimes for ICS
      ## In this case we should hide the 
      my $other_id_str = &check_existing_strain($strain_name, $epd_id, $colony_name, $mi_date, $id_labo, 0, 1);
      print STDERR "Found a different id str: $other_id_str\n" if ($other_id_str && $id_str != $other_id_str && $debug);

      ## rule for Sanger/patch when genotyped confirmed
      if ($genotype_confirmed) {
	if ($b_dist_center_differs == 0 && $current_repository eq 'SANG' && $str_access eq 'P' && $current_strain_status eq 'ARCHD') {
	  $sql = "UPDATE strains SET str_status = 'ARRD' WHERE id_str='$id_str'";
	  &execute_query($sql, $emma_dbh);
	}
		
	if ($current_strain_access ne 'N' && $current_strain_access ne 'C') {
	  $sql = "UPDATE strains SET charact_gen=\"$genetic_description\", pheno_text=\"$phenotype_description\", bg_id_bg=3270, per_id_per='$id_per',str_access='$str_access',maintenance='$maintenance', mutant_viable='$mutant_viable', mutant_fertile='$mutant_fertile', require_homozygous='$require_homozygous', immunocompromised='$immunocompromised' WHERE id_str='$id_str' AND str_access != 'N' AND str_access != 'C'";
	  &execute_query($sql, $emma_dbh);
	} elsif ($current_strain_access eq 'N' && $str_access eq 'P') {
	  print STDERR "Fixing strain access to P for $id_str\t$current_strain_name\n";
	  $sql = "UPDATE strains SET charact_gen=\"$genetic_description\", pheno_text=\"$phenotype_description\", bg_id_bg=3270, per_id_per='$id_per',str_access='$str_access',maintenance='$maintenance', mutant_viable='$mutant_viable', mutant_fertile='$mutant_fertile', require_homozygous='$require_homozygous', immunocompromised='$immunocompromised' WHERE id_str='$id_str'";
	  &execute_query($sql, $emma_dbh);
	}
      } else {
	## when not genotyped confirmed (and was genotype confirmed before)
	if ($current_strain_access eq 'P') {
	  print STDERR "ACCESS $current_strain_access REVERTED TO R for $id_str\t$marker_symbol\t$strain_name\t$colony_name\n";
	  $sql = "UPDATE strains SET str_access = 'R' WHERE id_str='$id_str'";
	  &execute_query($sql, $emma_dbh);
	}
      }
    }
  }
    
  if ($genotype_confirmed) {

    # CREATE/UPDATE GENE, ALLELE AND MUTATION RECORDS
    # All gene symbol updates should now come from the MGI BioMart so don't insert if a record already exists
    $log->debug("CREATE/UPDATE GENE, ALLELE AND MUTATION RECORDS\n");

    # We don't use the strain name but the strain ID
    my ($strain_id_gene) = &execute_query("SELECT id_gene FROM genes g, alleles a, mutations m, mutations_strains ms WHERE id_gene = gen_id_gene and alls_id_allel=id_allel and id=mut_id and ms.str_id_str='$id_str'", $emma_dbh);
    my ($id_gene) = &execute_query("SELECT id_gene FROM genes WHERE symbol ='$marker_symbol'", $emma_dbh);

    # if the gene does not exist
    # could be because it's an obsolete symbol?

    if (!$id_gene) {
      $sql = "INSERT INTO genes (symbol, name) VALUES ('$marker_symbol','$marker_symbol')";
      &execute_query($sql, $emma_dbh);
      $log->info("Gene entry didn't exists for strain $id_str: new marker symbol inserted $marker_symbol\n");
    }
    ($id_gene) = &execute_query("SELECT id_gene FROM genes WHERE symbol ='$marker_symbol'", $emma_dbh);
    $log->debug("Gene $marker_symbol id is $id_gene\n");

    # Check allele - first by looking at the association to the strain and then in the allele table 
    my $id_allele;
    my ($strain_id_allele) = &execute_query("SELECT id_allel FROM alleles a, mutations m, mutations_strains ms WHERE alls_id_allel=id_allel and id=mut_id and ms.str_id_str='$id_str'", $emma_dbh);
    my ($id_allele) = &execute_query("SELECT id_allel FROM alleles WHERE alls_form = '$allele'", $emma_dbh);
    $log->debug("SELECT id_allel FROM alleles WHERE alls_form = '$allele'\n");

    if (!$id_allele && $id_gene) {
      $sql = "INSERT INTO alleles (alls_form, name, gen_id_gene) VALUES ('$allele','$full_allele_name','$id_gene')";
      &execute_query($sql, $emma_dbh);
      ($id_allele) = &execute_query("SELECT id_allel FROM alleles WHERE alls_form = '$allele'", $emma_dbh);
      $log->debug("Allele entry didn't exists for strain $id_str: new allele symbol inserted $allele\n");
    }
    $log->debug("Allele $allele id is $id_allele\n");

    # is this step necessary?
    my ($id_labo) = &execute_query("SELECT id_labo FROM laboratories WHERE code = '$current_repository'", $emma_dbh);
    #my $id_str = &check_existing_strain($strain_name, $epd_id, $colony_name, $mi_date, $id_labo, 1, 1);
	
    if ($id_str) {		# to exclude UCD lines and BCM lines
      my ($mutation_id_allele,$mut_id) = &execute_query("SELECT mutations.alls_id_allel, mutations.id FROM mutations, mutations_strains WHERE mutations.id=mutations_strains.mut_id and mutations_strains.str_id_str ='$id_str' AND alls_id_allel = '$id_allele'", $emma_dbh);
      #&execute_query("SELECT alls_id_allel FROM mutations, mutations_strains WHERE mutations.id=mutations_strains.mut_id and mutations_strains.str_id_str ='$id_str' AND alls_id_allel = '$id_allele'");
      if (!$mutation_id_allele || $mutation_id_allele != $id_allele) {
	($id_bg) = &execute_query("SELECT id_bg FROM backgrounds WHERE name = '$cell_line_bg'", $emma_dbh);
	my @mutation_types = split(/\//,$mutation_type);
	$log->debug($mutation_types[1] . " " . $strain_name . "\n") if ($mutation_types[1] eq 'KOF');

	# check whether this mutation already existed
	my ($mut_id) = &execute_query("SELECT id from mutations WHERE alls_id_allel='$id_allele' AND main_type='$mutation_types[0]' AND sub_type='$mutation_types[1]' AND bg_id_bg='$id_bg' AND tm_esline='$tm_es_line'", $emma_dbh);
	if (!$mut_id) {
	  $log->debug("Mutation does not exist for this gene will create one\n");
	  $sql = "INSERT INTO mutations (alls_id_allel, main_type, sub_type, bg_id_bg, tm_esline) VALUES ('$id_allele','$mutation_types[0]','$mutation_types[1]','$id_bg','$tm_es_line')";
	  &execute_query($sql, $emma_dbh);
	  ($mut_id) = &execute_query("SELECT max(id) FROM mutations", $emma_dbh);
	}
	my ($mut_id_str) = &execute_query("SELECT str_id_str FROM mutations_strains WHERE str_id_str = '$id_str' AND mut_id = '$mut_id'", $emma_dbh);
		
	if (!$mut_id_str) {
	  $sql = "DELETE FROM mutations_strains WHERE str_id_str = '$id_str'";
	  &execute_query($sql, $emma_dbh);
	  $log->debug("$sql\n");
	  $sql = "INSERT INTO mutations_strains (str_id_str,mut_id) VALUES ('$id_str','$mut_id')";
	  &execute_query($sql, $emma_dbh);
	  $log->debug("$sql\n");
	}

      } else {
	# Check the correct mutation is associated

	my @mutation_types = split(/\//,$mutation_type);
	if ($mutation_types[1] eq 'KOF') {
	  $log->debug($mutation_types[1]. " update mutation " .$strain_name . "\n");
	  $sql = "UPDATE mutations set sub_type = '$mutation_types[1]' WHERE alls_id_allel = '$mutation_id_allele' AND id = '$mut_id'";
	  &execute_query($sql, $emma_dbh);
	}
      }
    }
  }
  $log->debug("=== EMMA_UPDATE PROCEDURE END ===\n");
}

# Update the strain impc_phenotype_data_exists field given a specific colony
# Double check the colony id is not empty and contains characters                                                                                                                                                                                               
# to avoid side effects with empty string

sub set_impc_phenotype_data_exists { 
  my ($colony_prefix, $impc_data_available,$phenotype_attempts) = @_;
  $log->debug("Setting IMPC data - Colony prefix arg= $colony_prefix and IMPC data availability is set to $impc_data_available");
  if (defined($colony_prefix) && length($colony_prefix) > 0 && $impc_data_available == 1) {
    $log->debug("Marking $colony_prefix as having IMPC data\n");
    my $sql = "UPDATE strains set impc_phenotype_data_exists = 'yes' WHERE colony_prefix='$colony_prefix'";
    &execute_query($sql, $emma_dbh);
  } elsif ($phenotype_attempts > 0) {
    #	my $sql = "UPDATE strains set impc_phenotype_data_exists = 'yes' WHERE colony_prefix='$colony_prefix'";
    #	&execute_query($sql, $emma_dbh);
    my $sql = "UPDATE strains set impc_phenotype_data_exists = NULL WHERE colony_prefix='$colony_prefix'";
    &execute_query($sql, $emma_dbh);
  } else {
    my $sql = "UPDATE strains set impc_phenotype_data_exists = NULL WHERE colony_prefix='$colony_prefix'";
    &execute_query($sql, $emma_dbh);
  }
}

# Update the strain ls_consortium field given a specific colony 
# Double check the colony id is not empty and contains characters
# to avoid side effects with empty string

sub update_strain_consortium_information {
  my ($colony_prefix, $consortium) = @_;

  if ($colony_prefix && length($colony_prefix) > 0 && $consortium && length($consortium) > 0) {

    $log->debug("Marking $colony_prefix as IMPC strain\n");
    my $sql = "UPDATE strains set ls_consortium = '$consortium' WHERE colony_prefix='$colony_prefix'";
    &execute_query($sql, $emma_dbh);

  }

}

# This method will check whether a strain exists already
# Parameters:
# b_same_lab: indicates whether we check if it's the same lab or not
#
sub check_existing_strain{
  my($strain_name, $epd_id, $colony_name, $mi_date, $id_labo, $b_same_lab, $b_strict_colony_name) = @_;
  my $id_str = undef;

  $log->debug("\tcheck_existing_strain:\n");
  $log->debug("\t\tStrain name: $strain_name\n");
  $log->debug("\t\tEPD ID: $epd_id\n");
  $log->debug("\t\tColony name: $colony_name\n");
  $log->debug("\t\tMI date: $mi_date\n");
  $log->debug("\t\tid_labo: $id_labo\n");

  my $operator = ($b_same_lab == 1) ? '=' : '!=';

  if (defined($id_labo) && $id_labo ne '') {
    my ($labo_code, $labo_name) = &execute_query("select code, name from laboratories where id_labo = $id_labo", $emma_dbh);
	
    # check first the colony name 
    if (defined($colony_name) && length($colony_name) > 0 ) {
      my ($id_q) = &execute_query("SELECT s.id_str FROM strains s, archive a WHERE s.archive_id=a.id and s.code_internal= '$epd_id' AND s.colony_prefix = '$colony_name' AND a.lab_id_labo $operator $id_labo AND a.submitted='$mi_date'", $emma_dbh);
      $id_str = $id_q;
      if ($id_str) {
	$log->debug("check_existing_strain ($b_same_lab): strain $id_str for colony prefix $colony_name - $strain_name - $epd_id - $mi_date - $id_labo ($labo_code) found.\n");
	return $id_str;
      }
    } 
	
    # otherwise select the strain with this mi date
    # &execute_query("SELECT id_str FROM strains, archive WHERE archive_id = id and code_internal = '$epd_id' AND submitted = '$mi_date' and lab_id_labo != $id_labo");
    if (!$b_strict_colony_name) {
      my ($id_q) = &execute_query("SELECT s.id_str FROM strains s, archive a WHERE s.archive_id=a.id and s.code_internal= '$epd_id' AND a.lab_id_labo $operator $id_labo AND a.submitted='$mi_date'", $emma_dbh);
      $id_str = $id_q;
      $log->debug("check_existing_strain ($b_same_lab): strain $id_str (no colony prefix) - $strain_name - $epd_id - $mi_date - $id_labo ($labo_code) found.\n") if ($id_str);
    }
  }
    
  return $id_str;
}

sub reject_strain_access_in_emma {
  my ($colony_prefix) = @_;

  my $rejected = undef;

  if (defined($colony_prefix) && length($colony_prefix) > 0) {
    my $sql = "select id_str, name, str_access from strains where colony_prefix = '$colony_prefix'";
    my ($id_str, $name, $str_access) = &execute_query($sql, $emma_dbh);
    if ($id_str && $str_access ne 'R') {
      $log->debug("Reject_strain_access:\t$colony_prefix\t$id_str\t$name\t$str_access\n");
      $sql = "UPDATE strains set str_access = 'R' WHERE id_str ='$id_str'";
      &execute_query($sql, $emma_dbh);
      $rejected = 1;
    }
  }
  return $rejected;
}

sub convert_center {
  my ($imits_center) = @_;
    
    $imits_center =~ s/WTSI/SANG/;
    $imits_center =~ s/MRC - Harwell/MRC/;
    $imits_center =~ s/Harwell/MRC/;
    $imits_center =~ s/GSF/HMGU/;
    $imits_center =~ s/Monterotondo/CNR/;
    $imits_center =~ s/Oulu/UNIOULU/;
    $imits_center =~ s/INFRAFRONTIER-Oulu/UNIOULU/;
    $imits_center =~ s/INFRAFRONTIER-VETMEDUNI/VETMEDUNI/;
  return $imits_center;
}


sub compute_strain_nomenclature {
  my ($escell_strain, $test_cross, $back_cross, $marker_symbol, $allele, $repository, $ilar_code) = @_;

  my $rule_number = 0;
  my $strain_name = undef;

  if ($escell_strain eq '129S7' && ($test_cross eq '129S5' || $test_cross eq '129S5/SvEvBrd/Wtsi')) {
    $rule_number = 1;
    $strain_name = "129-$allele$ilar_code";
  } elsif ($escell_strain eq '129P2/OlaHsd' && $test_cross eq '129P2/OlaHsd' && $back_cross eq 'C57BL/6J') {
    $rule_number = 2;
    $strain_name = "B6J;129P2-$allele$ilar_code";
  } elsif ($escell_strain eq '129P2/OlaHsd' && $test_cross eq '129P2/OlaHsd' && $back_cross eq '129P2/OlaHsd') {
    $rule_number = 3;
    $strain_name = "129P2/OlaHsd-$allele$ilar_code";
  } elsif ($escell_strain eq '129P2/OlaHsd' && $test_cross eq 'C57BL/6Dnk' && $back_cross eq 'C57BL/6Dnk') {
    $rule_number = 4;
    $strain_name = "B6Dnk;129P2-$allele$ilar_code";
  } elsif ($escell_strain eq 'C57BL/6N' && ($test_cross eq '129S5' || $test_cross eq '129S5/SvEvBrd/Wtsi')) {
    $rule_number = 5;
    $strain_name = "129S5;B6N-$allele$ilar_code";
  } elsif ($escell_strain eq 'C57BL/6N' && ($test_cross eq 'B6JTyr<c-Brd>' || $test_cross eq 'C57BL/6Brd-Tyr<c-Brd>')) {
    $rule_number = 6;
    if ($back_cross eq 'C57BL/6Dnk') {
      $strain_name = "B6Dnk;B6Brd;B6N-Tyr<c-Brd> $allele$ilar_code";
    } else {
      $strain_name = "B6Brd;B6N-Tyr<c-Brd> $allele$ilar_code";
    }
  }
  # rule applied for 2 WTSI strains
  # [C57BL/6N-A<tm1Brd>/a][C57BL/6NTac][129S5/SvEvBrd/Wtsi
  elsif (($escell_strain eq 'C57BL/6N-A<tm1Brd>/a' && $test_cross eq 'C57BL/6NTac' && $back_cross eq '129S5/SvEvBrd/Wtsi' && $repository eq 'SANG') ||
	 ($escell_strain eq 'C57BL/6N' && $test_cross eq 'C57BL/6NTac' && $back_cross eq '129S5/SvEvBrd/Wtsi' && $marker_symbol eq 'Ywhae')) {
    $rule_number = 7;
    $strain_name = "129S5;B6N-$allele$ilar_code";
  } elsif ($escell_strain eq 'C57BL/6N-A<tm1Brd>/a' && ($test_cross eq 'B6JTyr<c-Brd>' || $test_cross eq 'C57BL/6Brd-Tyr<c-Brd>' || ($test_cross eq 'C57BL/6J-Tyr<c-Brd>;C57BL/6N' && $back_cross eq $test_cross))) {
    if ($repository eq 'SANG' || $repository eq 'CNR' || $repository eq 'ICS' || $repository eq 'CIPHE') {
      $rule_number = 8;
      if ($back_cross eq 'C57BL/6Dnk') {
	$strain_name = "B6Dnk;B6Brd;B6N-Tyr<c-Brd> A<tm1Brd> $allele$ilar_code";
      } else {
	$strain_name = "B6Brd;B6N-Tyr<c-Brd> A<tm1Brd> $allele$ilar_code";
      }
    } else {
      $rule_number = 9;
      $strain_name = "B6Brd;B6N-Tyr<c-Brd> $allele$ilar_code";
    }
  } elsif ($escell_strain eq 'C57BL/6N' && $test_cross eq 'C57BL/6Dnk' && $back_cross eq 'C57BL/6Dnk') {
    $rule_number = 10;
    $strain_name = "B6Dnk;B6N-$allele$ilar_code";
  } elsif ($escell_strain eq 'C57BL/6N' && $test_cross eq 'C57BL/6N' && ($back_cross eq 'C57BL/6N' || $back_cross eq '')) {
    $rule_number = 11;
    $strain_name = "C57BL/6N-$allele$ilar_code";
  } elsif ($escell_strain eq 'C57BL/6N' && $test_cross eq 'C57BL/6NTac') {
    $rule_number = 12;
    $strain_name = "C57BL/6NTac-$allele$ilar_code";
  } elsif (($escell_strain eq 'C57BL/6N' || $escell_strain eq 'C57BL/6N-A<tm1Brd>/a') &&
	   ($test_cross eq 'C57BL/6NTac' || $test_cross eq 'B6N' || $test_cross eq 'C57BL/6Dnk') &&
	   ($back_cross eq '' || $back_cross eq 'C57BL/6NTac')) {
    if ($escell_strain eq 'C57BL/6N-A<tm1Brd>/a' && ($repository eq 'SANG' || $repository eq 'ICS')) {
      $rule_number = 13;
      $strain_name = "C57BL/6N-A<tm1Brd> $allele$ilar_code";
    } else {
      print "we are at l1567 rule 14:" if $debug;
      $rule_number = 14;
      $strain_name = "C57BL/6N-$allele$ilar_code";
      my $find = "tm1a";
      my $replace = "tm1b";
      $find = quotemeta $find;
      #$strain_name =~ s/$find/$replace/g; 
    }
  } elsif (($escell_strain eq 'C57BL/6N-A<tm1Brd>/a') && ($test_cross eq 'C57BL/6NTac' || $test_cross eq '' || $test_cross eq 'C57BL/6N') && $back_cross eq 'C57BL/6NTac') {
    if ($repository eq 'SANG' || $repository eq 'CNR' || $repository eq 'Monterotondo' || $repository eq 'ICS' || $repository eq 'MRC') {
      $rule_number = 15;
      $strain_name = "B6NTac;B6N-A<tm1Brd> $allele$ilar_code";
      print "we are at l1576 rule 15" if $debug;
    } else {
      $rule_number = 16;
      $strain_name = "B6NTac;B6N-$allele$ilar_code";
    }
  } elsif (($escell_strain eq 'C57BL/6N-A<tm1Brd>/a') && $test_cross eq 'C57BL/6Dnk' && $back_cross eq 'C57BL/6Dnk') {
    if ($repository eq 'SANG' || $repository eq 'CNR' || $repository eq 'ICS') {
      $rule_number = 17;
      $strain_name = "B6Dnk;B6N-A<tm1Brd> $allele$ilar_code";
    } else {
      $rule_number = 18;
      $strain_name = "B6Dnk;B6N-$allele$ilar_code";
    }
  } elsif ($escell_strain eq 'C57BL/6N-A<tm1Brd>/a' && ($test_cross eq 'C57BL/6NCrl' || $test_cross eq '') && $back_cross eq 'C57BL/6NCrl') {
    $rule_number = 19;
    $strain_name = "B6NCrl;B6N-A<tm1Brd> $allele$ilar_code";
  } elsif ($escell_strain eq 'C57BL/6N-A<tm1Brd>/a' && $back_cross eq 'C57BL/6NJ') {
    #$rule_number = 20;
    $strain_name = "B6NJ;B6N-$allele$ilar_code";
  } elsif ($escell_strain eq 'C57BL/6N' && $test_cross eq 'C57BL/6NCrl' && $back_cross eq 'C57BL/6NCrl') {
    $rule_number = 20;
    $strain_name = "B6NCrl;B6N-$allele$ilar_code";
  } elsif (($escell_strain eq 'C57BL/6N-A<tm1Brd>/a') && ($test_cross eq 'C57BL/6N') && ($back_cross eq 'C57BL/6N')) {
    if ($repository eq 'SANG' || $repository eq 'CNR' || $repository eq 'ICS') {
      $rule_number = 21;
      $strain_name = "C57BL/6N-A<tm1Brd> $allele$ilar_code";
    } else {
      $rule_number = 22;
      $strain_name = "C57BL/6N-$allele$ilar_code";
    }
  } else {
    $rule_number = 23;
    $strain_name = "STOCK $allele$ilar_code";
    $log->debug("STOCK nomenclature for:\t");
    $log->debug("ES Cell:[$escell_strain]\tTest Cross:[$test_cross]\tBack Cross:[$back_cross]\n");
    $countStock++;
    $stockStrains{"ES Cell:[$escell_strain]\tTest Cross:[$test_cross]\tBack Cross:[$back_cross]"}++;
    push @stocks, "$marker_symbol\tES Cell:[$escell_strain]\tTest Cross:[$test_cross]\tBack Cross:[$back_cross]\t$repository";
  }

  $log->debug("$rule_number\tES CELL\tTEST CROSS\tBACK CROSS\n");
  $log->debug("\t[$escell_strain]\t[$test_cross]\t[$back_cross]\n");
  $log->debug($strain_name . "\n");

  return $strain_name;
}
         
sub compute_mta_file {
  my ($pipeline, $repository, $distribution_center, $consortium, $phenotyping_attempt, $allele, $mi_date) = @_;
    
  my ($year, $month, $day) = ($mi_date =~ m/^(\d+)\-(\d+)\-(\d+)/);
  print STDOUT "YEAR $year\n";
  print STDOUT "compute_mta_file: $pipeline, $repository, $distribution_center, $consortium, $phenotyping_attempt, $allele\n";
  my $mta_file = "MTA_EUCOMM_for_EMMA_$repository.pdf";

  # BUSINESS RULE 1: IKMC (iMits) strains produced by HMGU 
  # tm1a/tm1e: IF MI date before 25 Oct 2011 THEN set strains.mta_file = 'MTA_EUCOMM_for_EMMA_HMGU.pdf'
  # ELSE (i.e. if MI date is 25 Oct 2011 or later) THEN set strains.mta_file = 'HMGU-SMTA-EMMA_2014.pdf'
  # tm1b/tm1c: I understand there is no MI date for tm1b strains (but rather a Cre excision complete date). 
  # Anyway, ALL Cre or Frt excised strains produced by HMGU need to be linked to 'HMGU-SMTA-EMMA_2014.pdf'.      
  if ($distribution_center eq 'HMGU') {
  	$log->debug($allele);
    # Check allele nomenclature
    if ($allele =~ /tm\d+a/ || $allele =~ /tm\d+e/) {
      if ( $year < 2011 || ($year == 2011 && ($month < 10 || ($month == 10 && $day < 25)))) { # 25 Oct 2011
      	if ($pipeline eq 'EUCOMM') {
			$mta_file = 'MTA_EUCOMM_for_EMMA_HMGU.pdf'
		}
      } else {
		$mta_file = 'HMGU-SMTA-EMMA_2014.pdf'
      }
    } elsif ($allele =~/tm\d+b/ || $allele =~ /tm\d+c/ || $allele =~ /tm\d\.\d{1}/) {
      # tm1b: Reporter-tagged deletion allele (post-Cre)
      $mta_file = 'HMGU-SMTA-EMMA_2014.pdf';
    }
    $log->info($mta_file);
    return $mta_file;
  }
  
  # Business Rule 2:
  # Link WTSI lines to the following MTA files (strains.mta_file):
  # EUCOMM lines --> MTA_EUCOMM_for_EMMA_SANG.pdf
  # EUCOMMTools lines (KO and Cre) --> MTA-Sanger-Standard-form_EUCOMMTools.pdf
  # KOMP ESC/CRISPR/ncRNA lines --> MTA-SangerMP-Standard-Form.pdf
  
  if ($repository eq 'SANG' || $distribution_center eq 'SANG') {
  	if ($pipeline eq 'EUCOMM') {
  		$mta_file = "MTA_EUCOMM_for_EMMA_SANG.pdf"
  	} 
  	elsif (($pipeline eq 'EUCOMMTools' || $pipeline eq 'EUCOMMToolsCre')) { 
  		$mta_file = "MTA-Sanger-Standard-form_EUCOMMTools.pdf";
  	}
  	elsif ($pipeline eq 'KOMP-CSD') {
  		$mta_file = "MTA-Sanger-Standard-Form_KOMP.pdf";
  	}
  	elsif ($pipeline eq 'Sanger MGP') {
  		$mta_file = "MTA-SangerMP-Standard-Form.pdf";
  	} else {
  		die "WTSI lines: There are no business rule to associate an MTA file for pipeline $pipeline";
  	}
  	$log->info($mta_file);
  	return $mta_file;
  }

  # For all other centers
  # When we have read a phenotyping attempt
  if ($phenotyping_attempt) {

    if ($pipeline eq 'EUCOMM') {
      $mta_file = "MTA_EUCOMM_for_EMMA_$repository.pdf";

      if ($distribution_center eq 'MRC' && $repository ne 'SANG') {
	    $mta_file = "Master_Mouse_MTA_for_IMPC_Jan-2013.doc";
      } elsif ($distribution_center eq 'HMGU') {
	    $mta_file = "MTA_EUCOMM-tm1b_for_EMMA_HMGU.pdf";
      } elsif ($distribution_center eq 'ICS') {
	    $mta_file = "MTA_EUCOMM-tm1b_for_EMMA_ICS.pdf";
      }
    }

  } else {
    if (($pipeline eq 'EUCOMMTools' || $pipeline eq 'EUCOMMToolsCre')) { 
      if ( $repository eq 'CNR' ) {
	    $mta_file = "MTA-CNR-Standard-form_EUCOMMTools.pdf";
      }
    } elsif ($pipeline eq 'KOMP-Regeneron' && $repository eq 'IMG') {
      $mta_file = "MTA_IMG.pdf";
    } elsif ($consortium eq 'BaSH' && $repository eq 'MRC') {
      $mta_file = "Master_Mouse_MTA_for_IMPC_Jan-2013.doc";
    } elsif ($repository eq 'SEAT' && $distribution_center eq '') {
    } else {
      $mta_file = "MTA_EUCOMM_for_EMMA_$repository.pdf";
    }
  }
  print $mta_file;
  return $mta_file;
}

sub correct_background_information {
  my ($background) = @_;

  # typos in iMits fixing
  $background =~ s/C57BL\/6N-A\<sup\>tm1Br/C57BL\/6N-A\<tm1Brd\>\/a/;
  $background =~ s/C57BL6\/NTac/C57BL\/6NTac/;
  $background =~ s/C57B6NTac/C57BL\/6NTac/;
  $background =~ s/C57BL\/6NTAC/C57BL\/6NTac/;
  $background =~ s/C57bl\/6ntac/C57BL\/6NTac/;
  $background =~ s/129\/P2/129P2\/OlaHsd/;
  $background =~ s/C57BL\/6Ntac/C57BL\/6NTac/;
  $background =~ s/C57BL6\/NCrl/C57BL\/6NCrl/;
    
  $background = 'C57BL/6NTacDen' if ($background eq 'C57BL/6NTac/Den' || $background eq 'B6N/TacDen');
  $background = 'C57BL/6NTacUSA' if ($background eq 'C57BL/6NTac/USA' || $background eq 'B6N/TacUSA');
  $background = 'C57BL/6Dnk' if ($background eq 'C57BL/6NTac/Den' || $background eq 'C57BL/6NTacDen');
  $background = 'C57BL/6NTac' if ($background eq 'C57BL/6NTac/USA' || $background eq 'C57BL/6NTacUSA');

  $background = 'C57BL/6NCrl' if ($background eq 'C57BL6/NCrl');

  # This is the correct nomenclature
  $background = 'C57BL/6Brd-Tyr<c-Brd>' if ($background eq 'C57BL/6J-Tyr<c-Brd>');
  $background = 'B6BrdTyr<c-Brd>' if ($background eq 'B6JTyr<c-Brd>');
    
  return $background;
}

sub check_impc_data_exists {
  my ($colony_prefix) = @_;
  my $sql = "select count(distinct ls.colony_id) from live_sample ls join biological_sample bs on ls.id = bs.id join external_db edb on bs.db_id = edb.id join organisation o on o.id = bs.organisation_id where bs.sample_group = 'experimental' and edb.short_name = 'IMPC' and ls.colony_id = '$colony_prefix'";
  my ($c) = execute_query($sql, $komp2_dbh);
  return $c;
}

sub check_pipeline_project {
  my ($pipeline, $allele) = @_;
  my $allele_project = undef;

  if (defined($pipeline)) {
    if ($pipeline =~ /KOMP/ || $pipeline =~ /EUCOMM/) {
      $allele_project = ($pipeline eq 'KOMP-CSD' || $pipeline eq 'KOMP-Regeneron') ? "KOMP" : "EUCOMM";
    } elsif ($pipeline eq "Sanger MGP" && $allele =~ /KOMP/) {
      $allele_project = "KOMP";
    } elsif ($pipeline eq "Sanger MGP" && $allele =~ /MGP/) {
      $allele_project = "MGP";      
    } else {
      print "ERROR: unknown allele project for pipeline '$pipeline' and allele '$allele' - PLEASE FIX THIS - and restart the script\n";
      exit 1;
    }
  } else {
    if ($allele =~ /KOMP/) {
      $allele_project = 'KOMP';
    } elsif ($allele =~ /EUCOMM/) {
      $allele_project = 'EUCOMM';
    } else {
      print "ERROR: unknown allele project for allele $allele - PLEASE FIX THIS - and restart the script\n";
      exit 1;	    
    }
  }

  return $allele_project;

}

sub convert_center_information {
    my ($center) = @_;
    # use CIPHE for CIPHE, IMG for IMG, CNRS for SEAT
    $center =~ s/WTSI/SANG/;
    $center =~ s/MRC - Harwell/MRC/;
    $center =~ s/Harwell/MRC/;
    $center =~ s/GSF/HMGU/;
    $center =~ s/Monterotondo/CNR/; 
    $center =~ s/SEAT/CNRS/;
    $center =~ s/Oulu/UNIOULU/;
    $center =~ s/INFRAFRONTIER-Oulu/UNIOULU/;
    $center =~ s/Fleming/BSRC/;
    $center =~ s/INFRAFRONTIER-VETMEDUNI/VETMEDUNI/;

    return $center;
}

sub get_ilar_code {
  my ($center) = @_;

  return $ilar_codes{$center};

}

sub build_ilar_codes {

  $ilar_codes{'SANG'} = "Wtsi";
  $ilar_codes{'CNR'} = "Cnrm";
  $ilar_codes{'MRC'} = "H";
  $ilar_codes{'ICS'} = "Ics";
  $ilar_codes{'HMGU'} = "Ieg";
  $ilar_codes{'CNB'} = "Cnbc";
  $ilar_codes{'UNIOULU'} = "Oulu";
  $ilar_codes{'VETMEDUNI'} = "Biat";
  $ilar_codes{'CNRS'} = "Orl"; # the reason being SEAT and TAAM are the same place!  
  $ilar_codes{'CIPHE'} = "Ciphe"; # unlikely to happen but who knows?
  $ilar_codes{'IMG'} = "Ph";
  $ilar_codes{'BSRC'} = "Flmg";

}

#
# Check whether it's an EMMA repository
# Includes: "Harwell,HMGU,ICS,Monterotondo,WTSI,CIPHE,SEAT,IMG"
# Excludes: "APN, BCM, JAX, MARC, RIKEN BRC, TCP, UCD from production"
#

sub is_emma_repository {
    my ($repository) = @_;

    return ($repository eq 'GSF' || 
	    $repository eq 'HMGU' || 
	    $repository eq 'ICS' || 
	    $repository eq 'Monterotondo' || 
	    $repository eq 'MRC - Harwell' || 
	    $repository eq 'Harwell' || 
	    $repository eq 'WTSI' || 
	    $repository eq 'CIPHE' || 
	    $repository eq 'SEAT' || 
	    $repository eq 'IMG' ||
            $repository eq 'INFRAFRONTIER-Oulu' ||
            $repository eq 'INFRAFRONTIER-VETMEDUNI');
}

#
# This is critical code
# It checks whether the line should be actually deposited in EMMA
# For instance, the EUCOMMTOOLSCre lines were in embargo but it was decided
# later to release them
# For KOMP lines, it's important to check that the distribution center is a 
# valid distribution center
#
sub emma_compliant {
    my ($marker_symbol, $colony_prefix, $pipeline, $repository, $distribution_center, $phenotyping_attempt, $genotype_confirmed, $emma) = @_;

    my $marker_exclusion_test = !($marker_symbol eq 'Wls' && $colony_prefix eq 'gpr177 albino' && $emma == 0) ? 1 : 0;
    my $pipeline_compliance_test = ($pipeline eq 'Sanger MGP' || $pipeline eq 'MGP' || $pipeline eq 'CRISPR/CAS9') ? 1 : 0;

    my $eucommtools_cre_compliance_test = ($pipeline eq 'EUCOMMToolsCre' && ( 
					   $repository eq 'WTSI' || $repository eq 'Monterotondo' ) &&
					   ($distribution_center eq 'WTSI' ||
					    $distribution_center eq 'Harwell' ||
                                                $distribution_center eq 'Monterotondo' ||
						$distribution_center  eq 'HMGU')
					   ) ? 1 : 0;

    my $eucomm_compliance_test = (($pipeline eq 'EUCOMM' || $pipeline eq 'EUCOMMTools') && 
				  ($repository eq 'GSF' || 
				   $repository eq 'HMGU' || 
				   $repository eq 'ICS' || 
				   $repository eq 'Monterotondo' || 
				   $repository eq 'MRC - Harwell' || 
				   $repository eq 'Harwell' || 
				   $repository eq 'WTSI' || 
				   $repository eq 'CIPHE' || 
				   $repository eq 'SEAT' || 
				   $repository eq 'IMG' ||
                                   $repository eq 'INFRAFRONTIER-Oulu' ||
				   $repository eq 'INFRAFRONTIER-VETMEDUNI')) ? 1 : 0;

    my $komp_compliance_test = (($pipeline eq 'KOMP-CSD') && 
				($repository eq 'HMGU' || 
				 $repository eq 'ICS' || 
				 $repository eq 'Monterotondo' || 
				 $repository eq 'MRC - Harwell' || 
				 $repository eq 'Harwell' || 
				 $repository eq 'WTSI' || 
				 $repository eq 'CIPHE' || 
				 $repository eq 'SEAT' || 
				 $repository eq 'IMG' ||
				 $repository eq 'INFRAFRONTIER-Oulu')) ? 1 : 0;

  my $komp_compliance_test = (($pipeline eq 'KOMP-CSD') && 
			      ($repository eq 'HMGU' || 
			       $repository eq 'ICS' || 
			       $repository eq 'Monterotondo' || 
			       $repository eq 'MRC - Harwell' || 
			       $repository eq 'Harwell' || 
			       $repository eq 'WTSI' || 
			       $repository eq 'CIPHE' || 
			       $repository eq 'SEAT' || 
			       $repository eq 'IMG') && 
			      ($distribution_center eq 'HMGU' || 
			       $distribution_center eq 'ICS' || 
			       $distribution_center eq 'Monterotondo' ||
			       $distribution_center eq 'CNRS' || 
			       $distribution_center eq 'Oulu' || 
			       $distribution_center eq 'Harwell' ||
			       $distribution_center eq 'WTSI' || 
			       $distribution_center eq 'Fleming' || 
			       $distribution_center eq 'IMG' || 
			       $distribution_center eq 'VETMEDUNI'
			      )) ? 1 : 0; 

  my $komp_regeneron_compliance_test = ($pipeline eq 'KOMP-Regeneron' && $repository eq 'IMG') ? 1 : 0;
  my $cre_excised_compliance_test =  ($phenotyping_attempt && $emma && $genotype_confirmed) ? 1 : 0;
  if ($debug) {
	
    print "COMPLIANCE TEST:\n";
    print "\tCOLONY:$colony_prefix\n";
    print "\tPIPELINE:$pipeline (EUCOMM,EUCOMMTools,KOMP-Regeneron,KOMP-CSD,Sanger MGP)\n";
    print "\tREPOSITORY:$repository (GSF,HMGU,ICS,Monterotondo,MRC - Harwell, WTSI, CIPHE, SEAT, IMG)\n";
    print "\tDISTRIBUTION:$distribution_center (CNRS,Oulu,WTSI,Harwell, Fleming)\n";
    print "\tDISTRIBUTED BY EMMA: $emma\n";
    print "\tPIPELINE_COMPLIANCE: $pipeline_compliance_test\n";
    print "\tEUCOMM_COMPLIANCE: $eucomm_compliance_test\n";
    print "\tEUCOMTOOLSCRE_COMPLIANCE: $eucommtools_cre_compliance_test\n";
    print "\tKOMP_COMPLIANCE: $komp_compliance_test\n";
    print "\tKOMP_REGENERON_COMPLIANCE: $komp_regeneron_compliance_test\n";
    print "\tCRE_EXCISED_COMPLIANCE: $cre_excised_compliance_test\n";
  }
    
  my $compliance = 
    $cre_excised_compliance_test ||
      ($marker_exclusion_test && !$phenotyping_attempt &&
       ($eucomm_compliance_test || $komp_compliance_test || 
	$komp_regeneron_compliance_test ||
	$eucommtools_cre_compliance_test ||
	$pipeline_compliance_test));
  #	($pipeline eq 'Sanger MGP' && $distribution_center ne 'IMG') ||
    
  print "COMPLIANCE:$compliance\n";
  return $compliance;
}

sub execute_query{
  my ($sql, $dbh) = @_;
  $sql =~ s/=''/=null/g;
  $sql =~ s/''/null/g;
  if ($sql =~ /DELETE/ || $sql =~ /INSERT/ || $sql =~ /UPDATE/) {
    #print STDERR "$sql\n" ;

    my $sth = $dbh->prepare($sql) or die "Couldn't prepare statement $sql: " . $dbh->errstr;
    my $result = eval {
      $sth->{PrintError} = 0;
      $sth->execute(); #  or die "Couldn't execute statement $sql: " . $sth->errstr;
      1;
    };

    unless ($result) {
      print STDERR $dbh->errstr(),"\n";
      print STDERR "Couldn't execute statement $sql\n";
    }
	
  } else {
    $sql =~ s/= null/is null/g;
    print STDERR "$sql\n" ;
    my $sth = $dbh->prepare($sql);
    $sth->execute;
    my @result = $sth->fetchrow_array;
    return @result;
  }
}

