#    Copyright (C) 2015  Anthony Bowen, see LICENSE.txt for more details
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

# =============================================================================== Scientist Search Version 1.7 ============================================================================================= #
# Anthony Bowen
# anthony.bowen@med.einstein.yu.edu
# Laboratory of Arturo Casadevall
# Albert Einstein College of Medicine, Bronx NY
# Last Modified 1/26/15; Added a function to query PubMed for only publications with a given affiliation
#                        Added a function to record the number of unique first authors, rather than all authors
#
# This will find the number of unique Pubmed published scientists, retractions, and publications per year over a given time period.
    
#!/usr/bin/perl -w
use strict;
use warnings;

use Time::Piece;
use LWP::Simple;
use Parallel::ForkManager; #Pt1
use Data::Dumper;
use List::Util qw( min max );
use Statistics::Descriptive;

# ======================================================================= Variable Selections  ======================================================================================= #
my $version = "1.7";
my $file_directory = "/Volumes/Travelstar 1TB/Documents/Perl Temp/";       #"/Users/Tony/Documents/Casadevall Lab/PERL/Scientist Search/";
my $data_directory = "/Users/Tony/Documents/Casadevall Lab/PERL/Scientist Search/Data Reference Files/";
 
my $data_retrieval_mode = 'medline';            #Pt1; Select which formats to retrieve PubMed data (csv, medline, all)
my $MAX_PROCESSES = 0;                          #Pt1; Maximum number of processes to use for Parallel::ForkManager, change to 0 for NO fork to be done
my @years = (1809..2014);                       #Pt1; List of years to search (1809..2014)
my @journal_years = (1809..2014);               #Pt1; List of years to analyze journal publication numbers and impact factors
my $max_ID_search = 100000;                     #Pt1; Maximum number of pubmed IDs to retrieve per web query
my $max_CSV_search = 500;                       #Pt1; Maximum number of pubmed ID CSV data to retrieve per web query
my $max_medline_search = 500;                   #Pt1; Maximum number of pubmed ID medline data to retrieve per web query
my $publication_type = "journal article";       #Pt1; Publication type requirement for citations to search (all, journal article)
my $affiliation_restriction = "none";            #Pt1; Affiliation requirement for citations (none, USA, ...)
my $app_output = "summary";                     #Pt1; Print full arrays of authors per publication data per year or just the summary stats (full, summary)
my $keyword_mode = "none";                       #Pt1; Analyze 'all' keywords or only those matching the 'query' list or 'none' for medline publications
my @keyword_query =                             #Pt1; Keywords to search medline citations for listed keywords
    ("human","etiology","epidemiology","shigella");

# for each section, change to 1 to run or 0 to skip:
my $part1 = 1; my $part2 = 0; 	#1: Get Yearly Pubmed Data				   #2: 

# ======================================================================== Script Sections ================================================================================== #
my $date = Time::Piece->new->strftime('%m/%d/%Y');
#PART 1: ########### Get Yearly Pubmed Data #############
# This part of the script searches all @years for the number of pubmed publications, retractions, authors, and journals
# Input @years
# Output $year_pubmed_authors.tab, ss_data.tab
if ($part1 == 1) {
    my $pm = Parallel::ForkManager->new($MAX_PROCESSES);

YEAR: foreach my $y (0..$#years){
        $pm->start and next; # do the fork
        
        ### Skip years that have already been finished
        my $output_filename = "$file_directory"."$years[$y]"."_pubmed_authors.tab";
        my $csv_pmid_check = 0;
        my $medline_pmid_check = 0;
        if (!-e $output_filename){
            
            ### Get XML file for the first $max_ID_search pubmed IDs
            my $pubmed_search_url;
            if ($publication_type eq "all"){
                if ($affiliation_restriction eq "none"){
                    $pubmed_search_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=\"$years[$y]\"[Date - Publication] : \"$years[$y]\"[Date - Publication]&retmax=$max_ID_search&retstart=0";
                }
                else {
                    $pubmed_search_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=(\"$affiliation_restriction\"[Affiliation]) AND (\"$years[$y]\"[Date - Publication] : \"$years[$y]\"[Date - Publication])&retmax=$max_ID_search&retstart=0";
                }
            }
            else {
                if ($affiliation_restriction eq "none"){
                    $pubmed_search_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=(\"$publication_type\"[Publication Type]) AND (\"$years[$y]\"[Date - Publication] : \"$years[$y]\"[Date - Publication])&retmax=$max_ID_search&retstart=0";
                }
                else {
                    $pubmed_search_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=(\"$affiliation_restriction\"[Affiliation]) AND (\"$publication_type\"[Publication Type]) AND (\"$years[$y]\"[Date - Publication] : \"$years[$y]\"[Date - Publication])&retmax=$max_ID_search&retstart=0";
                }
            }
            my $pubmed_search_file = "$file_directory"."pubmed_result_$years[$y].tab";
            open (my $PSF,'+>',$pubmed_search_file) or die "Cannot open $pubmed_search_file\n";
            my $retrieval = getstore($pubmed_search_url, $pubmed_search_file);
            if (is_error($retrieval)){
                close $PSF;
                unlink $pubmed_search_file;
                print "ERROR: getstore of $pubmed_search_file at $pubmed_search_url failed with $retrieval.\n";
                redo YEAR;
            }
            
            ### Determine the total count
            my $count;
            while (my $line = <$PSF>){
                $line =~ s/^\s+//;
                if ($line =~ m/^<Count>/){
                    $line =~ /<Count>(.*?)<\/Count>/;
                    $count = $1;
                }
                last if (defined $count);
            }
            
            ### Add the first $max_ID_search (or fewer) IDs to an array
            my @ID_list;
            while (my $line = <$PSF>){
                $line =~ s/^\s+//;
                if ($line =~ m/<Id>/){
                    $line =~ /<Id>(.*?)<\/Id>/;
                    my $line_ID = $1;
                    push @ID_list, $line_ID;
                }
                last if ($line =~ m/<\/IdList>/);
            }
            
            ### Add the remaining IDs to the list, if $count is greater than the set maximum
            my $search_max = int(($count/$max_ID_search)+1); #round up to the next integer for max numer of searches
            my $search_num = 1;
            
         PMID: until (($count <= scalar(@ID_list)) or ($search_num >= $search_max)){
                close $PSF;
                unlink $pubmed_search_file;
                ++$search_num;
                #print "\t      $years[$y] Pubmed search $search_num of $search_max.\n"; #test line
                ### Get XML file for the next $max_ID_search (or fewer) pubmed IDs
                if ($publication_type eq "all"){
                    if ($affiliation_restriction eq "none"){
                        $pubmed_search_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=\"$years[$y]\"[Date - Publication] : \"$years[$y]\"[Date - Publication]&retmax=$max_ID_search&retstart=".(scalar(@ID_list));
                    }
                    else {
                        $pubmed_search_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=(\"$affiliation_restriction\"[Affiliation]) AND (\"$years[$y]\"[Date - Publication] : \"$years[$y]\"[Date - Publication])&retmax=$max_ID_search&retstart=".(scalar(@ID_list));
                    }
                }
                else{
                    if ($affiliation_restriction eq "none"){
                        $pubmed_search_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=(\"$publication_type\"[Publication Type]) AND (\"$years[$y]\"[Date - Publication] : \"$years[$y]\"[Date - Publication])&retmax=$max_ID_search&retstart=".(scalar(@ID_list));
                    }
                    else {
                        $pubmed_search_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=(\"$affiliation_restriction\"[Affiliation]) AND (\"$publication_type\"[Publication Type]) AND (\"$years[$y]\"[Date - Publication] : \"$years[$y]\"[Date - Publication])&retmax=$max_ID_search&retstart=".(scalar(@ID_list));
                    }
                }
                $pubmed_search_file = "$file_directory"."pubmed_result_$years[$y].tab";
                open (my $PSF,'+>',$pubmed_search_file) or die "Cannot open $pubmed_search_file\n";
                $retrieval = getstore($pubmed_search_url, $pubmed_search_file);
                if (is_error($retrieval)){
                    close $PSF;
                    unlink $pubmed_search_file;
                    --$search_num;
                    print "ERROR: getstore of $pubmed_search_file at $pubmed_search_url failed with $retrieval.\n";
                    redo PMID;
                }
                
                ### Add the next $max_ID_search (or fewer) IDs to an array
                while (my $line = <$PSF>){
                    $line =~ s/^\s+//;
                    if ($line =~ m/<Id>/){
                        $line =~ /<Id>(.*?)<\/Id>/;
                        my $line_ID = $1;
                        push @ID_list, $line_ID;
                    }
                    last if ($line =~ m/<\/IdList>/);
                }
            }
            #print "Year: $years[$y]\tPubs: $count\tIDs: ",scalar(@ID_list),"\tPM Max: $search_max\tPM Num: ",$search_num,"\n"; #test line
            close $PSF;
            unlink $pubmed_search_file;
            
            ### Gather data from Pubmed according to selected $data_retrieval_modes
            my @csv_authors;
            my @csv_pub_authors;
            my @medline_authors;
            my @medline_authors_short;
            my @medline_firstau_short;
            my @medline_pub_authors;
            my @journals_short;
            my @journals_full;
            my @pmid_errors;
            my @keyword_list;
            
            ### Get necessary number of CSV files
            my $CSV_search_max = int(($count/$max_CSV_search)+1); #round up to the next integer for max numer of searches
            if (($data_retrieval_mode eq 'medline')){$CSV_search_max = 0;}
            my $CSV_search_num = 0;
            my @CSV_file_check;
            
         CSV: until ($CSV_search_num >= $CSV_search_max){ # if CSV_search_max = 0, this loop should skip
                ++$CSV_search_num;
                #print "\t      $years[$y] CSV search $CSV_search_num of $CSV_search_max.\n"; #test line
                
                ### Get CSV file for $max_CSV_search (or fewer) pubmed IDs
                my $search_list;
                if (($CSV_search_num*$max_CSV_search)-1 <= $count){
                    $search_list = join ",", @ID_list[($CSV_search_num-1)*$max_CSV_search..(($CSV_search_num*$max_CSV_search)-1)];
                    $CSV_file_check[$CSV_search_num] = $max_CSV_search;
                }
                else{
                    $search_list = join ",", @ID_list[($CSV_search_num-1)*$max_CSV_search..($count-1)];
                    $CSV_file_check[$CSV_search_num] = $count-(($CSV_search_num-1)*$max_CSV_search);
                }
                my $CSV_search_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=$search_list&rettype=csv&retmode=text";
                #print "$CSV_search_url\n"; #test print line
                my $CSV_search_file = "$file_directory"."pubmed_result_$years[$y]-$CSV_search_num.csv";
                
                ### Get files if they don't already exist
                if (!-e $CSV_search_file){
                    open (my $CSV_write,'>',$CSV_search_file) or die "Cannot open $CSV_search_file\n";
                    $retrieval = getstore($CSV_search_url, $CSV_search_file);
                    if (is_error($retrieval)){
                        close $CSV_write;
                        unlink $CSV_search_file;
                        --$CSV_search_num;
                        print "ERROR: getstore of $CSV_search_file failed with $retrieval.\n";
                        redo CSV;
                    }
                    else{close $CSV_write;}
                }
            }
            
            ### Get necessary number of Medline files
            my $medline_search_max = int(($count/$max_medline_search)+1); #round up to the next integer for max numer of searches
            if (($data_retrieval_mode eq 'csv')){$medline_search_max = 0;}
            my $medline_search_num = 0;
            my @medline_file_check;
            
         MEDLINE: until ($medline_search_num >= $medline_search_max){ # if $medline_search_max = 0, this loop should skip
                ++$medline_search_num;
                #print "\t      $years[$y] medline search $medline_search_num of $medline_search_max.\n"; #test line
                
                ### Get Medline file for $max_medline_search (or fewer) pubmed IDs
                my $search_list;
                if (($medline_search_num*$max_medline_search)-1 <= $count){
                    $search_list = join ",", @ID_list[($medline_search_num-1)*$max_medline_search..(($medline_search_num*$max_medline_search)-1)];
                    $medline_file_check[$medline_search_num] = $max_medline_search;
                }
                else{
                    $search_list = join ",", @ID_list[($medline_search_num-1)*$max_medline_search..($count-1)];
                    $medline_file_check[$medline_search_num] = $count-(($medline_search_num-1)*$max_medline_search);
                }
                my $medline_search_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=$search_list&rettype=medline&retmode=text";
                #print "$medline_search_url\n"; #test print line
                my $medline_search_file = "$file_directory"."pubmed_result_$years[$y]-$medline_search_num.tab";
                
                ### Get files if they don't already exist
                if (!-e $medline_search_file){
                    open (my $medline_write,'>',$medline_search_file) or die "Cannot open $medline_search_file\n";
                    $retrieval = getstore($medline_search_url, $medline_search_file);
                    if (is_error($retrieval)){
                        close $medline_write;
                        unlink $medline_search_file;
                        --$medline_search_num;
                        print "ERROR: getstore of $medline_search_file failed with $retrieval.\n";
                        redo MEDLINE;
                    }
                    else{close $medline_write;}
                }
            }
            
            ### Add data from CSV files to necessary arrays
            foreach my $i(1...$CSV_search_max){  # if $CSV_search_max = 0 this loop should skip;
                my $CSV_search_file = "$file_directory"."pubmed_result_$years[$y]-$i.csv";
                
                ### Add details from next $max_CSV_search (or fewer) pubmed IDs to necessary arrays
                if (-e $CSV_search_file){
                    open (my $CSV,'<',$CSV_search_file) or die "Cannot open $CSV_search_file\n";
                    my $file_check_counter = 0;
                    while (my $line = <$CSV>){
                        unless ($line =~ m/^Title/){
                            ++$file_check_counter;
                            $line = lc $line;
                            my @line_authors;
                            my @line_sections = split /","/,$line;
                            if ((defined $line_sections[2]) and ($line_sections[2] !~ m/^\\/) and ($line_sections[2] !~ m/^\[/) and ($line_sections[2] !~ m/et al/)
                                and ($line_sections[2] !~ m/cente(rs for|r for)/) and ($line_sections[2] !~ m/united states/) and ($line_sections[2] !~ m/united nations/)
                                and ($line_sections[2] !~ m/\(/)) {
                                @line_authors = split /,/, $line_sections[2];
                            }
                            push @csv_authors, @line_authors;
                            if (scalar(@line_authors) > 0) {
                                push @csv_pub_authors, scalar(@line_authors);
                            }
                        }
                    }
                    if ($file_check_counter < $CSV_file_check[$i]){
                        print "ERROR: $CSV_search_file should have $CSV_file_check[$i] entries, but only has $file_check_counter.\n";
                        close $CSV;
                        unlink $CSV_search_file;
                        $pm->finish; # do the exit in the child process
                        redo YEAR;
                    }
                    else{
                        #print "\tCSV file $i has $file_check_counter out of $CSV_file_check[$i] entries.\n"; #test line
                        $csv_pmid_check = $csv_pmid_check + $file_check_counter;
                        foreach my $k (@csv_authors){
                            $k =~ s/\.+//g;
                            $k =~ s/;.*//;
                            $k =~ s/^\s+|\s+$//;
                        }
                        close $CSV;
                    }
                }
                else{
                    $pm->finish; # do the exit in the child process
                    redo YEAR;
                }
            }
         
            ### Add details from Medline files to necessary arrays   
            foreach my $i (1..$medline_search_max){ # if medline_search_max = 0, this loop should skip
                my $medline_search_file = "$file_directory"."pubmed_result_$years[$y]-$i.tab";
                if (-e $medline_search_file){
                    #my $pmid_found_number = 0;
                    #my @pmid_search_list = split /,/,$search_list;
                    #my $pmid_search_number = scalar(@pmid_search_list);
                    my $file_check_counter = 0;
                    my @medline_file = get_file_content($medline_search_file);
                    foreach my $line (0..$#medline_file){
                        if ($medline_file[$line] =~ m/^PMID/){
                            ++$file_check_counter;
                            #++$pmid_found_number;
                            my $next_line = 1;
                            my $pub_author_num = 0;
                            my $first_shortauthor_num = 0;
                            until (($medline_file[$line+$next_line] =~/^PMID/) or ($line+$next_line == $#medline_file)){
                                if ($medline_file[$line+$next_line] =~/^FAU/){
                                    my $line_var = $medline_file[$line+$next_line];
                                    $line_var = lc $line_var;
                                    $line_var =~ s/^fau - //g;
                                    $line_var =~ s/\.+//g;
                                    #$line_var =~ s/;.*//;
                                    $line_var =~ s/^\s+|\s+$//;
                                    push @medline_authors, $line_var;
                                    ++$pub_author_num;
                                }
                                elsif ($medline_file[$line+$next_line] =~/^AU/){
                                    my $line_var = $medline_file[$line+$next_line];
                                    $line_var = lc $line_var;
                                    $line_var =~ s/^au  - //g;
                                    $line_var =~ s/\.+//g;
                                    #$line_var =~ s/;.*//;
                                    $line_var =~ s/^\s+|\s+$//;
                                    push @medline_authors_short, $line_var;
                                    if ($first_shortauthor_num == 0){
                                        push @medline_firstau_short, $line_var;
                                    }
                                    $first_shortauthor_num = 1;
                                }
                                elsif ($medline_file[$line+$next_line] =~/^MH  - /){
                                    my $line_var = $medline_file[$line+$next_line];
                                    $line_var = lc $line_var;
                                    $line_var =~ s/^mh  - //g;
                                    $line_var =~ s/\.+//g;
                                    $line_var =~ s/\*//g;
                                    $line_var =~ s/^\s+|\s+$//;
                                    push @keyword_list, $line_var;
                                }
                                elsif ($medline_file[$line+$next_line] =~/^TA/){
                                    my $line_var = $medline_file[$line+$next_line];
                                    $line_var = lc $line_var;
                                    $line_var =~ s/^ta  - //g;
                                    $line_var =~ s/\.+//g;
                                    #$line_var =~ s/;.*//;
                                    $line_var =~ s/^\s+|\s+$//;
                                    push @journals_short, $line_var;
                                }
                                elsif ($medline_file[$line+$next_line] =~/^JT/){
                                    my $line_var = $medline_file[$line+$next_line];
                                    $line_var = lc $line_var;
                                    $line_var =~ s/^jt  - //g;
                                    $line_var =~ s/\.+//g;
                                    #$line_var =~ s/;.*//;
                                    $line_var =~ s/^\s+|\s+$//;
                                    push @journals_full, $line_var;
                                }
                                ++$next_line;
                            }
                            if ($pub_author_num > 0){
                                push @medline_pub_authors, $pub_author_num;
                            }
                        }
                        elsif (($medline_file[$line] =~ m/^id:/) and ($medline_file[$line] =~ m/Error occurred: The following PMID is not available:/)){
                            ++$file_check_counter;
                            print "ERROR: $medline_search_file has a missing PMID:\n\t$medline_file[$line]\n";
                            push @pmid_errors, "ERROR: $medline_search_file has a missing PMID: $medline_file[$line]";
                        }
                    }
                    if ($file_check_counter < $medline_file_check[$i]){
                        print "ERROR: $medline_search_file should have $medline_file_check[$i] entries, but only has $file_check_counter.\n";
                        unlink $medline_search_file;
                        $pm->finish; # do the exit in the child process
                        redo YEAR;
                    }
                    else{
                        #print "\tMedline file $i has $file_check_counter out of $medline_file_check[$i] entries.\n"; #test line
                        $medline_pmid_check = $medline_pmid_check + $file_check_counter;
                    }
                    #print "PMID search num:\t$pmid_search_number\n"; #test
                    #print "PMID found num:\t$pmid_found_number\n"; #test
                    #print "Authors\n@medline_authors\n\n"; #test
                    #print "Pub Author Nums\n@medline_pub_authors\n\n"; #test
                }
                else{
                    $pm->finish; # do the exit in the child process
                    redo YEAR;
                }
            }
            
            ### Get retraction count for the given year (retracted papers originally published in the given year)
            my $retraction_url;
            if ($affiliation_restriction eq "none"){
                $retraction_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=(\"retracted publication\"[Publication Type]) AND (\"$years[$y]\"[Date - Publication] : \"$years[$y]\"[Date - Publication])&rettype=count";
            }
            else {
                $retraction_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=(\"retracted publication\"[Publication Type]) AND (\"$affiliation_restriction\"[Affiliation]) AND (\"$years[$y]\"[Date - Publication] : \"$years[$y]\"[Date - Publication])&rettype=count";
            }
            my $retraction_search_file = "$file_directory"."retraction_result_$years[$y].tab";
            open (my $RET,'+>',$retraction_search_file) or die "Cannot open $retraction_search_file\n";
            $retrieval = getstore($retraction_url, $retraction_search_file);
            if (is_error($retrieval)){
                close $RET;
                unlink $retraction_search_file;
                print "ERROR: getstore of $retraction_search_file at $retraction_url failed with $retrieval.\n";
                redo YEAR;
            }
            my $retraction_count;
            while (my $line = <$RET>){
                $line =~ s/^\s+//;
                if ($line =~ m/^<Count>/){
                    $line =~ /<Count>(.*?)<\/Count>/;
                    $retraction_count = $1;
                }
                last if (defined $retraction_count);
            }
            close $RET;
            unlink $retraction_search_file;
            
            ### Get count of retractions made in the given year (originally published in any prior year)
            my $made_retraction_url;
            if ($affiliation_restriction eq "none"){
                $made_retraction_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=(\"retraction of publication\"[Publication Type]) AND (\"$years[$y]\"[Date - Publication] : \"$years[$y]\"[Date - Publication])&rettype=count";
            }
            else {
                $made_retraction_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=(\"retraction of publication\"[Publication Type]) AND (\"$affiliation_restriction\"[Affiliation]) AND (\"$years[$y]\"[Date - Publication] : \"$years[$y]\"[Date - Publication])&rettype=count";
            }
            my $made_retraction_search_file = "$file_directory"."made_retraction_result_$years[$y].tab";
            open (my $madeRET,'+>',$made_retraction_search_file) or die "Cannot open $made_retraction_search_file\n";
            $retrieval = getstore($made_retraction_url, $made_retraction_search_file);
            if (is_error($retrieval)){
                close $madeRET;
                unlink $made_retraction_search_file;
                print "ERROR: getstore of $made_retraction_search_file at $made_retraction_url failed with $retrieval.\n";
                redo YEAR;
            }
            my $made_retraction_count;
            while (my $line = <$madeRET>){
                $line =~ s/^\s+//;
                if ($line =~ m/^<Count>/){
                    $line =~ /<Count>(.*?)<\/Count>/;
                    $made_retraction_count = $1;
                }
                last if (defined $made_retraction_count);
            }
            close $madeRET;
            unlink $made_retraction_search_file;
            
            ###Find the average, max, min, and variance authors per publication in the given year
            my $csv_pub_author_sum = 0;
            my $csv_authors_per_pub;
            my $medline_pub_author_sum = 0;
            my $medline_authors_per_pub;
            my $medline_app_max;
            my $medline_app_min;
            my $medline_app_var;
            my $csv_app_max;
            my $csv_app_min;
            my $csv_app_var;
            foreach my $k (@csv_pub_authors){
                $csv_pub_author_sum = $k + $csv_pub_author_sum;
            }
            foreach my $k (@medline_pub_authors){
                $medline_pub_author_sum = $k + $medline_pub_author_sum;
            }
            if (scalar(@csv_pub_authors) > 0){
                ($csv_authors_per_pub, $csv_app_max, $csv_app_min, $csv_app_var) = array_stats(\@csv_pub_authors);
            }
            else{
                $csv_app_var = 0;
                $csv_app_max = 0;
                $csv_app_min = 0;
                $csv_authors_per_pub = 0;
            }
            if (scalar(@medline_pub_authors) > 0){
                ($medline_authors_per_pub, $medline_app_max, $medline_app_min, $medline_app_var) = array_stats(\@medline_pub_authors);
            }
            else{
                $medline_app_var = 0;
                $medline_app_max = 0;
                $medline_app_min = 0;
                $medline_authors_per_pub = 0;
            }
            #print "CSV stats: Array @csv_pub_authors\n\t\tMean $csv_authors_per_pub, Max $csv_app_max, Min $csv_app_min, Var $csv_app_var\n"; #test
            #print "Medline stats: Array @medline_pub_authors\n\t\tMean $medline_authors_per_pub, Max $medline_app_max, Min $medline_app_min, Var $medline_app_var\n"; #test
            
            ### When all the data has been gathered, print the status and remove the search files
            print "\tYear: $years[$y]\tPubs: $count\tIDs: ",scalar(@ID_list),"\tCSV Author Pubs: ",scalar(@csv_pub_authors),"\tMedline Author Pubs: ",scalar(@medline_pub_authors),"\n"; #test line
            print "\t\t\tPM Max: $search_max\tPM Num: $search_num\tCSV Max: $CSV_search_max\tCSV Num: $CSV_search_num\tMedline Max: $medline_search_max\tMedline Num: $medline_search_num\n"; #test line
            print "\t\t\tCSV IDs: $csv_pmid_check\tMedline IDs: $medline_pmid_check\n";
            if ($CSV_search_num > 0){
                foreach my $k (0..$CSV_search_num){
                    my $CSV_search_file = "$file_directory"."pubmed_result_$years[$y]-$k.csv";
                    unlink $CSV_search_file;
                }
            }
            if ($medline_search_num > 0){
                foreach my $k (0..$medline_search_num){
                    my $medline_search_file = "$file_directory"."pubmed_result_$years[$y]-$k.tab";
                    unlink $medline_search_file;
                }
            }
            
            ### Separate duplicate authors and find frequencies
            my ($unique_array_ref, $duplicate_array_ref, $frequency_ref) = unique_array(\@csv_authors);
            my @csv_unique_authors = @$unique_array_ref;
            my @csv_duplicate_authors = @$duplicate_array_ref;
            my %csv_author_frequency = %$frequency_ref;
            my $csv_total_duplicates = 0;
            foreach my $k (@csv_duplicate_authors){
                $csv_total_duplicates = $csv_total_duplicates + $csv_author_frequency{$k};
            }
            #print "CSV Author num: ",scalar(@csv_authors),"\n@csv_authors\n\n"; #test print
            #print "CSV Unique author num: ",scalar(@csv_unique_authors),"\n@csv_unique_authors\n\n"; #test print
            #print "CSV Duplicate author num: ",scalar(@csv_duplicate_authors),"\n@csv_duplicate_authors\n\n"; #test print
            
            ($unique_array_ref, $duplicate_array_ref, $frequency_ref) = unique_array(\@medline_authors);
            my @medline_unique_authors = @$unique_array_ref;
            my @medline_duplicate_authors = @$duplicate_array_ref;
            my %medline_author_frequency = %$frequency_ref;
            my $medline_total_duplicates = 0;
            foreach my $k (@medline_duplicate_authors){
                $medline_total_duplicates = $medline_total_duplicates + $medline_author_frequency{$k};
            }
            #print "Medline Author num: ",scalar(@medline_authors),"\n@medline_authors\n\n"; #test print
            #print "Medline Unique author num: ",scalar(@medline_unique_authors),"\n@medline_unique_authors\n\n"; #test print
            #print "Medline Duplicate author num: ",scalar(@medline_duplicate_authors),"\n@medline_duplicate_authors\n\n"; #test print
            
            ($unique_array_ref, $duplicate_array_ref, $frequency_ref) = unique_array(\@medline_authors_short);
            my @medline_unique_shortauthors = @$unique_array_ref;
            my @medline_duplicate_shortauthors = @$duplicate_array_ref;
            my %medline_shortauthor_frequency = %$frequency_ref;
            #print "Medline Short Author num: ",scalar(@medline_authors_short),"\n@medline_authors_short\n\n"; #test print
            #print "Medline Unique Short author num: ",scalar(@medline_unique_shortauthors),"\n@medline_unique_shortauthors\n\n"; #test print
            #print "Medline Duplicate Short author num: ",scalar(@medline_duplicate_shortauthors),"\n@medline_duplicate_shortauthors\n\n"; #test print
            
            ($unique_array_ref, $duplicate_array_ref, $frequency_ref) = unique_array(\@medline_firstau_short);
            my @medline_unique_shortfirstau = @$unique_array_ref;
            my @medline_duplicate_shortfirstau = @$duplicate_array_ref;
            my %medline_shortfirstau_frequency = %$frequency_ref;
            
            ### Separate duplicate journal titles and find frequencies
            ($unique_array_ref, $duplicate_array_ref, $frequency_ref) = unique_array(\@journals_full);
            my @journals_unique = @$unique_array_ref;
            my @journals_duplicate = @$duplicate_array_ref;
            my %journals_frequency = %$frequency_ref;
            ($unique_array_ref, $duplicate_array_ref, $frequency_ref) = unique_array(\@journals_short);
            my @journals_unique_short = @$unique_array_ref;
            #print "PMIDs with Journal Titles: ",scalar(@journals_short),"\n@journals_short\n\n"; #test print
            #print "Unique Journal Titles: ",scalar(@journals_unique),"\n\n"; #test print
            #foreach my $k (@journals_unique){print "$journals_frequency{$k}\t$k\n";}  #test print
            
            ### Separate duplicate keywords and find frequencies
            ($unique_array_ref, $duplicate_array_ref, $frequency_ref) = unique_array(\@keyword_list);
            my @keywords_unique = @$unique_array_ref;
            my %keywords_frequency = %$frequency_ref;
            #print "Unique Keywords: ",scalar(@keywords_unique),"\n\n"; #test print
            #foreach my $k (@keywords_unique){print "$keywords_frequency{$k}\t$k\n";}  #test print
            
            ### Create output files, and print
            open (my $output,'>',$output_filename) or die "Cannot open $output_filename\n";
            if ($count > 0){
                
                ###Print Header
                print $output "Date: $date\t Scientist Search $version\n";
                print $output "Data Retrieval Mode: $data_retrieval_mode, reflects file type source of PubMed data in following analysis.\n";
                print $output "Publication Type Selection: $publication_type, Affiliation Restriction: $affiliation_restriction, reflects search constraints of PubMed data in following analysis.\n";
                print $output "Keyword Mode Selection: $keyword_mode\n\n";
                print $output "Year: $years[$y]\tPubs: $count\tIDs: ",scalar(@ID_list),"\tCSV Author Pubs: ",scalar(@csv_pub_authors),"\tMedline Author Pubs: ",scalar(@medline_pub_authors),"\n";
                print $output "\t\tPM Max: $search_max\tPM Num: $search_num\tCSV Max: $CSV_search_max\tCSV Num: $CSV_search_num\tMedline Max: $medline_search_max\tMedline Num: $medline_search_num\n";
                print $output "\t\tCSV IDs: $csv_pmid_check\tMedline IDs: $medline_pmid_check\n\n";
                
                ###Print Error List
                if (scalar (@pmid_errors) > 0){
                    print $output "List of errors:\n";
                    foreach my $k (@pmid_errors){print $output "$k\n";}
                    print $output "\n";
                }
                
                ###Print Data
                print $output "Publication PubMed Count:<$count>\nCSV Publications with Authors:<",scalar(@csv_pub_authors),">\nMedline Publications with Authors:<",scalar(@medline_pub_authors),">\n";
                print $output "Retracted Publications:<$retraction_count>\nRetractions Made:<$made_retraction_count>\nCSV Unique authors:<",scalar(@csv_unique_authors),">\nMedline Unique authors:<",scalar(@medline_unique_authors),">\nMedline Unique short authors:<",scalar(@medline_unique_shortauthors),">\nMedline Unique short first authors:<",scalar(@medline_unique_shortfirstau),">\n";
                print $output "CSV Duplicate authors:<",scalar(@csv_duplicate_authors),">\nMedline Duplicate authors:<",scalar(@medline_duplicate_authors),">\nMedline Duplicate short authors:<",scalar(@medline_duplicate_shortauthors),">\n";
                print $output "CSV Duplicate author publications:<$csv_total_duplicates>\n";
                print $output "Medline Duplicate author publications:<$medline_total_duplicates>\n";
                print $output "CSV Authors per publication:<$csv_authors_per_pub>\nCSV APP max:<$csv_app_max>\nCSV APP min:<$csv_app_min>\nCSV APP var:<$csv_app_var>\n";
                print $output "Medline Authors per publication:<$medline_authors_per_pub>\nMedline APP max:<$medline_app_max>\nMedline APP min:<$medline_app_min>\nMedline APP var:<$medline_app_var>\n";
                print $output "Unique Journal Titles:<",scalar(@journals_unique),">\n";
                print $output "Unique Short Journal Titles:<",scalar(@journals_unique_short),">\n";
                unless (scalar(@csv_unique_authors) == 0){
                    print $output "CSV Publications per unique author:\t",$count/scalar(@csv_unique_authors),"\n";
                }
                unless (scalar(@csv_duplicate_authors) == 0){
                    print $output "CSV Duplicate publications per duplicate author:\t",$csv_total_duplicates/scalar(@csv_duplicate_authors),"\n";
                }
                unless (scalar(@csv_unique_authors) == 0){
                    print $output "CSV Duplicate authors / Unique authors\t",scalar(@csv_duplicate_authors)/scalar(@csv_unique_authors),"\n";
                }
                unless (scalar(@medline_unique_authors) == 0){
                    print $output "Medline Publications per unique author:\t",$count/scalar(@medline_unique_authors),"\n";
                }
                unless (scalar(@medline_duplicate_authors) == 0){
                    print $output "Medline Duplicate publications per duplicate author:\t",$medline_total_duplicates/scalar(@medline_duplicate_authors),"\n";
                }
                unless (scalar(@medline_unique_authors) == 0){
                    print $output "Medline Duplicate authors / Unique authors\t",scalar(@medline_duplicate_authors)/scalar(@medline_unique_authors),"\n\n";
                }
                print $output "\n*****JOURNAL PUB FREQUENCIES*****\n";
                foreach my $k (sort {$journals_frequency{$b} <=> $journals_frequency{$a}} keys %journals_frequency){
                    print $output "$journals_frequency{$k}\t$k\n";
                }
                
                ### Print Keyword frequences for the year's publications based on $keyword_mode
                print $output "\n*****KEYWORDS LIST*****\n";
                if ($keyword_mode eq "all"){
                    print $output "\n*****KEYWORDS FREQUENCIES*****\n<";
                    foreach my $k (sort {$keywords_frequency{$b} <=> $keywords_frequency{$a}} keys %keywords_frequency){
                        print $output "$keywords_frequency{$k}*$k|";
                    }
                    print $output ">\n";
                }
                elsif ($keyword_mode eq "query"){
                    foreach my $k (@keyword_list){
                        print $output "$k,";
                    }
                    print $output "\n\n*****KEYWORDS FREQUENCIES*****\n<";
                    foreach my $k (sort {$keywords_frequency{$b} <=> $keywords_frequency{$a}} keys %keywords_frequency){
                        foreach my $j (@keyword_list){
                            if ($j eq $k){
                                print $output "$keywords_frequency{$k}*$k|";
                            }
                        }
                        print $output ">\n";
                    }
                }
                else {
                    print $output "\n*****KEYWORDS FREQUENCIES*****\n<>\n";
                }
                
                print $output "\n*****CSV PUBMED AUTHOR FREQUENCIES*****\n";
                foreach my $k (sort {$csv_author_frequency{$b} <=> $csv_author_frequency{$a}} keys %csv_author_frequency){
                    print $output "$csv_author_frequency{$k}\t$k\n";
                }
                print $output "\n*****CSV PUBMED PUB FREQUENCIES*****\n<";
                foreach my $k (@csv_pub_authors){
                    print $output "$k,";
                }
                print $output ">\n";
                print $output "\n*****MEDLINE PUBMED AUTHOR FREQUENCIES*****\n";
                foreach my $k (sort {$medline_author_frequency{$b} <=> $medline_author_frequency{$a}} keys %medline_author_frequency){
                    print $output "$medline_author_frequency{$k}\t$k\n";
                }
                print $output "\n*****MEDLINE PUBMED SHORT AUTHOR FREQUENCIES*****\n";
                foreach my $k (sort {$medline_shortauthor_frequency{$b} <=> $medline_shortauthor_frequency{$a}} keys %medline_shortauthor_frequency){
                    print $output "$medline_shortauthor_frequency{$k}\t$k\n";
                }
                print $output "\n*****MEDLINE PUBMED PUB FREQUENCIES*****\n<";
                foreach my $k (@medline_pub_authors){
                    print $output "$k,";
                }
                print $output ">\n";
                close $output;
            }
            else {
                
                ###Print Header
                print $output "Date: $date\t Scientist Search $version\n";
                print $output "Data Retrieval Mode: $data_retrieval_mode, reflects file type source of PubMed data in following analysis.\n";
                print $output "Publication Type Selection: $publication_type, Affiliation Restriction: $affiliation_restriction, reflects search constraints of PubMed data in following analysis.\n";                print $output "Keyword Mode Selection: $keyword_mode\n\n";
                print $output "Year: $years[$y]\tPubs: $count\tIDs: ",scalar(@ID_list),"\tCSV Author Pubs: ",scalar(@csv_pub_authors),"\tMedline Author Pubs: ",scalar(@medline_pub_authors),"\n";
                print $output "\t\tPM Max: $search_max\tPM Num: $search_num\tCSV Max: $CSV_search_max\tCSV Num: $CSV_search_num\tMedline Max: $medline_search_max\tMedline Num: $medline_search_num\n";
                print $output "\t\tCSV IDs: $csv_pmid_check\tMedline IDs: $medline_pmid_check\n\n";
                
                ###Print Error List
                if (scalar (@pmid_errors) > 0){
                    print $output "List of errors:\n";
                    foreach my $k (@pmid_errors){print $output "$k\n";}
                    print $output "\n";
                }
                
                ###Print Data
                print $output "Publication PubMed Count:<$count>\nCSV Publications with Authors:<",scalar(@csv_pub_authors),">\nMedline Publications with Authors:<",scalar(@medline_pub_authors),">\n";
                print $output "Retracted Publications:<$retraction_count>\nRetractions Made:<$made_retraction_count>\nCSV Unique authors:<",scalar(@csv_unique_authors),">\nMedline Unique authors:<",scalar(@medline_unique_authors),">\nMedline Unique short authors:<",scalar(@medline_unique_shortauthors),">\nMedline Unique short first authors:<",scalar(@medline_unique_shortfirstau),">\n";
                print $output "CSV Duplicate authors:<",scalar(@csv_duplicate_authors),">\nMedline Duplicate authors:<",scalar(@medline_duplicate_authors),">\nMedline Duplicate short authors:<",scalar(@medline_duplicate_shortauthors),">\n";
                print $output "CSV Duplicate author publications:<$csv_total_duplicates>\n";
                print $output "Medline Duplicate author publications:<$medline_total_duplicates>\n";
                print $output "CSV Authors per publication:<$csv_authors_per_pub>\nCSV APP max:<$csv_app_max>\nCSV APP min:<$csv_app_min>\nCSV APP var:<$csv_app_var>\n";
                print $output "Medline Authors per publication:<$medline_authors_per_pub>\nMedline APP max:<$medline_app_max>\nMedline APP min:<$medline_app_min>\nMedline APP var:<$medline_app_var>\n";
                print $output "Unique Journal Titles:<",scalar(@journals_unique),">\n";
                print $output "Unique Short Journal Titles:<",scalar(@journals_unique_short),">\n";
                close $output;
            }
            if (($data_retrieval_mode eq 'all') or ($data_retrieval_mode eq 'medline')){
                if ($medline_pmid_check < $count) {
                    print "ERROR: $years[$y] should have $count medline PMIDs, but only has $medline_pmid_check.\n";
                    unlink $output_filename;
                    $pm->finish; # do the exit in the child process
                    redo YEAR;
                }
            }
            if (($data_retrieval_mode eq 'all') or ($data_retrieval_mode eq 'csv')){
                if ($csv_pmid_check < $count) {
                    print "ERROR: $years[$y] should have $count CSV PMIDs, but only has $csv_pmid_check.\n";
                    unlink $output_filename;
                    $pm->finish; # do the exit in the child process
                    redo YEAR;
                }
            }
        }
        elsif (-e $output_filename){print"\tYear: $years[$y] has already been completed.\n";}
        $pm->finish; # do the exit in the child process
    }
    $pm->wait_all_children;
    #print "starting data compilation\n"; #test line
    
    
=for
    ### Edit yearly result files if necessary before tabulating ss_data file.
    foreach my $y (0..$#years){
        my $authored_publications;
        my $temp_file = "$file_directory"."$years[$y]temp.tab";
        my $output_filename = "$file_directory"."$years[$y]"."_pubmed_authors.tab";
        open (my $output,'<',$output_filename) or die "Cannot open $output_filename\n";
        open (my $temp, '>', $temp_file) or die "Cannot open $temp_file\n";
        while (my $line = <$output>){
            if ($line =~ m/^Year:/){
                $line =~ /Author Pubs: (.*?)\t/;
                my $line_num = $1;
                $authored_publications = $line_num;
            }
            if ($line =~ m/^Publication number:/){
                $line =~s/Publication number/Publication PubMed Count/;
                $line =~s/>/>\nPublications with Authors:<$authored_publications>/;
            }
            print $temp "$line";
        }
        close $output;
        close $temp;
        rename $temp_file,$output_filename;
    }
=cut
    
    # Fill hashes with data from each year
    my %pubmed_count;
    my %csv_authored_publication_num;
    my %csv_unique_author_num;
    my %csv_duplicate_author_num;
    my %csv_duplicate_publication_num;
    my %csv_authors_per_pub;
    my %csv_app_max;
    my %csv_app_min;
    my %csv_app_var;
    my %medline_authored_publication_num;
    my %medline_unique_author_num;
    my %medline_unique_shortauthor_num;
    my %medline_unique_shortfirstau_num;
    my %medline_duplicate_author_num;
    my %medline_duplicate_shortauthor_num;
    my %medline_duplicate_publication_num;
    my %medline_authors_per_pub;
    my %medline_app_max;
    my %medline_app_min;
    my %medline_app_var;
    my %retraction_num;
    my %made_ret_num;
    my %medline_pub_list;
    my %csv_pub_list;
    my %unique_journals;
    my %unique_journals_short;
    my %journal_publication_num;
    my %journal_publication_total;
    my @comp_journal_list;
    my %yearly_journal_data;
    my $journal_impact_header;
    my %journal_if_matchlist;
    #my @comp_keyword_list;
    my %yearly_keyword_data;
    my %keyword_numbers;
    
    foreach my $y (0..$#years){
        
        ### Open the completed output file and read data into output hashes
        #print "starting data1 $years[$y]\n"; #test line
        my $output_filename = "$file_directory"."$years[$y]"."_pubmed_authors.tab";
        open (my $output,'<',$output_filename) or die "Cannot open $output_filename\n";
        while (my $line = <$output>){
            if ($line =~ m/^Publication PubMed Count:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $pubmed_count{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^Retracted Publications:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $retraction_num{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^Retractions Made:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $made_ret_num{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^CSV Publications with Authors:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $csv_authored_publication_num{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^CSV Unique authors:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $csv_unique_author_num{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^CSV Duplicate authors:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $csv_duplicate_author_num{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^CSV Duplicate author publications:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $csv_duplicate_publication_num{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^CSV Authors per publication:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $csv_authors_per_pub{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^CSV APP max:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $csv_app_max{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^CSV APP min:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $csv_app_min{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^CSV APP var:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $csv_app_var{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^Medline Publications with Authors:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $medline_authored_publication_num{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^Medline Unique authors:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $medline_unique_author_num{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^Medline Unique short authors:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $medline_unique_shortauthor_num{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^Medline Unique short first authors:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $medline_unique_shortfirstau_num{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^Medline Duplicate authors:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $medline_duplicate_author_num{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^Medline Duplicate short authors:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $medline_duplicate_shortauthor_num{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^Medline Duplicate author publications:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $medline_duplicate_publication_num{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^Medline Authors per publication:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $medline_authors_per_pub{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^Medline APP max:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $medline_app_max{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^Medline APP min:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $medline_app_min{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^Medline APP var:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $medline_app_var{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^Unique Journal Titles:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $unique_journals{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^Unique Short Journal Titles:/){
                $line =~ /<(.*?)>/;
                my $line_num = $1;
                $unique_journals_short{$years[$y]} = $line_num;
            }
            elsif ($line =~ m/^\*\*\*\*\*MEDLINE PUBMED PUB FREQUENCIES\*\*\*\*\*/){
                my $next_line = <$output>;
                $next_line =~ /<(.*?)>/;
                my $line_num = $1;
                my @medline_pub_list = split /,/,$line_num;
                $medline_pub_list{$years[$y]} = [@medline_pub_list];
            }
            elsif ($line =~ m/^\*\*\*\*\*CSV PUBMED PUB FREQUENCIES\*\*\*\*\*/){
                my $next_line = <$output>;
                $next_line =~ /<(.*?)>/;
                my $line_num = $1;
                my @csv_pub_list = split /,/,$line_num;
                $csv_pub_list{$years[$y]} = [@csv_pub_list];
            }
            elsif ($line =~ m/^\*\*\*\*\*KEYWORDS FREQUENCIES\*\*\*\*\*/){
                my $next_line = <$output>;
                $next_line =~ /<(.*?)>/;
                my $line_num = $1;
                @{$yearly_keyword_data{$years[$y]}} = split /\|/,$line_num;
                foreach my $k (@{$yearly_keyword_data{$years[$y]}}){
                    my @k_vars = split /\*/, $k;
                    $k_vars[0] =~ s/^\s+|\s+$//;
                    $k_vars[1] =~ s/^\s+|\s+$//;
                    $keyword_numbers{$k_vars[1]}[$y] = $k_vars[0];
                }
            }
            
        }
        close $output;
    }
    
    foreach my $k (keys %keyword_numbers){
        foreach my $y (0..$#years){
            if (!defined $keyword_numbers{$k}[$y]){
                $keyword_numbers{$k}[$y] = 0;
            }
        }
    }
    
=for 
    unless ($keyword_mode eq "none"){
        foreach my $y (0..$#years){
            
            ### Open the completed output file and read data into output hashes
            #print "starting dataKEY $years[$y]\n"; #test line
            my $output_filename = "$file_directory"."$years[$y]"."_pubmed_authors.tab";
            open (my $output,'<',$output_filename) or die "Cannot open $output_filename\n";
            my $keyword_remember_var = 0;
            while (my $line = <$output>){
                if ($line =~ m/^\*\*\*\*\*CSV PUBMED AUTHOR FREQUENCIES\*\*\*\*\*/){
                    $keyword_remember_var = 0; #stop pushing subsequent lines to the journal array
                }
                if (($keyword_remember_var ==1) and ($line =~ m/^\d+/)){
                    chomp $line;
                    push @{$yearly_keyword_data{$years[$y]}}, $line;
                    my @line_vars = split /\t/,$line;
                    $line_vars[1] =~ s/^\s+|\s+$//; #remove spaces at the beginning or end of the journal title
                    push @comp_keyword_list, $line_vars[1];
                }
                if ($line =~ m/^\*\*\*\*\*KEYWORDS FREQUENCIES\*\*\*\*\*/){
                    $keyword_remember_var = 1; #start pushing subsequent lines to the journal array
                }
            }
            close $output;
            #print "keywords read for $years[$y]\n"; #test line
        }
        my ($array_unique, $array_duplicates, $array_frequency) = unique_array(\@comp_keyword_list);
        my @unique_keyword_list = @$array_unique;
        foreach my $keyword(@unique_keyword_list){
            foreach my $y (@years){
                my $keyword_match = 0;
                my $yearly_keyword_hits = 0;
                foreach my $i (@{$yearly_keyword_data{$y}}){
                    my @line_vars = split /\*/, $i;
                    $line_vars[1] =~ s/^\s+|\s+$//; #remove spaces at the beginning or end of the keyword
                    if ($line_vars[1] eq $keyword) {
                        ++$keyword_match;
                        $yearly_keyword_hits = $yearly_keyword_hits + $line_vars[0];
                    }
                }
                if ($keyword_match == 0) {push @{$keyword_numbers{$keyword}}, 0;}
                elsif ($keyword_match == 1){push @{$keyword_numbers{$keyword}}, $yearly_keyword_hits;}
                elsif ($keyword_match >1){
                    push @{$keyword_numbers{$keyword}}, $yearly_keyword_hits;
                    print "\tNote: $keyword had $keyword_match matches in $y\n";
                }
            }
            #print "$keyword data found.\n"; #test line
        }
    }
    #print Dumper %yearly_keyword_data; #test line
    #print Dumper %keyword_numbers; #test line
=cut
    
    foreach my $y (0..$#journal_years){
        
        ### Open the completed output file and read data into output hashes
        my $output_filename = "$file_directory"."$journal_years[$y]"."_pubmed_authors.tab";
        open (my $output,'<',$output_filename) or die "Cannot open $output_filename\n";
        my $journal_remember_var = 0;
        while (my $line = <$output>){
            if ($line =~ m/^\*\*\*\*\*KEYWORDS LIST\*\*\*\*\*/){
                $journal_remember_var = 0; #stop pushing subsequent lines to the journal array
            }
            if (($journal_remember_var ==1) and ($line =~ m/^\d+/)){
                chomp $line;
                push @{$yearly_journal_data{$journal_years[$y]}}, $line;
                my @line_vars = split /\t/,$line;
                $line_vars[1] =~ s/\((.*)\)//g; #remove anything within (...) in the journal title
                $line_vars[1] =~ s/\((.*)//g; #remove anything following (... in the journal title
                $line_vars[1] =~ s/^\s+//;
                $line_vars[1] =~ s/^the//; #remove "the" at the beginning of a journal title
                $line_vars[1] =~ s/^\s+//;
                $line_vars[1] =~ s/[^\w\s]//g; #remove non-word characters except spaces
                $line_vars[1] =~ s/"//g;
                $line_vars[1] =~ s/^\s+|\s+$//; #remove spaces at the beginning or end of the journal title
                push @comp_journal_list, $line_vars[1];
            }
            if ($line =~ m/^\*\*\*\*\*JOURNAL PUB FREQUENCIES\*\*\*\*\*/){
                $journal_remember_var = 1; #start pushing subsequent lines to the journal array
            }
        }
        close $output;
    }
    
    ### Get unique list of journal titles for all analyzed years and populate %journal_publication_num hash of arrays with data
    if (scalar(@journal_years) > 0){
        my ($array_unique, $array_duplicates, $array_frequency) = unique_array(\@comp_journal_list);
        my @journal_list = @$array_unique;
        foreach my $journal (@journal_list){
            $journal_publication_total{$journal} = 0;
            foreach my $y (@journal_years){
                my $journal_match = 0;
                my $yearly_journal_pubs = 0;
                foreach my $i (@{$yearly_journal_data{$y}}){
                    my @line_vars = split /\t/, $i;
                    $line_vars[1] =~ s/\((.*)\)//g; #remove anything within (...) in the journal title
                    $line_vars[1] =~ s/\((.*)//g; #remove anything following (... in the journal title
                    $line_vars[1] =~ s/^\s+//;
                    $line_vars[1] =~ s/^the//; #remove "the" at the beginning of a journal title
                    $line_vars[1] =~ s/^\s+//;
                    $line_vars[1] =~ s/[^\w\s]//g; #remove non-word characters except spaces
                    $line_vars[1] =~ s/"//g;
                    $line_vars[1] =~ s/^\s+|\s+$//; #remove spaces at the beginning or end of the journal title
                    if ($journal eq $line_vars[1]){
                        ++$journal_match;
                        $yearly_journal_pubs = $yearly_journal_pubs + $line_vars[0];
                    }
                }
                $journal_publication_total{$journal} = $journal_publication_total{$journal} + $yearly_journal_pubs;
                if ($journal_match == 0){push @{$journal_publication_num{$journal}}, "";}
                elsif ($journal_match == 1){push @{$journal_publication_num{$journal}}, $yearly_journal_pubs;}
                elsif ($journal_match > 1){
                    push @{$journal_publication_num{$journal}}, $yearly_journal_pubs;
                    print "\tNote: $journal had $journal_match matches in $y\n";
                    }
            }
        }
        
        ### Get information on journal impact factors for relevant journals
        my $journal_if_file = "$data_directory"."SNIP_SJR_IPP_complete_1999_2013_v2_Number.csv";
        my @content = get_file_content($journal_if_file);
        my @journal_if_data = split /</, $content[0];
        my %journal_if_data;
        $journal_impact_header = $journal_if_data[1];
        $journal_impact_header =~ s/^,(.*?),>,//g;
        $journal_impact_header =~ s/,/\t/g;
        foreach my $i (2..$#journal_if_data) {
            if ($journal_if_data[$i] =~ m/^,(.*?),>/){
                my $line = $journal_if_data[$i];
                $line =~ /^,(.*?),>/;
                my $line_journal = $1;
                $line_journal = lc $line_journal;
                $line_journal =~ s/\((.*)\)//g;
                $line_journal =~ s/\((.*)//g;
                $line_journal =~ s/^\s+//;
                $line_journal =~ s/^the//; #remove "the" at the beginning of a journal title
                $line_journal =~ s/^\s+//;
                $line_journal =~ s/[^\w\s]//g;
                $line_journal =~ s/"//g;
                $line_journal =~ s/\.+//g;
                $line_journal =~ s/^\s+|\s+$//;
                my $line_vars = $journal_if_data[$i];
                $line_vars =~ s/^,(.*?),>,//g;
                my @line_vars = split /,/,$line_vars;
                if (!defined $journal_if_data{$line_journal}){
                    push @{$journal_if_data{$line_journal}}, @line_vars;
                }
                else {
                    print "ERROR: $line_journal has a duplicate set of impact data:\n\t1:",@{$journal_if_data{$line_journal}},"\n\t2:",@line_vars,"\n";
                    if ((scalar @line_vars) > (scalar @{$journal_if_data{$line_journal}})){
                        $journal_if_data{$line_journal} = @line_vars;
                    }
                    elsif (((scalar @line_vars) == (scalar @{$journal_if_data{$line_journal}})) and ($#line_vars > $#{$journal_if_data{$line_journal}})){
                        $journal_if_data{$line_journal} = @line_vars;
                    }
                }
            }
        }
        #print Dumper (%journal_if_data); #test line
        foreach my $journal (keys(%journal_publication_num)){
            my $match1 = 0;
            foreach my $j_full (keys %journal_if_data){
                if (($j_full eq $journal) and ($match1 == 0)){
                    $journal_if_matchlist{$journal} = $journal_if_data{$j_full};
                }
            }
        }
        #print Dumper %journal_if_matchlist; #test line
    }
    
    ### Create consolidated output files
    my $data_comparison_filename = "$file_directory"."data_modes.tab";
    my $ss_filename = "$file_directory"."ss_data.tab";
    my $app_filename = "$file_directory"."app_data.tab";
    open my $SS, '>', $ss_filename;
    open my $APP, '>', $app_filename;
    open my $MODES, '>', $data_comparison_filename;
    
    ### Print headers
    print $SS "Date: $date\t Scientist Search $version\n";
    print $SS "Data Retrieval Mode: $data_retrieval_mode, reflects file type source of PubMed data in following analysis.\n";
    print $SS "Publication Type Selection: $publication_type, Affiliation Restriction: $affiliation_restriction, reflects search constraints of PubMed data in following analysis.\n";
    print $SS "Keyword Mode Selection: $keyword_mode\n\n";
    print $SS "Year\tPubMed Result Count\tAuthored Publications\tUnique Author No.\tUnique Short Author No.\tDuplicate Author No.\tDuplicate Short Author No.\tDuplicate Author Publications\t";
    print $SS "Retracted Pubs\tRetractions Made\tUnique Journals\tUnique Short Journals\tAuthors Per Publication\tAPP Max\tAPP Min\t APP Var\tUnique First Authors\n";
    
    if ($data_retrieval_mode eq 'all'){
        print $MODES "Date: $date\t Scientist Search $version\n";
        print $MODES "Data Retrieval Mode: $data_retrieval_mode, reflects file type source of PubMed data in following analysis.\n";
        print $MODES "Publication Type Selection: $publication_type, Affiliation Restriction: $affiliation_restriction, reflects search constraints of PubMed data in following analysis.\n";
        print $MODES "Keyword Mode Selection: $keyword_mode\n\n";
        print $MODES "Year\tPubMed Result Count\tCSV Authored Publications\tMedline Authored Publications\tCSV Unique Author No.\tMedline Unique Author No.\tMedline Unique Short Author No.\t";
        print $MODES "CSV Duplicate Author No.\tMedline Duplicate Author No.\tMedline Duplication Short Author No.\tCSV Duplicate Author Publications\tMedline Duplicate Author Publications\t";
        print $MODES "CSV Authors Per Publication\tMedline Authors Per Publication\tCSV APP Max\tMedline APP Max\tCSV APP Min\tMedline APP Min\tCSV APP Var\tMedline APP Var\n";
    }
    
    print $APP "Date: $date\t Scientist Search $version\n";
    print $APP "Data Retrieval Mode: $data_retrieval_mode, reflects file type source of PubMed data in following analysis.\n";
    print $APP "Publication Type Selection: $publication_type, Affiliation Restriction: $affiliation_restriction, reflects search constraints of PubMed data in following analysis.\n";
    print $APP "Keyword Mode Selection: $keyword_mode\n";
    print $APP "Authors per Publication output type: $app_output\n\n";
    if ($app_output eq "summary"){print $APP "Year,Q0,Q1,Q2,Q2,Q2,Q3,Q4,Mean,Count,Sum,Variance\n";}
    
    ### Print data
    foreach my $year (@years){
        if  (($data_retrieval_mode eq 'all') or ($data_retrieval_mode eq 'medline')){
            print $SS "$year\t$pubmed_count{$year}\t$medline_authored_publication_num{$year}\t$medline_unique_author_num{$year}\t$medline_unique_shortauthor_num{$year}\t$medline_duplicate_author_num{$year}\t$medline_duplicate_shortauthor_num{$year}\t$medline_duplicate_publication_num{$year}\t";
            print $SS "$retraction_num{$year}\t$made_ret_num{$year}\t$unique_journals{$year}\t$unique_journals_short{$year}\t$medline_authors_per_pub{$year}\t$medline_app_max{$year}\t$medline_app_min{$year}\t$medline_app_var{$year}\t$medline_unique_shortfirstau_num{$year}\n";
            print $APP "$year,";
            if (defined $medline_pub_list{$year}){
                if ($app_output eq "summary"){
                    my $stat = Statistics::Descriptive::Full->new();
                    $stat -> add_data(@{$medline_pub_list{$year}});
                    my $sum_squares = 0;
                    foreach my $i (@{$medline_pub_list{$year}}){$sum_squares = $sum_squares + ($i*$i);}
                    my $variance = ($sum_squares / ($stat->count()))-(($stat->mean())*($stat->mean()));
                    print $APP $stat->quantile(0),",",$stat->quantile(1),",",$stat->median(),",",$stat->median(),",",$stat->median(),",",$stat->quantile(3),",",$stat->quantile(4),",",$stat->mean(),",",$stat->count(),",",$stat->sum(),",$variance";
                }
                elsif ($app_output eq "full"){
                    foreach my $k (0..$#{$medline_pub_list{$year}}){print $APP "$medline_pub_list{$year}[$k],";}
                }
            }
            print $APP "\n";
        }
        elsif  ($data_retrieval_mode eq 'csv'){
            print $SS "$year\t$pubmed_count{$year}\t$csv_authored_publication_num{$year}\t$csv_unique_author_num{$year}\t$csv_unique_author_num{$year}\t$csv_duplicate_author_num{$year}\t$csv_duplicate_author_num{$year}\t$csv_duplicate_publication_num{$year}\t";
            print $SS "$retraction_num{$year}\t$made_ret_num{$year}\t$unique_journals{$year}\t$csv_authors_per_pub{$year}\t$csv_app_max{$year}\t$csv_app_min{$year}\t$csv_app_var{$year}\n";
            print $APP "$year,";
            if (defined $csv_pub_list{$year}){
                foreach my $k (0..$#{$csv_pub_list{$year}}){
                    print $APP "$csv_pub_list{$year}[$k],";
                }
            }
            print $APP "\n";
        }
        if ($data_retrieval_mode eq 'all'){
            print $MODES "$year\t$pubmed_count{$year}\t$csv_authored_publication_num{$year}\t$medline_authored_publication_num{$year}\t$csv_unique_author_num{$year}\t$medline_unique_author_num{$year}\t$medline_unique_shortauthor_num{$year}\t";
            print $MODES "$csv_duplicate_author_num{$year}\t$medline_duplicate_author_num{$year}\t$medline_duplicate_shortauthor_num{$year}\t$csv_duplicate_publication_num{$year}\t$medline_duplicate_publication_num{$year}\t";
            print $MODES "$csv_authors_per_pub{$year}\t$medline_authors_per_pub{$year}\t$csv_app_max{$year}\t$medline_app_max{$year}\t$csv_app_min{$year}\t$medline_app_min{$year}\t$csv_app_var{$year}\t$medline_app_var{$year}\n";
        }
    }
    close $SS;
    close $MODES;
    
    unless ($keyword_mode eq "none"){
        my $keyword_filename = "$file_directory"."keywords_"."$keyword_mode"."_$years[0]to$years[$#years].tab";
        open my $KEY, '>', $keyword_filename;
        
        print $KEY "Date: $date\t Scientist Search $version\n";
        print $KEY "Data Retrieval Mode: $data_retrieval_mode, reflects file type source of PubMed data in following analysis.\n";
        print $KEY "Publication Type Selection: $publication_type, Affiliation Restriction: $affiliation_restriction, reflects search constraints of PubMed data in following analysis.\n";
        print $KEY "Keyword Mode Selection: $keyword_mode\n\n";
        print $KEY "Keyword\t";
        foreach my $year (@years){print $KEY "$year\t";}
        print $KEY "\n";
        
        if ($keyword_mode eq "all"){
            foreach my $keyword (keys(%keyword_numbers)){
                print $KEY "$keyword\t";
                foreach my $y (0..$#years){
                    print $KEY "$keyword_numbers{$keyword}[$y]\t";
                }
                print $KEY "\n";
            }
        }
        elsif ($keyword_mode eq "query"){
                foreach my $k (@keyword_query){
                    print $KEY "$k\t";
                    foreach my $y (0..$#years){
                        my $result_sum = 0;
                        foreach my $keyword (keys(%keyword_numbers)){
                            if ($keyword =~ m/$k/){
                                $result_sum = $result_sum + $keyword_numbers{$keyword}[$y];
                            }
                        }
                        print $KEY "$result_sum\t";
                    }
                    print $KEY "\n";
                }
            }
        close $KEY;
    }
    
    ### Create Journal Data Output files and print headers and data
    if (scalar(@journal_years) > 0){
        my $journal_pub_filename = "$file_directory"."journal_pubs_$journal_years[0]to$journal_years[$#journal_years].tab";
        my $journal_impact_filename = "$file_directory"."journal_impact_$journal_years[0]to$journal_years[$#journal_years].tab";
        open my $JPUBS, '>', $journal_pub_filename;
        open my $JIF, '>', $journal_impact_filename;
        
        print $JPUBS "Date: $date\t Scientist Search $version\n";
        print $JPUBS "Data Retrieval Mode: $data_retrieval_mode, reflects file type source of PubMed data in following analysis.\n";
        print $JPUBS "Publication Type Selection: $publication_type, Affiliation Restriction: $affiliation_restriction, reflects search constraints of PubMed data in following analysis.\n";
        print $JPUBS "Keyword Mode Selection: $keyword_mode\n\n";
        print $JPUBS "Journal\t";
        foreach my $year (@journal_years){print $JPUBS "$year\t";}
        print $JPUBS "Total\n";
        
        print $JIF "Date: $date\t Scientist Search $version\n";
        print $JIF "Data Retrieval Mode: $data_retrieval_mode, reflects file type source of PubMed data in following analysis.\n";
        print $JIF "Publication Type Selection: $publication_type, Affiliation Restriction: $affiliation_restriction, reflects search constraints of PubMed data in following analysis.\n";
        print $JIF "Keyword Mode Selection: $keyword_mode\n\n";
        print $JIF "Journal\t$journal_impact_header\n";
        
        foreach my $journal (keys(%journal_publication_num)){
            print $JPUBS "$journal\t";
            print $JIF "$journal\t";
            foreach my $y (0..$#journal_years){
                print $JPUBS "$journal_publication_num{$journal}[$y]\t";
            }
            foreach my $j_match (keys %journal_if_matchlist){
                if ($journal eq $j_match){
                    my $data = join "\t", @{$journal_if_matchlist{$journal}};
                    print $JIF "$data";
                }
            }
            print $JPUBS "$journal_publication_total{$journal}\n";
            print $JIF "\n";
        }
        close $JPUBS;
        close $JIF;
    }
}
# ================================================================== FUNCTIONS ======================================================================== #
#~~~~~~~~~~~~~~~~~~~ 1 ~~~~~~~~~~~~~~~~~~~~~#
# Open a file and return an array containing content
# Input: $filename
# Output: @content
sub get_file_content {
   use strict;
   use warnings;
   my ($name) = @_;
   open (my $FILE,'<',$name) or die "Cannot open file = $name\n";
   my @data;
   while ( my $line = <$FILE>) {
	chomp $line;
	push @data, $line;
   }  
   close $FILE;
   return @data;
}
#~~~~~~~~~~~~~~~~~~~ 2 ~~~~~~~~~~~~~~~~~~~~~#
# Given an array, return only unique values of the array in @unique, and duplicates in @duplicates, preserving original order. The %seen hash lists the frequencey of each array value
# Input: \@array
# Output: \@unique, \@duplicates, \%seen
sub unique_array {
    my ($array_ref) = @_;
    my @array = @$array_ref;
    my %seen = ();
    my @unique;
    my @duplicates;
    foreach my $k (@array) {
        if (!defined $seen{$k}) {
            push @unique, $k;
            $seen{$k} = 1;
        }
        elsif ($seen{$k} == 1) {
            push @duplicates, $k;
            $seen{$k} = $seen{$k} + 1;
        }
        elsif ($seen{$k} > 1) {
            $seen{$k} = $seen{$k} + 1;
        }
    }
    return (\@unique, \@duplicates, \%seen);
}
#~~~~~~~~~~~~~~~~~~~ 3 ~~~~~~~~~~~~~~~~~~~~~#
# Given an array of numbers, return the mean, max, min, and variance
# Input: @array
# Output: $mean, $max, $min, $variance
sub array_stats {
    use strict;
    use warnings;
    use List::Util qw (max min);
   
    my ($name) = @_;
    my @data = @$name;
    my $max = max (@data);
    my $min = min (@data);
    my $data_sum = 0;
    foreach my $i (@data){
        $data_sum = $data_sum + $i;
    }
    my $mean;
    if (scalar (@data) > 0){
        $mean = $data_sum / scalar(@data);
    }
    else {$mean = 0;}
    my $var_sum = 0;
    foreach my $i (@data){
        $var_sum = $var_sum + (($i-$mean)*($i-$mean));
    }
    my $var;
    if (scalar (@data) > 0){
        $var = $var_sum / scalar(@data);
    }
    else {$var = 0;}
    
    return ($mean, $max, $min, $var);
}
