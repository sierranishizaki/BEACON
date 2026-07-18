
#USAGE: perl run_BEACON_candidates.pl proteins_of_interest.txt known_proteins.txt cache_directory 
#!/usr/bin/env perl                                                            
                                                                              
use warnings;                                                                  
use strict;
use Getopt::Long qw(GetOptions);


my $start_time = time();

my $help_opt;
my $poi_opt = '';
my $known_opt = '';
my $cache_opt = 'cache';
my $output_opt = 'output';
GetOptions(
    'help|h'    => \$help_opt,
    'input|i=s' => \$poi_opt,
    'known|k=s' => \$known_opt,
    'cache|c=s' => \$cache_opt,
    'output|o=s' => \$output_opt,
) || die "Usage: $0 --i proteins_of_interest.txt  --cache cache_directory --o output_directory\n";

if ($help_opt) {
    print "
▒█▀▀█ ▒█▀▀▀ ░█▀▀█ ▒█▀▀█ ▒█▀▀▀█ ▒█▄░▒█
▒█▀▀▄ ▒█▀▀▀ ▒█▄▄█ ▒█░░░ ▒█░░▒█ ▒█▒█▒█
▒█▄▄█ ▒█▄▄▄ ▒█░▒█ ▒█▄▄█ ▒█▄▄▄█ ▒█░░▀█
version 0.1
Written by Dr. Sierra Sky Nishizaki and Pelle Hall
9/11/25
BEACON is a novel method to query the relatedness of a candidate gene network by measuring the degrees of separation of protein-protein functional associations from the STRING database

Run guide:
1. First, please download all protein links from STRING. Ex:
wget \"https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz\"                                                           
2. Next, run beaconPathfind.py to generate a cache of the shortest number of functional links between each pair of proteins in the database. Ex: 
python beaconPathfind.py --input 9606.protein.links.v12.0.txt.gz --threshold 400 --output cache                                                                
3. Generate a BEACON score for each protein in a candidate list alongside known disease-associated proteins. Ex:

perl run_BEACON_candidates.pl --input proteins_of_interest.txt --known known_disease_proteins.txt --cache cache_directory --output output_directory";
}

if ($poi_opt eq '' || $known_opt eq '') {
    print "Usage: $0 --i proteins_of_interest.txt  --cache cache_directory --o output_directory\n";
}

if(-e $output_opt){
    print "Output directory $output_opt already exists, please remove $output_opt to continue\n";
    exit;
}
else{
    my $run_out = "mkdir ./".$output_opt;
    my $mkoutdir = `$run_out`;
}

open(my $POI, '<',$poi_opt);
my @poi = '';
while(<$POI>){
    $_ =~ /(.+)/g;
    my $poi = $1;
    push @poi, $poi;
}

open(my $KNOWN, '<',$known_opt);
my @known = '';
while(<$KNOWN>){
    $_ =~ /(.+)/g;
    my $known = $1;
    push @known, $known;
}


my $temp_out = $output_opt."/unsorted_beacon_candidate_numbers.txt";
open(my $temp_file, '>>',$temp_out);

foreach my $p (1..$#poi){
    my @full_list = @known[1..$#known];
    push @full_list, $poi[$p];

    my $cp_out = "./candidate_proteins.txt"; 
    open(my $cp_file, '>',$cp_out);
    print $cp_file join("\n",@full_list);

    `mkdir $output_opt/$poi[$p]`;
    
    my $beacon_out = `perl BEACON.pl candidate_proteins.txt $cache_opt $output_opt/$poi[$p]`;
    print "$beacon_out\n";
    print $temp_file "$poi[$p] Wu2020 BEACON score:\t$beacon_out\n";
}    

my $out_out = $output_opt."/beacon_candidate_numbers.txt";
#open(my $out_file, '>' $out_out);
`sort -k 2 $temp_out > $out_out`;


`rm ./candidate_proteins.txt`;
`rm $temp_out`;
