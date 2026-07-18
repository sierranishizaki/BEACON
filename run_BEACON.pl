
#USAGE: perl run_BEACON.pl proteins_of_interest.txt cache_directory 
#!/usr/bin/env perl                                                            
                                                                              
use warnings;                                                                  
use strict;
use Getopt::Long qw(GetOptions);


my $start_time = time();

my $help_opt;
my $poi_opt = '';
my $cache_opt = 'cache';
my $output_opt = 'output';
GetOptions(
    'help|h'    => \$help_opt,
    'input|i=s' => \$poi_opt,
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
3. Generate a BEACON score and compare it to random protein lists using run_BEACON.pl. Ex:                                                                     
perl run_BEACON.pl --input proteins_of_interest.txt --cache cache_directory --output output_directory";
}

if ($poi_opt eq '') {
    print "Usage: $0 --i proteins_of_interest.txt --cache cache_directory --o output_directory\n";
    exit;
}

if(-e $output_opt){
    print "Output directory $output_opt already exists, please remove $output_opt to continue\n";
    exit
}
else{
    print "Generating $output_opt directory\n";
    my $run_out = "mkdir ./".$output_opt;
    my $mkoutdir = `$run_out`;
}

my $bs_file = $output_opt."/BEACON_scores.txt";
open(my $bs, '>',$bs_file);

#print "Running: perl BEACON.pl $poi_opt $cache_opt $output_opt\n";
#my $beacon_out = `perl BEACON.pl $poi_opt $cache_opt $output_opt`;
my $l = `wc -l $poi_opt | cut -f1 -d\\s`;
$l =~ /(\d+)/g;
my $length = $1;
print "Running: perl simulations.pl $cache_opt $length $output_opt\n";
my $sim_out = `perl simulations.pl $cache_opt $length $output_opt`;

print "Running: perl BEACON.pl $poi_opt $cache_opt $output_opt\n";
my $beacon_out = `perl BEACON.pl $poi_opt $cache_opt $output_opt`;

my $rand_out_file = $output_opt."/random_protein_scores.txt";
open(my $sims, '<',$rand_out_file);

my $sim_count = 0;
my @sim_array = ();
while(<$sims>){
    $_ =~ /(.+)/g;
    my $sim_score = $1;
    push @sim_array, $sim_score;
    if($sim_score <= $beacon_out){
	$sim_count++;
    }
}
my $sim_avg = avg(@sim_array);

print $bs "Proteins of interest BEACON score: $beacon_out\n";
print $bs "Random proteins BEACON score: $sim_avg\n";

my $p_val = ($sim_count+1)/20001;

print $bs "Monte Carlo Approximation P-value: $p_val\n";

sub avg {
    my $total;
    $total += $_ foreach @_;
    return $total / @_;
}

close $bs;

my $end_time = time();
my $run_duration = ($end_time - $start_time)/60;
print "Run Duration: $run_duration minutes\n";
