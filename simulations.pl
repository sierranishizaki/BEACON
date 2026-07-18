
#USAGE: perl simulations.pl cache_directory length output

#!/usr/bin/env perl

use warnings;
use strict;

my $cache_dir = $ARGV[0];
opendir my $CACHE_DIR, $cache_dir || die "Cannot open directory: $!";
my @cache_list = readdir $CACHE_DIR;
closedir $CACHE_DIR;

my $outdir = $ARGV[2];

#my @count_array = (10,20,30,40,50,60,70,80,90,100,200,300,400,500); 
my @count_array = $ARGV[1];
foreach my $ca (0..$#count_array){
    my $start_time = time();
    print "Generating $count_array[$ca] random proteins";
    my @all_beacon_scores = ();
    foreach my $n (0..20000){
	my @random_genes;
	my $count = $count_array[$ca]; 
	my $max_value = $#cache_list;
	
	for (1..$count) {
	    push @random_genes, $cache_list[int(rand($max_value))];
	}
	
	my $tp = $outdir."/test_proteins.txt";
	open(my $out, '>', $tp) or die "Could not open file '$tp': $!";
	print $out join("\n", @random_genes);
	close $tp;
	
	my $beacon_out = `perl BEACON.pl $outdir/test_proteins.txt $cache_dir $outdir`;
	push @all_beacon_scores, $beacon_out;
    }
#    my $outfile = "simulations/".$count_array[$ca];
    my $outfile = $outdir."/random_protein_scores.txt";
    open(my $out, '>', $outfile) or die "Could not open file '$outfile': $!"; 
    print $out join("\n",@all_beacon_scores);
    close $out;

    my $end_time = time();
    my $run_duration = $end_time - $start_time;
    print "\t$run_duration\n";
}    




