
#USAGE: perl BEACON.pl proteins_of_interest.txt cache_directory 

use warnings;
use strict;

my %cache_files;
opendir my $CACHE_DIR, $ARGV[1] || die "Cannot open directory: $!";
my @cache_list = readdir $CACHE_DIR;
foreach my $c (0..$#cache_list){
    $cache_files{$cache_list[$c]} = 1;
}
closedir $CACHE_DIR;

my @poi_list = ();

open(my $POI, $ARGV[0]) || die "Cannot open protein list of interest: $!";
while(<$POI>){
    $_ =~ /(.+)/g;
    my $poi = $1;
    push(@poi_list, $poi);
}
close $POI;

my $out_dir = $ARGV[2];

my $outfile = $out_dir."/protein_of_interest_interactions.txt";
open(my $poi_out, '>', $outfile) or die "Could not open file '$outfile': $!";


my @all_scores = ();
foreach my $p (0.. $#poi_list-1){
    if ($cache_files{$poi_list[$p]}){
	my $protein_file = $ARGV[1]."/".$poi_list[$p];
	open(my $PF, $protein_file);
	
	my %poi2 = ();
	while(<$PF>){
	    $_ =~ /(.+)\s(\d+)/g;
	    my $poi2 = $1;
	    my $score = $2;
	    $poi2{$poi2} = $score;
	}
	close $PF;

	foreach my $q ($p+1..$#poi_list){
	    if (exists($poi2{$poi_list[$q]})){
		print $poi_out "$poi_list[$p]\t$poi_list[$q]\t$poi2{$poi_list[$q]}\n";
		push(@all_scores, $poi2{$poi_list[$q]});
	    }
	    else{
		print $poi_out "$poi_list[$p]\t$poi_list[$q]\t20\n";
		push(@all_scores, 20);
		
	    }
	}
    }
    
    else {
	print "Protein of interest $poi_list[$p] is not in the protein cache, please verify it is in proper ENSEMBL ID format. Ex: 9606.ENSP00000501317\n";
    }    
}

close $poi_out;

#print "@all_scores\n";

my $average = avg(@all_scores);
#print "$average\n";
print "$average";

sub avg {
    my $total;
    $total += $_ foreach @_;
    return $total / @_;
}

