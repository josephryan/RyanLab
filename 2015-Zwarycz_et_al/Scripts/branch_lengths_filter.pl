#!/usr/bin/perl

use strict;
use warnings;

MAIN: {
   my $tree_branch_lengths = $ARGV[0] or die "usage: $0 TREEFILE\n";
   my $branch_length_results = filter_branch_lengths($tree_branch_lengths);
}

sub filter_branch_lengths {
   my $file = shift;
   open IN, $file or die "cannot open $file:$!";
   my $max_length = 2.3;
   while (my $line = <IN>) {
      my @fields = split /,/, $line;
      foreach my $F (@fields) {
         $F =~m/([^\(\):]+):([0-9E\.\-]+)/;
         if ($2 > $max_length) {
            print "$1:$2\n";
         }
      }
   }
}   
