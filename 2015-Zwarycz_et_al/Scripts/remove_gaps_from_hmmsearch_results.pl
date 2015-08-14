#!/usr/bin/perl

###########################################################################
# IMPORTANT TIP
###########################################################################
# ALWAYS MAKE SURE YOU RUN stockholm2fasta.pl WITH THE -g OPTION
###########################################################################


our $VERSION = 0.01;

use lib qw(/home/s9/jfryan/lib);
use strict;
use warnings;
use JFR::Fasta;
use Data::Dumper;


MAIN: {
    my $fasta_file = $ARGV[0] or die "usage: $0 FASTA_ALN\n";
    my %seqs = ();
    my $ra_gaps = get_gap_positions($fasta_file,\%seqs);
    print_sequences($ra_gaps,\%seqs);
}

sub print_sequences {
    my $ra_gaps = shift;
    my $rh_seqs = shift;
    foreach my $def (keys %{$rh_seqs}) {
        print "$def\n";
        for (my $i=0; $i < @{$rh_seqs->{$def}}; $i++) {
            print "$rh_seqs->{$def}->[$i]" unless ($ra_gaps->[$i]);
        }
        print "\n";
    }
}

sub get_gap_positions {
    my $file = shift;
    my $rh_seqs = shift;
    my @gaps = ();
    my $fp = JFR::Fasta->new($file);
    while (my $rec = $fp->get_record()) {
        my @nts = split /|/, $rec->{'seq'};
        $rh_seqs->{$rec->{'def'}} = \@nts;
        for (my $i = 0; $i < @nts; $i++) {
            $gaps[$i]++ if ($nts[$i] =~ m/[a-z*]/);
        }
#        if (scalar(@gaps)) {
#print "$gaps[10]\n";
#print "$nts[10]\n";
#            print Dumper \@gaps;
#            print Dumper \@nts;
#            exit;
#        }
    }
    return \@gaps;
}

