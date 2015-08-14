#!/usr/bin/perl

use strict;
use warnings; 
use JFR::Fasta;

MAIN: {
    my $hmm   = $ARGV[0] or die "usage: $0 HMMSEARCHFILE FASTA\n";
    my $fasta = $ARGV[1] or die "usage: $0 HMMSEARCHFILE FASTA\n";
    my $rh_ids = get_ids($hmm);
    print_if_id_not_in_hash($fasta,$rh_ids);
}

sub print_if_id_not_in_hash {
    my $file = shift;
    my $rh_ids = shift;
    my $fp = JFR::Fasta->new($file);
    while (my $rec = $fp->get_record()) {
        my $id = JFR::Fasta->get_def_w_o_gt($rec->{'def'});
        next if ($rh_ids->{$id});
        print "$rec->{'def'}\n";
        print "$rec->{'seq'}\n";
    }
}

sub get_ids {
    my $file = shift;
    my %ids = ();
    open IN, $file or die "cannot open $file:$!";
    while (my $line = <IN>) {
        next if ($line =~ m/^#/);
        next if ($line =~ m/^\s*$/);
        next if ($line =~ m/^Query:/);
        next if ($line =~ m/^Scores/);
        next if ($line =~ m/^   ---/);
        next if ($line =~ m/^    E-value/);
        next if ($line =~ m/^    -------/);
        last if ($line =~ m/inclusion threshold/);
        last if ($line =~ m/^Domain annotation/);
        chomp $line;
        my @fields = split /\s+/, $line;
        $ids{$fields[9]} = 1;
    }
    return \%ids;
}

