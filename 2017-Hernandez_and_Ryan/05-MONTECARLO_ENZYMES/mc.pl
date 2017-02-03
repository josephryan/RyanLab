#!/usr/bin/perl

use strict;
use warnings;
use List::Util 'shuffle';
use FileHandle;
use IO::Uncompress::Gunzip qw($GunzipError);
use Data::Dumper;

our $VERSION = 0.02;

our $SRAND = 42;
our $REPS = 1000;

our $DIR = '00-DATA';
our @FILES = ("$DIR/ecoli.pep_v_uniprot_inclusive.blastp.gz",
              "$DIR/bsub.pep_v_uniprot_inclusive.blastp.gz",
              "$DIR/ccre.pep_v_uniprot_inclusive.blastp.gz",
              "$DIR/mgen.pep_v_uniprot_inclusive.blastp.gz",
              "$DIR/afis.pep_v_uniprot_inclusive.blastp.gz",
              "$DIR/spcc.pep_v_uniprot_inclsuive.blastp.gz",
              "$DIR/pflu.pep_v_uniprot_inclusive.blastp.gz"
             );
# %TOTAL = (species => [NUM_ENZ_IN_BACTERIAL_HGT => TOTAL_BACTERIAL_HGT]);
our %TOTAL = ('Mnemiopsis_leidyi' => [7,8],
              'Amphimedon_queenslandica' => [18,20]);
# Amphimedon is a placeholder

MAIN: {
    srand($SRAND);
    my %gte = ();  # (species => count gte $TOTAL{species}->[0]);
    my @alldata = ();
    foreach my $file (@FILES) {
        my %seen = ();
        my $fh = IO::Uncompress::Gunzip->new($file)
            or die "IO::Uncompress::Gunzip of $file failed: $GunzipError\n";
        my @data = ();
        while(my $line = $fh->getline() ) {
            my @f = split /\t/, $line;
            push @data, $f[1] unless ($seen{$f[0]}); 
            $seen{$f[0]}++;
        }
        push @alldata, \@data;
    }
    foreach my $sp (keys %TOTAL) {
        for (my $h = 0; $h < @alldata; $h++) {    
            for (my $i = 0; $i < $REPS; $i++) {
                my @shuffled = shuffle(@{$alldata[$h]});
                my $count = 0;
                for (my $j = 0; $j < $TOTAL{$sp}->[1]; $j++) {
                    $count++ if ($shuffled[$j] =~ m/^EC/ ||
                                 $shuffled[$j] =~ m/^ASE/);
                }
                $gte{$sp}++ if $count >= $TOTAL{$sp}->[0];
            }
        }
    }
    foreach my $sp (sort keys %gte) {
        my $denominator = scalar(@FILES) * $REPS;
        my $pval = $gte{$sp} / $denominator;
        print "$sp: P=$pval ($gte{$sp}/$denominator)\n";
    }
}

