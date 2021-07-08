#!usr/bin/perl
# tophit_ectc.pl
use warnings;
use strict;
use Data::Dumper;

our $BLAST = $ARGV[0] or die "usage: $0 BLAST\n";

MAIN: {
    my $rh_th = get_top_hit ($BLAST);
    my $rh_ectc = get_ectc ($rh_th);
    my $count= get_counts ($rh_ectc);
}

sub get_counts {
    my $rh_ectc = shift;
    my %count = ();
    foreach my $val (values %{$rh_ectc}){
        if ($val eq "ectc"){
            $count{$val}++;
        }elsif ($val eq "ec"){
            $count{$val}++;
        }elsif ($val eq "ase"){
            $count{$val}++;
        }elsif ($val eq 'tc'){
            $count{$val}++;
        }elsif ($val eq 'tra'){
            $count{$val}++;
        }elsif ($val eq 'unchar'){
            $count{$val}++;
        }elsif ($val eq 'other'){
            $count{$val}++;
        }
    }
        foreach my $val (sort keys %count){
        print "$val, $count{$val}\n";
    }
}

sub get_top_hit {
    my $blast_file = shift;
    my %tophits = ();
    open IN, $blast_file or die "cannot open $blast_file:$!";
    while (my $line = <IN>){
        my @fields = split /\s+/, $line;
        next if ($tophits{$fields[0]});
        $tophits{$fields[0]} = $fields[1];

    }
    return \%tophits;
}
sub get_ectc {
    my $rh_th = shift;
    my %ectc = ();
    foreach my $key (keys %{$rh_th}){
        if ($rh_th->{$key}=~m/^ECTC\|/){
            $ectc{$key} = "ectc";
        }elsif ($rh_th->{$key}=~m/^EC\|/){
            $ectc{$key} = "ec";
        }elsif ($rh_th->{$key}=~m/^ASE\|/){
            $ectc{$key} = "ase";
        }elsif ($rh_th->{$key}=~m/^TC\|/){
            $ectc{$key} = "tc";
        }elsif ($rh_th->{$key}=~m/^TRA\|/){
            $ectc{$key} = "tra";
        }elsif ($rh_th->{$key}=~m/^UNCHAR\|/){
            $ectc{$key} = "unchar";
        }else{
            $ectc{$key} = "other";
        }
    }
    return \%ectc;
}

