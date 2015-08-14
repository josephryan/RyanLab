#!/usr/bin/perl

use strict;
use warnings; 
use Data::Dumper;
use JFR::Fasta;

MAIN: {
    my $blast = $ARGV[0] or die "usage: $0 BLASTFILE FASTA CUTOFF OUTDIR DATABASE\n";
    my $fasta = $ARGV[1] or die "usage: $0 BLASTFILE FASTA CUTOFF OUTDIR DATABASE\n";
    my $cutoff = $ARGV[2] or die "usage: $0 BLASTFILE FASTA CUTOFF OUTDIR DATABASE\n";
    my $outdir = $ARGV[3] or die "usage: $0 BLASTFILE FASTA CUTOFF OUTDIR DATABASE\n";
    my $db = $ARGV[4] or die "usage: $0 BLASTFILE FASTA CUTOFF OUTDIR DATABASE\n";
    my $rh_evalues = get_evalues($blast);
    my $rh_all_seqs = get_all_seqs($fasta);
    if ($outdir) {
        my $seq_str = get_seqs($rh_evalues, $cutoff, $rh_all_seqs);
        my $blast_results = run_blast($seq_str, $outdir, $db);
    } else {
            print_evalues($rh_evalues);
    }
}

sub run_blast {
    my $seq_str = shift;
    my $outdir = shift; 
    my $db = shift;
    mkdir $outdir unless(-d $outdir);
    my $file = "$outdir/seqs.fa";
    open OUT, ">$file" or die "cannot open $file:$!";
    print OUT $seq_str;
    my $cmd = "blastp -query $file -evalue 10 -num_descriptions 3 -num_alignments 3 -db $db -out $outdir/blastp.txt";
    system $cmd; 
    print "blast report is available in $outdir/blastp.txt\n";
}

sub get_all_seqs {
    my $fasta = shift;
    my %seqs = ();
    my $fp = JFR::Fasta->new($fasta);
    while (my $rec = $fp->get_record()) {
        my $id = JFR::Fasta->get_def_w_o_gt($rec->{'def'});
        $seqs{$id} = $rec->{'seq'};
    }
    return \%seqs;
}

sub get_seqs {
    my $rh_evalues = shift;
    my $cutoff = shift;
    my $rh_all_seqs = shift;
    my $seq_str = '';
    my @sorted = sort {$rh_evalues->{$a} <=> $rh_evalues->{$b}} keys %{$rh_evalues}; 
    foreach my $id (@sorted) {
        last if ($rh_evalues->{$id} > $cutoff);
        #print "$id = $rh_evalues->{$id}\n";
        $seq_str .= ">$id $rh_evalues->{$id}\n $rh_all_seqs->{$id}\n";
    }
    return $seq_str;
}

sub print_evalues {
    my $rh_evalues = shift;
    my @sorted = sort {$rh_evalues->{$a} <=> $rh_evalues->{$b}} keys %{$rh_evalues};
    foreach my $id (@sorted) {
        print "$id = $rh_evalues->{$id}\n";
    }
}

sub get_evalues {
    my $file = shift;
    open IN, $file or die "cannot open $file:$!";
    my %evalues = ();
	while (my $line = <IN>) {
        my @fields = split /\s+/, $line;
            if ($evalues{$fields[0]}) {
		if ($fields[10] < $evalues{$fields[0]}) {
			$evalues{$fields[0]} = $fields[10];
		}
	    }else {
		$evalues{$fields[0]} = $fields[10];
	    }
	}
	return \%evalues;
}
