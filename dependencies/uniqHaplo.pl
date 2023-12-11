#!/usr/bin/perl

# Copyright 2015, Naoki Takebayashi <ntakebayashi@alaska.edu>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

# Version: 20161129

my $usage =
    "Usage: $0 [-ha] [inputfile]\n".
    "  -h: help\n" .
    "  -a: inputfile is amino acid seq, instead of DNA\n".
    "\n" .
    "This program read in the sequence file, which may contain sequences " .
    "with opposite orientation.  It identifies identical alleles by going " .
    "through all pairwise comparisons.  When the shorter sequence of the two ".
    "is identical to the substring of the longer, they are considered as an " .
    "allele.  Gaps '-' will be removed before the comparison. The longest " .
    "sequences of each allele will be printed in STDOUT.  These unique ".
    "alleles are formated in FASTA.  Information about " .
    "which alleles are identical and the difference in the lengths are ".
    "printed in STDERR. When the sequences with the opposite direction are" .
    "included, it makes the complement of the sequences" .
    "and the comparison is made.\n";

# Changelog
# * 20190207
# - bug fix: RevCompIUPAC was not making a copy of Seq object, so some sequences
#   in the output were ending up in reverse complelement of the input.
# - STDIN can be the input file now.
# * 20161129
# - RevCompIUPAC was supposed to return Bio::Seq object, but it was returning DNA
#   sequence string. This bug was fixed. Thank you for Veronica El Mujtar for
#   noticing this problem.
# * 20151028
# - QsubSeq wasn't taking into account of regex special character.
#   More specifically if the sequence have * or ?, some sequences which
#   differ around these special characters were considered identical.
#   This was fixed by using \Q ... \E in the pattern match.
# - removed EMBOSS reliance to make reverse compelement
# - Changed the STDERR output about the information of identity, now the selected
#   one is printed in the beginning.

use Bio::AlignIO;
use Bio::SeqIO;
# use Bio::Factory::EMBOSS;
use IO::File;
#use POSIX qw(tmpnam);
use Getopt::Std;

getopts('ah') || die "$usage\n";
die "$usage\n" if (defined($opt_h));

@ARGV = ('-') unless @ARGV;
my $infile = shift;

# read in the data
$in =($infile eq '-') ?
    Bio::SeqIO->new(-fh => \*STDIN, '-format' => 'fasta') :
    Bio::SeqIO->new(-file => $infile, '-format' => 'fasta');

my $outFH = Bio::SeqIO->newFh(-format => 'fasta');

my @seqArray =();
while(my $seq = $in->next_seq()) {
    my $tmpSeq = $seq->seq();
#    $tmpSeq =~ s/-//g;  # remove gaps
    $seq->seq($tmpSeq);

    push @seqArray, $seq;
}

# compare seq strings
my @uniq = ();
my ($i, $j);
for ($i=0; $i < @seqArray; $i++) {
    my $keep = $seqArray[$i] ;
    my @identical = ($seqArray[$i]->display_id() . "(ref)");
    for ($j=$i+1; $j < @seqArray; $j++) {
	my $comp = QsubSeq($seqArray[$i], $seqArray[$j]);
	if ($comp eq "different" && (!defined($opt_a))) {
	    # check the complement
	    $comp = QsubSeq($seqArray[$i], RevCompIUPAC($seqArray[$j]));
	}
	next if ($comp eq "different");
	
	# They matches.
	
	if ($comp < 0) {
	    $keep = $seqArray[$j];  # keeping the longer seq
	    unshift @identical, $seqArray[$j]->display_id() . " (" . -$comp . ")";	    
	} else {
	    push @identical, $seqArray[$j]->display_id() . " (" . -$comp . ")";
	}
	
	splice(@seqArray, $j, 1); # remove this
	$j = $j - 1;
	
    }
    
    if (@identical > 1) {
	print STDERR join " = ", @identical;
	print STDERR "\n";
    }
    
    push @uniq, $keep;
}


# output
foreach my $i (@uniq) {
    print $outFH $i;
}

exit(0);


# takes two sequences, and returns an signed integer when one seq is
# a substring of the other.  The returned value is the difference 
# in the lengths of the seq
sub QsubSeq {
    my ($seq1, $seq2) = @_;

    my $s1 = lc($seq1->seq());  # setting them to lower case
    my $s2 = lc($seq2->seq());

    # remove gaps '-'
    $s1 =~ s/-//g;
    $s2 =~ s/-//g;

    my $len1 = CharLen($s1);
    my $len2 = CharLen($s2);

    if ($len1 > $len2) {
	unless ($s1 =~ /\Q$s2\E/) {
	    return ("different");
	}
    } elsif ($len2 > $len1) {
	unless ($s2 =~ /\Q$s1\E/) {
	    return ("different");
	}
    } else {
	unless ($s1 eq $s2) {
	    return ("different");
	}
    }
    # when reached here, 1 seq is substring of the other, so return
    # the length diff of the two seq (seq1 - seq2)
    return ($len1 - $len2)
}

# returns reverse complement
sub RevCompIUPAC {
        my $seqObj = shift;

	my $seq = Bio::Seq->new(-seq        => $seqObj->seq,
				-display_id => $seqObj->display_id,
				-desc       => $seqObj->desc,
				-alphabet   => $seqObj->alphabet);
	my $dnaSeq = $seq->seq();
	# reverse the DNA sequence
        my $revcomp = reverse($dnaSeq);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ABCDGHKMNRSTUVWYabcdghkmnrstuvwy/TVGHCDMKNYSAABWRtvghcdmknysaabwr/;

	$seq->seq($revcomp);
        return $seq;
}

# I used to use this to make reverse complement.
# sub Complement {
#     my $seq = shift;

#     my $tempOutfile;
#     my $fh;
#     do {$tempOutfile = tmpnam()} 
#     until $fh = IO::File->new($tempOutfile, O_RDWR|O_CREAT|O_EXCL);
#     $fh->close;

#     my $factory = new Bio::Factory::EMBOSS;
#     my $revseq = $factory->program('revseq');

#     my %input = (-sequence => $seq,
# 		 -outseq => $tempOutfile );

#     $revseq->run(\%input);
#     my $seqio = Bio::SeqIO->new (-file => $tempOutfile);

#     unlink $tempOutfile || die "ERROR: Unable to unlink $tempOutfile\n";

#     return($seqio->next_seq())
# }

sub CharLen {
    my $string = shift;
    my @charString = split (//, $string);
    return scalar(@charString);
}

