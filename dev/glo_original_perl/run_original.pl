#!/usr/bin/perl
# Runs the ORIGINAL, unchanged OptimizerTools::optimize of D. J. Dickinson's GLO code (https://github.com/dannyhmg/germline, MIT licence;
# frozen copy https://github.com/elifesciences-publications/germline) on a protein, with stand-ins (in ./stubs) for the three Perl modules that need
# a database, a random-number library or a progress indicator. The stand-in for BerkeleyDB reads the scores from a flat file of 4^12 doubles.
#
#   SCOREFILE=scores.f64 perl -I dev/glo_original_perl/stubs dev/glo_original_perl/run_original.pl <protein> <seed> <path to the repository's cgi-bin/libs>
# Prints the DNA sequence it designs and its score (the mean score per codon).
# See dev/verify_glo_against_original.R and docs/GLO_VERIFICATION.md.
use strict; use warnings;
my ($aa, $seed, $libs) = @ARGV;
unshift @INC, $libs;
require BerkeleyDB; require OptimizerTools; require Seqscore;
srand($seed);
my $lib = BerkeleyDB::Btree->new();
my $r = OptimizerTools::optimize($lib, $aa);
printf "%s\t%.6f\n", $r->{Sequence}, $r->{Sequence_score};
