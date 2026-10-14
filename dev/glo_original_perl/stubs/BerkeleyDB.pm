package BerkeleyDB;
package BerkeleyDB::Btree;
# Stand-in for the BerkeleyDB module: db_get reads the score of a 12-mer from a flat file of 4^12 doubles (alphabetical order of the 12-mers).
use strict;
my %v = (A=>0, C=>1, G=>2, T=>3);
sub new { my ($class, %a) = @_; open(my $fh, '<:raw', $ENV{SCOREFILE}) or die "no score file"; return bless { fh => $fh }, $class; }
sub db_get {
    my $self = shift;
    my $key = $_[0];
    return 1 unless $key =~ /^[ACGT]{12}$/;
    my $idx = 0;
    $idx = $idx * 4 + $v{$_} for split //, $key;
    sysseek($self->{fh}, $idx * 8, 0);
    my $buf; sysread($self->{fh}, $buf, 8);
    $_[1] = unpack("d<", $buf);
    return 0;
}
1;
