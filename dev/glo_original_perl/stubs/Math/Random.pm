package Math::Random;
use strict; use Exporter 'import';
our @EXPORT = qw(random_uniform_integer random_uniform);
sub random_uniform_integer { my ($n, $lo, $hi) = @_; return $lo + int(rand($hi - $lo + 1)); }
sub random_uniform { return rand(); }
1;
