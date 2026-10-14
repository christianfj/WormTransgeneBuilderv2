package Progress::Any;
sub get_indicator { return bless {}, 'Progress::Any' }
sub update { return 1 }
1;
