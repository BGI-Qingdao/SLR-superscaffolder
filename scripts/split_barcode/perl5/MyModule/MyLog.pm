package MyModule::MyLog;

use warnings;
use strict;
use POSIX;
use Exporter;
use vars qw(@ISA @EXPORT);


@ISA = qw(Exporter);
@EXPORT = qw(&printInfo2ERR);


sub printInfo2ERR {
    my($msg, ) = @_;
    
    my $logTime =  strftime("[%Y-%m-%d %H:%M:%S]", localtime());
    print STDERR "$logTime INFO: $msg\n";
}


1;