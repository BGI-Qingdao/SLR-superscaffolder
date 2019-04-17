package MyModule::GlobalVar;

use warnings;
use strict;
use File::Spec;
use File::Basename;
use Cwd qw(abs_path);
use Config::IniFiles;
use Exporter;
use vars qw(@ISA @EXPORT);


@ISA = qw(Exporter);
@EXPORT = qw($PIPELINE_PATH $BIN_PATH $DB_PATH $TOOL_PATH $CONFIG_FILE $DB_LIST_BED
             $REF_H_CONFIG $PYTHON3 $MAX_AVAILEBLE_MEM $REF_H_DBPATH);


# Path
my $modulePath = dirname(File::Spec->rel2abs(__FILE__));
our $PIPELINE_PATH  = abs_path("$modulePath/../../..");
our $BIN_PATH       = "$PIPELINE_PATH/bin";
our $DB_PATH        = "$PIPELINE_PATH/db";
our $TOOL_PATH      = "$PIPELINE_PATH/tools";
our $CONFIG_FILE    = "$PIPELINE_PATH/conf/config.ini";
our $DB_LIST_BED    = "$DB_PATH/db.list";

# Config file
&getCfg($CONFIG_FILE);


# Kraken db info
&getBEDbInfo($DB_LIST_BED);


#===============================================================================
#   Subroutine
#===============================================================================
sub getCfg {
    my ($cfg_file, ) = @_;
    tie my %ini, "Config::IniFiles", (-file, => $cfg_file);

    our $REF_H_CONFIG      = \%ini;
    our $MAX_AVAILEBLE_MEM = $ini{"resource"}{"MaxAvailableMem"};
}

sub getBEDbInfo{
    my ($dbList, ) = @_;
    $dbList = abs_path($dbList);
    my $dbDir = dirname($dbList);
    
    my %h_dbPath = ();
    my %h_dbMem = ();
    open IN, "$dbList" or die $!;
    while (<IN>) {
        chomp;
        next unless /\S+/;
        
        my($dbName, $dbPath) = split /\s+/;
        $dbPath = "$dbDir/$dbPath"; 
        $h_dbPath{$dbName} = $dbPath;
    }
    close IN;
    
    our $REF_H_DBPATH = \%h_dbPath;
}


1;
