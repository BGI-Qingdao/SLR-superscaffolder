#!/usr/bin/env perl

use Getopt::Long;
Getopt::Long::Configure("bundling");
use strict;

###############################################################################
#
# Usage
#
###############################################################################
=pod
Usage : bin_barcode.pl -i barcode_on_ref_file -b bin_size -o out_file_prefix.
=cut

my $Usage=
"Usage :bin_barcode.pl -i barcode_on_ref_file -b bin_size -o out_file_prefix.";

###############################################################################
#
# Parse arguments
#
###############################################################################

our $barcode_on_ref="";
our $bin_size="";
our $out_prefix="";

GetOptions(
        'in|i=s' => \$barcode_on_ref, 
        'bin|b=s' => \$bin_size, 
        'out|o=s' => \$out_prefix
        ) or die $!;
printf "Arguments : --in %s --bin %s --out %s\n",
                                   $barcode_on_ref ,$bin_size,$out_prefix ;
if ( $barcode_on_ref eq "" || int($bin_size) < 1 || $out_prefix eq "" )
{
    die($Usage);
}
###############################################################################
#
# Logic 
#
###############################################################################

open ( IN , "<".$barcode_on_ref ) or die ( "Failed to open input file" );

printf "Open file $barcode_on_ref succ. \n";
printf "Set bin size as $bin_size . \n";

our %bin_map;
our %barcode_index;
our $barcode_count = 1;
our $bin_max;
#open ( OUT1 , ">d1");
while ( <IN> )
{
    chomp;
    my @lines=split(/\t/,$_);
    my $pos=@lines[0];
    my $curr_bin = int( int($pos) / int($bin_size)) ;
    if ( $curr_bin > $bin_max ) 
    {
        $bin_max = $curr_bin;
    }
    my @barcodes = split(/\|/,$lines[1]);
  #  print @barcodes;
    foreach my $barcode ( @barcodes ) 
    {
        $barcode =~ s/^\s+//;
        $barcode =~ s/\s+$//;
        #print OUT1 "$barcode\n";
        if ( $barcode eq "0") 
        {
            next;
        }
        if ( ! exists ( $barcode_index{$barcode} ) )
        {
            $barcode_index{ $barcode} = $barcode_count ;
            $barcode_count +=1 ;
        }
        my $barcode_id = $barcode_index{$barcode}  ;
        if ( ! exists ( $bin_map{$curr_bin}{$barcode_id} ) )
        {
            $bin_map{$curr_bin}{$barcode_id} = 1;
        }
        else
        {
            $bin_map{$curr_bin}{$barcode_id} += 1;
        }
    }
}
close IN;
#for my $base ( keys %barcode_index )
#{
#    print "$base : $barcode_index{$base}\n";
#}
open( OUT , ">".$out_prefix);
printf "Open $out_prefix for write succ\n";
our $bin_index=0;
while ( $bin_index <= $bin_max )  
{
    
    if ( exists( $bin_map{$bin_index} ) )
    {
        my $size = keys %{$bin_map{$bin_index}};
        print OUT "$bin_index,$size\t"; 
        my %bin_hash = %{$bin_map{$bin_index}};
        for my $barcode ( sort keys( %bin_hash ) )
        {
            print OUT "$barcode,$bin_map{$bin_index}{$barcode}\t";
        }
    }
    else
    {
        print OUT "$bin_index,0"; 
    }
    print OUT "\n";
    $bin_index +=1;
}

