#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use Data::Dumper;
use Getopt::Long;
use lib "$Bin/perl5";
use MyModule::GlobalVar qw($REF_H_DBPATH);


#######################
##USAGE
#######################

sub usage{
    print STDERR <<USAGE;
Version:1.0
2018-06-8       zengpeng\@genomics.cn
ReVersion:1.01
Fri Jun  8 11:21:38 HKT 2018
        Usage: perl $0
        where: 
        Options
            -i1  <c>: fq1 (must be gave)
            -i2  <c>: fq2 (must be gave)
            -o  <c>: dir for output
            -b  <c>: barcode list of two columns:sequence and barcodeID
            -r  <n>: barcode regions: 101_10,117_10,145_10 or 101_10,117_10,133_10
            -n  <c>: 
            -h     : show this help message
USAGE
}

my %opt;
GetOptions(\%opt,"i1=s","i2=s","o=s","b=s","r=s","n=s","h");

#########check parameters and write into memory ######
if(exists $opt{h}){
    usage;
    exit;
}
if(!$opt{i1} || !$opt{i2}){
    usage;
    print STDERR "\n\t\t\tERROR: -i1 or -i2 option missing!\n";
    exit;
}

my $temp=`pwd`;  chomp $temp ; 
my %DATABASE = %{$REF_H_DBPATH};

$opt{o} ||=$temp;
$opt{o}   ="$temp/$opt{o}" if($opt{o} !~/^\//);
$opt{b} ||= $DATABASE{"barcodelist"};
$opt{r} ||= "101_10,117_10,145_10";
$opt{n} ||= 1;

mkdir  $opt{o} unless (-d $opt{o});

my @barcode_loci;
foreach my $re(split(/\,/,$opt{r})){
    my ($st,$le)=split(/\_/,$re);
    $st-=1;
    push @barcode_loci,[$st,$le];
}
@barcode_loci=sort {$a->[0]<=>$b->[0]} @barcode_loci;

my $read2len=$barcode_loci[0][0];

##############what you want to do #################


my %barcode_hash;
open B,$opt{b} || die "can't find barcode file\n";
my $n=0;
while(<B>){
    chomp;
    $n ++;
    my @f=split;
    my @arr=barcode_1mismatch($f[0]);
    foreach my $aa(@arr){
        if($barcode_hash{$aa}){
            print $aa,"\t",$barcode_hash{$aa},"\t",$f[1],"\n" if($barcode_hash{$aa} > $f[1]);
        }
        $barcode_hash{$aa}=$f[1];
    }
}
close B;

#exit;

my $barcode_types = $n * $n *$n;
my $barcode_each = $n;

open (F1,"gzip -dc $opt{i1} | ") || die "the Input Can't open $!";
open (F2,"gzip -dc $opt{i2} | ") || die "the Input Can't open $!";

open O1," | gzip > $opt{o}/split_read.1.fq.gz";
open O2," | gzip > $opt{o}/split_read.2.fq.gz";


my $reads_num;
my $progress;
my %index_hash;
my %index_hash_reverse;
my $split_barcode_num;
my $T;

my @line;
my @Read_num;
$Read_num[0] = 0;
my $split_reads_num;


while(<F1>){
    chomp ;

    $reads_num ++;

    my $id1=$_;
    my $seq1=<F1>;
    my $plus1=<F1>;
    my $qual1=<F1>;

    my $id2=<F2>;
    my $seq2=<F2>;
    my $plus2=<F2>;
    my $qual2=<F2>;

    if($id1 !~/^\@/ || $id2 !~/^\@/){
        last "Fastq error at reads $reads_num\n$id1\n$id2";
    }

    if($id1 !~/\/1/ || $id2 !~/\/2/){
        last "Reads1 and Reads2 ID error in /1 or /2 at reads $reads_num\n";
    }

    my $le=length $id1;
    $le-=2;
    if(substr($id1,0,$le) ne substr($id2,0,$le)){
       last "Fastq reads ID unequal at reads $reads_num\n$id1\n$id2";
    }

    my ($b1,$b2,$b3)=cut_barcode($seq2);

    if((exists $barcode_hash{$b1}) && (exists $barcode_hash{$b2}) && (exists $barcode_hash{$b3})){
        my $hash = $barcode_hash{$b1}."_".$barcode_hash{$b2}."_".$barcode_hash{$b3};
        if(!(exists $index_hash{$hash})){
            $split_barcode_num ++;
            $index_hash{$hash} = $split_barcode_num;
            $index_hash_reverse{$split_barcode_num} = $hash;
            $Read_num[$index_hash{$hash}] = 0;
        }

        $split_reads_num ++;
        $Read_num[$index_hash{$hash}] ++;

        my $id=$id1;
        $id=~s/\/1//;
        print O1 "$id\#$hash\/1\t$index_hash{$hash}\t1\n";
        print O1 $seq1;
        print O1 $plus1;
        print O1 $qual1;

        print O2 "$id\#$hash\/2\t$index_hash{$hash}\t1\n";
        print O2 substr($seq2,0,$read2len),"\n";
        print O2 $plus2;
        print O2 substr($qual2,0,$read2len),"\n";

    }
    else{
        $Read_num[0] ++;

        my $id=$id1;
        $id=~s/\/1//;

        print O1 "$id\#0_0_0\/1\t0\t1\n";
        print O1 $seq1;
        print O1 $plus1;
        print O1 $qual1;

        print O2 "$id\#0_0_0\/2\t0\t1\n";
        print O2 substr($seq2,0,$read2len),"\n";
        print O2 $plus2;
        print O2 substr($qual2,0,$read2len),"\n";
    }
}

open OUT2, ">$opt{o}/split_stat_read1.log" or die "Can't write file";
print OUT2 "Barcode_types = $barcode_each * $barcode_each * $barcode_each = $barcode_types\n";
my $r;
$r = 100 *  $split_barcode_num/$barcode_types;
print OUT2 "Real_Barcode_types = $split_barcode_num ($r %)\n";
$r = 100 *  $split_reads_num/$reads_num;
print OUT2 "Reads_pair_num  = $reads_num \n";
print OUT2 "Reads_pair_num(after split) = $split_reads_num ($r %)\n";
for(my $i=1;$i<=$split_barcode_num;$i++){
    print OUT2 "$i\t$Read_num[$i]\t$index_hash_reverse{$i}\n";
}


close F1;
close F2;
close O1;
close O2;

### swimming in the sky and flying in the sea ####


sub barcode_1mismatch{

    my $bb=shift;
    my @out;
    my $le=length $bb;

    my $rev_bb=reverse $bb;
    $rev_bb=~tr/ATCGatcg/TAGCtagc/;

    for(my $i=0;$i<$le;$i++){
        for my $base('A','T','C','G'){
            my $bb2=$bb;	
            substr($bb2,$i,1)=$base;
            my $bb2_rev=$rev_bb;
            substr($bb2_rev,$i,1)=$base;
            push @out,$bb2;
        }
    }

    return @out;
}

sub cut_barcode{
    my $seq=shift;
    my @out;
    for(my $i=0;$i<=$#barcode_loci;$i++){
        push @out,substr($seq,$barcode_loci[$i][0],$barcode_loci[$i][1]);
    }
    return @out;
}
