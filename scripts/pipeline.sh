#!/bin/bash
########################################################
#
# PRINT USAGE
#
########################################################

if [[ $# -lt 1 || $1 == "help" ||  $1 == "-h" || $1 == "--help" ]] ; then 
    echo "Usage : $0 conf step_nums"
    echo ""
    echo "example 1  , run all pipeline :"
    echo "  $0 conf 1 2 3 4 5 6 7 8 9 10 11 12 13" 
    echo "example 2  , only run mercontig step "
    echo "  $0 conf 13"
    echo ""
    echo "explain each step :"
    echo "  1.  soap pregraph"
    echo "  2.  soap contig"
    echo "  3.  use bwa index to make index of soap contig"
    echo "  4.  use bwa mem to map stLFT reads into soap contig"
    echo "  5.  covert sam file into read2contig format"
    echo "  6.  generate seed contig" 
    echo "  7.  generate barcode on contig data"
    echo "  8.  chop bin"
    echo "  9.  cluster bin"
    echo "  10. merge bin cluster data into contig cluster data"
    echo "  11. generate super contig skeleton"
    echo "  12. fill super contig skeleton into super contig"
    echo "  13. merge soap contig and super contig"
fi

########################################################
#
# check conf file
#
########################################################
function check_file_exsit()
{
    local file_2_check=$1
    if [[ ! -e $file_2_check ]] ; then 
        echo "ERROR : file $1 is invalid ! exit ... "
        exit
    fi
}

if [[ -f $1 ]] ; then 
    source $1
    if [[ ! -d $STLFR_BIN_PATH  || ! -e $STLFR_BIN_PATH"/SuperContig" ]] ; then
        echo "ERROR : $STLFR_BIN_PATH is invalid . exit ..."
        exit
    fi
    file_2_check $SOAP
    file_2_check $LIB_FILE
    file_2_check $BWA
else
    echo "ERROR : conf file $1 invalid . exit ..."
    exit
fi
shift

########################################################
#
# FILE NAMES.
# DO NOT MODIFY BELOW VALUE
#
########################################################

contig=$PREFIX".contig"
samread2contig=$PREFIX".read2contig.sam"
read2contig=$PREFIX".read2contig"
seeds=$PREFIX".seeds"
barcode2contig=$PREFIX".barcode2contig"
barcode2bin=$PREFIX".barcode2bin"
bin_cluster=$PREFIX".bin_cluster"
contig_cluster=$PREFIX".cluster"
super_contig=$PREFIX".contigroad"
super_contig_fill=$PREFIX".contigroadfill"

#######################################################
#
# RUN STEP ONE BY ONE
#
########################################################

while [[ $# -gt 0 ]] ;
do
    echo "run Step $1 ... now"
    date

    if [[ $1 -eq  1 ]] ; then 
        # Step 1. pregrap
        echo "$SOAP pregraph -s $LIB_FILE -R -o $PREFIX -K $K  -p $THREADS -d 1                    1>out_pregraph  2>&1"
              $SOAP pregraph -s $LIB_FILE -R -o $PREFIX -K $K  -p $THREADS -d 1                    1>out_pregraph  2>&1
    elif [[ $1 -eq  2 ]] ; then 
        # Step 2. contig
        echo "$SOAP contig -g $PREFIX -R -p $THREADS  -D 1                                    1>out_contig    2>&1"
              $SOAP contig -g $PREFIX -R -p $THREADS  -D 1                                    1>out_contig    2>&1
    elif [[ $1 -eq  3 ]] ; then 
        # Step 3. bwa index 
        echo "bwa index -a bwtsw $contig                                                      1>out_bwt       2>&1"
              bwa index -a bwtsw $contig                                                      1>out_bwt       2>&1
    elif [[ $1 -eq  4 ]] ; then 
        # Step 4. bwa align
        echo "$BWA mem $contig -t $THREADS -k $K                                              >$samread2contig 2>out_bwtmem"
              $BWA mem $contig -t $THREADS -k $K                                              >$samread2contig 2>out_bwtmen
    elif [[ $1 -eq  5 ]] ; then 
        # Step 5. sam to txt
        echo "$STLFR_BIN_PATH/Sam2ReadInContig -c 1000000000 -o chr19 <$samread2contig        >$read2contig 2>out_sam2readincontig"
              $STLFR_BIN_PATH/Sam2ReadInContig -c 1000000000 -o chr19 <$samread2contig        >$read2contig 2>out_sam2readincontig
    elif [[ $1 -eq  6 ]] ; then 
        # Step 6. find seeds
        echo "$STLFR_BIN_PATH/StaticsticUnique -K $K -m $BIN_SIZE <$contig                    >$seeds  2>out_staticstic"
              $STLFR_BIN_PATH/StaticsticUnique -K $K -m $BIN_SIZE <$contig                    >$seeds  2>out_staticstic
    elif [[ $1 -eq  7 ]] ; then 
        # Step 7. link barcode and contig information
        echo "$STLFR_BIN_PATH/BarcodeOnContig_NoRef -s $seeds <$read2contig                   >$barcode2contig 2>out_barcode_on_contig"
              $STLFR_BIN_PATH/BarcodeOnContig_NoRef -s $seeds <$read2contig                   >$barcode2contig 2>out_barcode_on_contig
    elif [[ $1 -eq  8 ]] ; then 
        # Step 8. chop bin
        echo "$STLFR_BIN_PATH/BarcodeOnBin -b $BIN_SIZE -i $barcode2contig -o $barcode2bin    1>out_chopbin 2>&1"
              $STLFR_BIN_PATH/BarcodeOnBin -b $BIN_SIZE -i $barcode2contig -o $barcode2bin    1>out_chopbin 2>&1
    elif [[ $1 -eq  9 ]] ; then 
        # Step 9. cluster
        echo "$STLFR_BIN_PATH/BinCluster -i $barcode2bin -t $THREADS \
                  -o $bin_cluster -s $SIMULARITY_THRESHOLD                                    >out_cluster 2>&1"
              $STLFR_BIN_PATH/BinCluster -i $barcode2bin -t $THREADS \
                  -o $bin_cluster -s $SIMULARITY_THRESHOLD                                    >out_cluster 2>&1
    elif [[ $1 -eq  10 ]] ; then 
        # Step 10. merge bin cluster to contig cluster
        echo "$STLFR_BIN_PATH/MergeClusterResult <$bin_cluster                                >$contig_cluster 2>out_mergebincluster"
              $STLFR_BIN_PATH/MergeClusterResult <$bin_cluster                                >$contig_cluster 2>out_mergebincluster
    elif [[ $1 -eq  11 ]] ; then 
        # Step 11. SuperContitg
        echo "$STLFR_BIN_PATH/SuperContig -o $PREFIX -K $K  -t $THREADS -i -s -l 10000        >out_super.txt 2>log_super"
              $STLFR_BIN_PATH/SuperContig -o $PREFIX -K $K  -t $THREADS -i -s -l 10000        >out_super.txt 2>log_super
    elif [[ $1 -eq 12 ]] ; then
        # Step 12 MergeContig
        echo "$STLFR_BIN_PATH/FillContigRoad -K $K  -t $THREADS -o $PREFIX -l 10000            >$super_contig_fill  2>log_fillsuper"
              $STLFR_BIN_PATH/FillContigRoad -K $K  -t $THREADS -o $PREFIX -l 10000            >$super_contig_fill  2>log_fillsuper
    elif [[ $1 -eq 13 ]] ; then 
        # Step 12 MergeContig
        echo "$STLFR_BIN_PATH/MergeContig -K $K -o $PREFIX -l                                 1>out_mergecontig 2>&1"
              $STLFR_BIN_PATH/MergeContig -K $K -o $PREFIX -l                                 1>out_mergecontig 2>&1
    else
        echo "WARN : invalid step num : $1 . skip it now ... "
    fi

    date
    echo "run Step $1 ... end"

    shift
done
