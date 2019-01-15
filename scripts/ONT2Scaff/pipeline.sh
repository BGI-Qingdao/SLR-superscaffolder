#!/bin/bash

##################################################
#
# MODIFY  BELOW CONFIG .
#
##################################################

OUT_PREFIX="test"
# you must give a ont read.
INPUT_ONT_READ_FA=""
#INPUT_ONT_READ_FQ=""

# USE_SCAFF_SEQ = "yes" or "no". 
#   * yes use INPUT_SCAFFOLD_SEQ 
#   * no  use INPUT_CONTIG && INPUT_SCAFF_INFOS
USE_SCAFF_SEQ="yes"
INPUT_SCAFFOLD_SEQ=""
#INPUT_CONTIG=""
#INPUT_SCAFF_INFOS=""

# USE_PA F= "yes" or "no". 
#   * yes use INPUT_PAF
#   * no  use MINMAP 2 align.
USE_PAF="no"
#INPUT_PAF=""

BIN_DIR="/home/guolidong/stLFR/debug/stLFR_Assembler/src/bin/"
MINMAP="/ldfssz1/ST_OCEAN/USER/xumengyang/software/minimap2/minimap2"
CPU=30


USE_ORDER="no"
USE_ORIENTATION="no"
USE_GAP_FILLER="no"
USE_GAP_CORRETOR="yes"

GAP_MODE=1  # 1 for shortest ; 3 for median ; 2 for random

######################################################
#
# DO NOT MODIFY BELOW CODE UNLESS YOU KNOW THE LOGIC .
#
######################################################

## prepare check ...
# check ont read
if [[  $INPUT_ONT_READ_FA == "" && $INPUT_ONT_READ_FQ == "" ]] ; then 
    echo "Please give a ont read. exit ..."
    exit
fi
if [[  $INPUT_ONT_READ_FA == ""  ]] ; then 
    if [[ ! -e $INPUT_ONT_READ_FQ ]] ; then 
        echo "$INPUT_ONT_READ_FQ is not exsist exit ..."
        exit
    fi
    FQ="yes"
else
    if [[ ! -e $INPUT_ONT_READ_FA ]] ; then 
        echo "$INPUT_ONT_READ_FA is not exsist exit ..."
        exit
    fi
    FQ="no"
fi
# check scaff*
if [[ $USE_SCAFF_SEQ ]] ; then  
    if [[ ! -e $INPUT_SCAFFOLD_SEQ ]] ; then 
        echo "$INPUT_SCAFFOLD_SEQ is not exsist exit ..."
        exit
    fi
else
    if [[ ! -e $INPUT_CONTIG || ! -e $INPUT_SCAFF_INFOS ]] ; then 
        echo "$INPUT_CONTIG or $INPUT_SCAFF_INFOS is not exsist exit ..."
        exit
    fi
fi
# check paf or minmap2
if [[ $USE_PAF ]] ; then 
    if [[ ! -e $INPUT_PAF ]] ; then 
        echo "$INPUT_PAF is not exsist exit ..."
        exit
    fi
else
    if [[ ! -x $MINMAP ]] ; then 
        echo "$MINMAP is not exsist or executable exit ..."
        exit
    fi
fi


## step 0 . split seqs if needed
if [[ $USE_SCAFF_SEQ == "yes" ]] ; then 
    $BIN_DIR/SplitSeq <$INPUT_SCAFFOLD_SEQ  >$OUT_PREFIX.contig 2>$OUT_PREFIX.scaff_infos
else
    ln -s $INPUT_CONTIG $OUT_PREFIX.contig
    ln -s $INPUT_SCAFF_INFOS $OUT_PREFIX.scaff_infos
fi

## step 1 . aligned if needed
if [[ $USE_PAF == "yes" ]] ; then 
    ln -s INPUT_PAF $OUT_PREFIX.paf
else
    if [[ $FQ == "yes" ]] ; then 
         $MINMAP -x ava-ont -t $CPU $INPUT_ONT_READ_FQ $OUT_PREFIX.contig  1>$OUT_PREFIX.paf 2>log_minmap2
    else
         $MINMAP -x ava-ont -t $CPU $INPUT_ONT_READ_FA $OUT_PREFIX.contig  1>$OUT_PREFIX.paf 2>log_minmap2
    fi
fi

## step 2. ONT process

mv $OUT_PREFIX.scaff_infos $OUT_PREFIX.scaff_infos_0
ln -s $OUT_PREFIX.scaff_infos_0 $OUT_PREFIX.scaff_infos
if [[ $USE_ORDER == "yes " ]] ; then 
    echo " Order is not valid now . continue ..."
fi

if [[ $USE_ORIENTATION == "yes " ]] ; then 
    echo " Order is not valid now . continue ..."
fi

if [[ $USE_GAP_CORRETOR  == "yes" ]] ; then
    $BIN_DIR/ONTGapCorrecter --contig2ont_paf $OUT_PREFIX.paf --force --work_mode $GAP_MODE <$OUT_PREFIX.scaff_infos >$OUT_PREFIX.scaff_infos_3
    rm $OUT_PREFIX.scaff_infos
    ln -s $OUT_PREFIX.scaff_infos_3 $OUT_PREFIX.scaff_infos
fi

if [[ $USE_GAP_FILLER  == "yes" ]] ; then
    if [[ $FQ == "yes" ]] ; then 
        $BIN_DIR/ONTGapFiller --ont_reads_fq $INPUT_ONT_READ_FQ --contig2ont_paf $OUT_PREFIX.paf --force --work_mode $GAP_MODE <$OUT_PREFIX.scaff_infos >$OUT_PREFIX.scaff_infos_4
    else
        $BIN_DIR/ONTGapFiller --ont_reads_fq $INPUT_ONT_READ_FA --contig2ont_paf $OUT_PREFIX.paf --force --work_mode $GAP_MODE <$OUT_PREFIX.scaff_infos >$OUT_PREFIX.scaff_infos_4
    fi
    rm $OUT_PREFIX.scaff_infos
    ln -s $OUT_PREFIX.scaff_infos_4 $OUT_PREFIX.scaff_infos
fi

## step 3. re-generate seq

$BIN_DIR/ScaffInfo2Seq --prefix  $OUT_PREFIX --min_n 1 --min_c 1

