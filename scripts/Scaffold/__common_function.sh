#
# this file is NOT auto-generated from __common_function.m4
#

###############################################################################
#
# check xxx_#1 xxx_#2 ...
# return the valid name
#
function get_backup_name()
{
    tag=1;
    while [[ -e $1"_#"$tag ]] ;
    do
        ((tag=$tag+1))
    done
    echo $1"_#"$tag
}

###############################################################################
#
# check each file in list and backup the file if it exsist.
#

function try_backup_list()
{
    while [[ $# -gt 0 ]] ;
    do
        if [[ -e $1 ]] ; then
            to=`get_backup_name $1`
            echo "info      : $STEP mv $1 to $to"
            mv $1 $to
        fi
        shift
    done
}

###############################################################################
#
# check each file in list
#   echo "yes" if all file exsist and valid for read.
#   echo "no" otherwise.
#
function check_file_read_list()
{
    while [[ $# -gt 0 ]] ;
    do
        if [[ !  -r $1 ]] ; then
            echo "no"
            return 1
        fi
        shift
    done
    echo 'yes'
    return 0;
}

function check_input()
{
    check_file=`check_file_read_list $@`
    if test $check_file = 'no' ; then
        echo "ERROR : $STEP check input file failed !!! exit ..."
        echo "ERROR : $STEP $@"
        echo "ERROR : $STEP please double check above files exist and readable !!!"
        exit 1
    fi
}

function check_output()
{
    check_file=`check_file_read_list $@`
    if test $check_file = 'no' ; then
        echo "ERROR : $STEP check output file failed !!! exit ..."
        echo "ERROR : $STEP $@"
        echo "ERROR : $STEP please double check above files exist and readable !!!"
        exit 1
    fi
}
