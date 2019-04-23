BEGIN{
    correct=0;
    wrong=0;
    total=0;
    refLoaded=0;
    o1=0;
    o2=0;
    order=0;
    correct1=0;
    wrong1=0;
    scaffNum = 0;
    pos=0;
}
{
    if($1 == ">scaffold")
    {
        if( refLoaded == 0 )
        {
            refLoaded = 1 ;
        }
        else
        {
            if( order > 0 )
            {
                correct1 = o1;
                wrong1 = o2 ;
            }
            else
            {
                correct1 = o2 ;
                wrong1 = o1 ;
            }
            correct += correct1;
            wrong += wrong1;
            printf("%s has %d err in %d \n",scaff,wrong1,wrong1+correct1);
        }
        o1 = 0 ;
        o2 = 0 ;
        order = 0 ;
        prev = 0;
        scaffNum+=1;
        scaff=$1;
    }
    else
    {
        if ( refLoaded == 0 )
        {
            pos+=1
            posArray[$1]=pos;
            oArray[$1]=$2;
        }
        else
        {
            if( prev != 0 )
            {
                if( posArray[$1] > prev )
                    order +=1 ;
                else
                    order -= 1;
            }
            prev = posArray[$1];
            if( oArray[$1] == $2 )
                o1 +=1 ;
            else
                o2 +=1 ;
        }
    }
}
END{
    total = correct + wrong ;
    printf (" %d wrong orientation in total %d seeds of %d scaffold\n"\
            ,wrong \
            , total\
            , scaffNum \
            );
}
