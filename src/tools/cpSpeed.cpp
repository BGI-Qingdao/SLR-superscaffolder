#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <sys/time.h>
#include <errno.h>
#include <string.h>

char * buffer;
size_t read_block_size = 32768  ; 
size_t block_size = 1024 * 1024 * 1024;
size_t big_size = (1024UL * 1024 *1024 * 10) ; 
size_t true_read_size = 0 ;
char file_name[1024*10];
int no_cache = 1 ;

void resize()
{
    if( buffer == NULL )
    {
        buffer = (char*)malloc(block_size);
    }
    else
    {
        if( block_size < big_size )
        {
            block_size *= 2 ;
        }
        else
        {
            block_size += big_size ;
        }
        buffer = (char*)realloc(buffer, block_size);
    }
}

int main(int argc , char ** argv)
{
    int file_seted = 0;
    int opt ;
    while( (opt = getopt(argc,argv,"cb:f:")) != -1 )
    {
        switch(opt)
        {
            case 'c':
                no_cache = 0;
                break ;
            case 'b':
                read_block_size = atoi(optarg);
                {
                    if( read_block_size < 0 )
                        read_block_size = 1 ;
                    else if ( read_block_size > block_size) 
                        read_block_size = block_size ;
                }
                break;
            case 'f':
                strcpy(file_name,optarg);
                file_seted = 1;
                break;
            default:
                ;
        }

    }

    if( file_seted == 0 )
    {
        fprintf(stderr,"Usage : %s [-c] [-b xxx ] -f file\n",argv[0]);
        fprintf(stderr,"                    -c cache file contents or not . default no. if -c is seted , file contents will be stored in memory. otherwise be droped\n");
        fprintf(stderr,"                    -b the read block size . default is 32768. can be [ 1 , 1000000000 ]\n");
        fprintf(stderr,"                    -f file to read\n");
        exit(1);
    }
    fprintf(stderr,"Now reading %s ...\n",file_name);
    fprintf(stderr,"Block size of reading %lu ...\n",read_block_size);
    struct timeval wall_start , wall_end;
    gettimeofday(&wall_start,NULL);
    clock_t cpu_start = clock();
    resize() ;
    FILE *file = fopen(file_name,"r");
    if( file == NULL )
    {
        fprintf(stderr,"Failed to open file  %s for read . Exit ... \n",file_name);
        fprintf(stderr,"errno is : %2d\t%s\n",errno,strerror(errno));
        exit(1);
    }
    char * curr_buffer = buffer ;
    while(1)
    {
        size_t rec = fread( curr_buffer , sizeof(char) , read_block_size , file );
        if( rec == 0 )
            break ;
        if( rec > read_block_size )
        {
            fprintf(stderr,"Read file %s err happends . Exit ... \n",argv[1]);
            fprintf(stderr,"errno is : %2d\t%s\n",errno,strerror(errno));
            exit(1);
        }
        if( no_cache  )
            continue ;
        true_read_size += rec ;
        if( true_read_size + read_block_size > block_size )
            resize() ;
        curr_buffer = buffer + true_read_size ;
    }
    gettimeofday(&wall_end,NULL);
    clock_t cpu_end= clock();
    fprintf( stderr , "Done , total %lu bytes of data loaded !\n\
            %lu of memory used !\n\
            %f seconds of cpu time used!\n\
            %lu seconds of wall time used!\n",true_read_size , block_size 
            , (cpu_end - cpu_start) * 1.0f / CLOCKS_PER_SEC
            , (wall_end.tv_sec - wall_start.tv_sec) 
           );
    return 0 ;
}


