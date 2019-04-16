#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string>
int main (int argc, char ** argv) {
    pid_t pId = fork();
    if (pId == -1) {
        perror("fork error");
        exit(EXIT_FAILURE);
    } else if (pId == 0) {
        // Run command 
        std::string command;
        for( int i = 1 ; i < argc ; i++ )
        {
            command +=( std::string(argv[i]) +" ");
        }
        system(command.c_str());
        return 0 ;
    }
    printf("Parent:SelfID=%d MyChildPID=%d \n", getpid(), pId);
    do{
        sleep(1);
    }while (1);
    return 0;
}
