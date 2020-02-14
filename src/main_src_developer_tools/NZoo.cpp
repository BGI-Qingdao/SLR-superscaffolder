#include "biocommon/fasta/fasta.h"
#include <iostream>
#include <vector>
#include <regex>

typedef BGIQD::FASTA::Id_Desc_Head Header;
typedef BGIQD::FASTA::Fasta<Header> Fasta;
typedef BGIQD::FASTA::FastaReader<Fasta> Reader ;

void printNBed(const Fasta & fa ){
    const std::string & seq = fa.seq.atcgs;
    bool N_start =false ;
    int prev = 0 ;
    int index = -1 ;
    for( char c: seq ){
        index ++ ;
        if(c!= 'n' && c!='N'){
            if(N_start) {
                N_start = false ;
                std::cout<<fa.head.Head()
                    <<'\t'<<prev
                    <<'\t'<<(index)<<'\n';
            }
        } else {
            if( ! N_start ){
                N_start = true ;
                prev = index ;
            }
        }
    }
    if(N_start) {
        N_start = false ;
        std::cout<<fa.head.Head()
            <<'\t'<<prev
            <<'\t'<<(index-1)<<'\n';
    }
}

int main()
{
    std::vector<Fasta> buffer ;
    Reader reader;
    reader.LoadAllFasta(std::cin , buffer);
    for( const auto & fa : buffer) 
        printNBed(fa);
    return 0 ;
}
