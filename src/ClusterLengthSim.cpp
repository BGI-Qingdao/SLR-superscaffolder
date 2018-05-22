#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include <vector>
#include <map>
#include <sstream>
#include <queue>
#include <tuple>

int main(int argc , char ** argv)
{

    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string , cluster, "seed contig linear info");
        DEFINE_ARG_REQUIRED(int  , step,"step length ");
    END_PARSE_ARGS

    // Load cluster first
    std::string line ;
    std::map< unsigned int , std::map<unsigned int , float > > connections;
    unsigned int contigId , to;
    float cov ;
    auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName( cluster.to_string() );
    while( ! std::getline( *in , line ).eof() )
    {
        std::istringstream ist(line);
        ist>>contigId;

        while(! ist.eof() )
        {
            ist>>to>>cov;
            connections[contigId][to] = cov;
            connections[to][contigId] = cov ;
        }
    }
    delete in ;

    std::queue< std::tuple< unsigned int , int , int , int > > pres;
    unsigned int id ;
    int len ; 
    int start , end ;
    std::cout<<step.to_int()<<std::endl;
    // process all seed contigs
    while( ! std::getline( std::cin , line ).eof() )
    {
        std::istringstream ist(line);
        ist>>id>>len>>start>>end;
        if( (int)pres.size() >= step.to_int() )
        {
            unsigned int id1;
            int len1 ; 
            int start1 , end1 ;
            std::tie( id1 , len1 , start1 , end1 ) = pres.front() ;
            pres.pop();
            if( connections.find( id ) != connections.end() && connections.at(id).find(id1) != connections.at(id).end() )
                std::cout<<start  - end1 <<"\t"<< connections.at(id).at(id1)<<std::endl;
        }
        pres.push( std::make_tuple( id , len ,start ,end ) );
    }
}
