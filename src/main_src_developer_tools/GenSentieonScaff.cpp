#include "common/args/argsparser.h"
#include "biocommon/fasta/fasta.h"

#include <sstream>
#include <string>
#include <cassert>
#include <stdio.h>

struct SentieonContigHead
{
    int contigId ;

    int scaffoldId ;

    int len ;

    void Init( const std::string & line )
    {
        sscanf(line.c_str() ,"C%d\tSID=S%d\tLEN=%d",&contigId , &scaffoldId,&len);
        assert( contigId >= 0);
        assert( scaffoldId >= 0);
        assert( len > 0);
    }

    std::string Head() const { return "";} 

    void Reset(){ contigId = -1 ; scaffoldId = -1 ; len = 0;  } ;

};

typedef BGIQD::FASTA::Fasta<SentieonContigHead> Fa;
typedef BGIQD::FASTA::Fasta<BGIQD::FASTA::NormalHead> Scaff;

int main(int argc , char **argv)
{

    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(int         , n_size," the n size between each contig");
    END_PARSE_ARGS
    // Load ref 
    std::vector<Fa> AllRef;
    BGIQD::FASTA::FastaReader<Fa> Reader;
    Reader.LoadAllFasta(std::cin ,AllRef);

    //       scaff_id          contig_id
    std::map< int , std::vector< int > > AllScaff;

    if( AllRef.size() < 1 )
    {
        std::cerr<<"no ref sequence loaded !!! error !!! exit ... "<<std::endl;
        return -1 ;
    }

    for(int i = 0 ; i <(int)AllRef.size() ; i++ )
    {
        const auto & a_contig = AllRef.at(i);
        AllScaff[a_contig.head.scaffoldId].push_back(i);
    }
    for( const auto & pair : AllScaff )
    {
        Scaff tmp ;
        const auto scaff_id = pair.first ;
        const auto contig_ids = pair.second ;

        std::ostringstream head_ost;
        std::ostringstream seq_ost;

        head_ost<<">SID="<<scaff_id<<"\tcontigs:";
        for(int i = 0 ; i <(int)contig_ids.size() ; i++ )
        {
            const auto & contig = AllRef.at(contig_ids.at(i));
            head_ost<<'\t'<<contig.head.contigId;
            seq_ost<<contig.seq.atcgs;
            if( i != (int)contig_ids.size() - 1 )
                seq_ost<<std::string(n_size.to_int(),'N');
        }
        tmp.AddHead(head_ost.str());
        tmp.AddSeq(seq_ost.str());
        std::cout<<tmp.head.Head()<<std::endl;
        std::cout<<tmp.seq.Seq(100);
    }
    return 0;
}
