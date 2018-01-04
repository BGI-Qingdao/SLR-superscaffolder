#include "contig_barcode.h"
#include "argsparser.h"
#include "file_reader.h"
#include "file_writer.h"
#include "sam_parser.h"
#include <cassert>

using namespace BGIQD;
using namespace BGIQD::JOB01;
using namespace BGIQD::ARGS;
using namespace BGIQD::LOG;
using namespace BGIQD::FILES;
using namespace BGIQD::SAM;


typedef std::map<int , int > contigTypes;

typedef std::vector<SAM::MatchData> contigMatchData;

bool HasPecfectMatch( const contigMatchData & v)
{
    for( const auto & i : v)
    {
        bool noHS = true ;
        for( const  auto & d : i.detail.infos )
        {
            if( d.type == SAM::CIGAR::H || d.type ==SAM::CIGAR::S )
                noHS = false;
        }
        if( noHS ) 
            return true ;
    }
    return false;
}

bool ErrConnect( const contigMatchData & v, int thresold ,  int & err)
{
    contigMatchData pre , next ;
    int a , b;
    for( const auto & i : v)
    {
        if( ( i.detail.infos[0].type == CIGAR::H || i.detail.infos[0].type == CIGAR::S )&& (i.detail.infos.rbegin()->type != CIGAR::H  && i.detail.infos.rbegin()->type != CIGAR::S ))
        {
            a = i.detail.infos[0].end_position_on_read - i.detail.infos[0].start_position_on_read ;
            pre.push_back(i);
        }
        else if ( (i.detail.infos[0].type != CIGAR::H && i.detail.infos[0].type != CIGAR::S) &&( i.detail.infos.rbegin()->type == CIGAR::H  || i.detail.infos.rbegin()->type == CIGAR::S ))
        {
            b = i.detail.infos.rbegin()->end_position_on_read - i.detail.infos.rbegin()->start_position_on_read ;
            next.push_back(i);
        }
    }
    if(pre.size() <1 || next.size() < 1)
        return false ;
    for(const auto & p : pre )
    {
        for( const auto n : next )
        {
            if( abs ((int)p.first_match_position - (int)n.first_match_position) <= thresold )
                return false ;
        }
    }
    err = (a < b ? a : b ) +1;
    if ( err < 6 ) 
        return false;
    return true;
}

void loadContigTypes( const std::string & file , contigTypes & types )
{
    timer t(log1,"loadContigTypes");
    auto in = FileReaderFactory::GenerateReaderFromFileName(file);
    int error_contig  = 0;
    int err_len = 0;
    size_t errcontig_len = 0;
    contigMatchData contigCache;

    float ttt = 4.0f;
    while( ! in->eof() )
    {
        std::string line;
        std::getline(*in,line);
        if( in->eof() )
            break;
        LineParser p(line);
        if(p.IsHead())
            continue;
        auto d0 = p.ParseAsMatchData();
        int read = std::stoi(d0.read_name);
        auto itr = types.find(read);
        if (d0.ref_name != "chr19" )
        {
            assert ( d0.ref_name == "*" );
            assert ( itr == types.end() );
            types[read] = 0 ; // means unmatched !
        }
        else
        {
            if( itr != types.end() )
            {
                assert( itr->second != 0);
                itr->second ++ ;
            }
            else
            {
                int err ;
                // deal last contig first
                if( contigCache.size() >1 && ! HasPecfectMatch(contigCache) && ErrConnect(contigCache , contigCache.at(0).read_len * ttt , err) )
                {
                    error_contig ++ ;
                    errcontig_len += contigCache.at(0).read_len ;
                    log1<<lstart()<<contigCache.at(0).read_name<<"\t"<<contigCache.at(0).read_len<<"\t"<<err<<lend();
                    err_len += err;
                    if( contigCache.size() > 2)
                    {
                        types[read] = 2;
                    }
                    else
                    {
                        types[read] = 1 ;
                    }
                }
                // for new contig
                contigCache.clear();
                types[read] = 1;
            }
            contigCache.push_back(d0);
        }
    }
    log1<<lstart()
        <<"Total "<<types.size()<<" and err "<<error_contig
        <<" and total err len "<<errcontig_len
        <<" and err  len "<<err_len
        <<lend();
}

void printContigTypes( const std::string & file , const contigTypes & types)
{
    auto out = FileWriterFactory::GenerateWriterFromFileName(file);
    timer t(log1,"printContigTypes");
    for ( const auto & pair : types )
    {
        (*out)<<pair.first<<"\t"<<pair.second<<std::endl;
    }
    delete out;
}

int main(int argc ,char **argv)
{
    initLog("JOB05");

    START_PARSE_ARGS
    DEFINE_ARG(std::string , input , 'i');
    DEFINE_ARG(std::string , output, 'o');
    END_PARSE_ARGS

    contigTypes types;
    loadContigTypes(input.to_string() , types);
    printContigTypes(output.to_string() , types );
    return 0;
}
