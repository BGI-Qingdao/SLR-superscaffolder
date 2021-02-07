/**********************************************************
 *
 * @Brief   :
 *      Generate final scaffolding sequnces and AGP format
 *      data by scaff_info format data and original contig
 *      sequences.
 *
 * ********************************************************/
#include "utils/args/argsparser.h"
#include "utils/log/log.h"
#include "utils/files/file_writer.h"
#include "utils/files/file_reader.h"
#include "utils/misc/Error.h"
#include "utils/agp/agp.h"
#include "utils/misc/fileName.h"
#include "stLFR/ScaffInfo.h"
#include "utils/misc/TagId.h"

#include "utils/misc/flags.h"

#include <map>
#include <set>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>

// Convert the informations if one scaffold into an array of AGPItems
struct Scaff2AGPItem
{
    public:
        void InitName(const std::string & n);
        // sbegin && send is 1base index
        void AddSeq(const std::string & sname ,
                int sbegin ,
                int send , 
                char orientation );
        void AddN(int n_size);

        const std::vector<BGIQD::AGP::AGP_Item> & Items() const { return data ; }
    private:
        long long length ;
        int part_number ;
        std::string scaff_name ;
        std::vector<BGIQD::AGP::AGP_Item> data;
};

void Scaff2AGPItem::InitName(const std::string & n)
{
    scaff_name = n ;
    length = 1 ;
    part_number = 0 ;
}

void Scaff2AGPItem::AddSeq(
        const std::string & sname ,
        int sbegin ,
        int send , 
        char orientation
        )
{
    long long add_length = send - sbegin +1 ;
    BGIQD::AGP::AGP_Item tmp ;
    tmp.object = scaff_name ;
    tmp.object_beg = length ;
    tmp.object_end = length + add_length -1;
    tmp.part_number = ++part_number ;
    tmp.component_type = BGIQD::AGP::AGP_Item::ComponentType::W ;
    tmp.lefta.component_id = sname ;
    tmp.lefta.orientation = orientation ;
    tmp.lefta.component_beg = sbegin ;
    tmp.lefta.component_end = send ;

    length += add_length ;
    data.emplace_back(std::move(tmp));
}

void Scaff2AGPItem::AddN(int n_size)
{
    BGIQD::AGP::AGP_Item tmp ;
    tmp.object = scaff_name ;
    tmp.object_beg = length;
    tmp.object_end = length + n_size -1;
    tmp.part_number = ++part_number ;
    tmp.component_type = BGIQD::AGP::AGP_Item::ComponentType::N ;
    tmp.leftb.gap_length = n_size ;
    tmp.leftb.gap_type = "scaffold";
    tmp.leftb.linkage = true ;
    tmp.leftb.linkage_evidence = "map";

    length += n_size;
    data.emplace_back(std::move(tmp));
}
// Print sequnce in the same width.
std::string blockSeq(const std::string & atcgs , int weight );
// reverse complete
std::string seqCompleteReverse(const std::string & line);

// Sequence struct.
// print seq in uniform width
struct Aseq
{
    std::string atcgs;

    std::string Seq(int weight = -1 ) const ;

    std::string ReverseCompleteSeq( int weight = -1 ) const ;

    int Len() const { return atcgs.size() ; }

    void AddPartSeq( const std::string & line ) { atcgs += line ; }

    void Reset() { atcgs.clear() ; }
};
std::string seqCompleteReverse(const std::string & line)
{
    std::string ret ;
    ret.resize(line.size(),'N');
    int index = 0;
    for( auto i = line.rbegin() ; i!= line.rend() ; i++)
    {
        if( *i == 'A' || *i == 'a' )
            ret[index++] = 'T';
        else if( *i == 'G' || *i == 'g' )
            ret[index++] = 'C';
        else if( *i == 'C' || *i == 'c' )
            ret[index++] = 'G';
        else if( *i == 'T' || *i == 't' )
            ret[index++] = 'A';
        else
            ret[index++] = 'N';
    }
    return ret;
}

std::string blockSeq(const std::string & atcgs , int weight )
{
    if( weight < 1 )
        return atcgs ;
    else
    {
        std::ostringstream ost;
        int i = 1 ;
        for( char c : atcgs)
        {
            ost<<c;
            if( i % weight == 0 || i == (int)atcgs.size() )
                ost<<'\n';
            i++ ;
        }
        return ost.str();
    }
}
std::string Aseq::ReverseCompleteSeq(int weight) const
{
    std::string r_atcgs = seqCompleteReverse(atcgs);
    if( weight < 1 )
        return r_atcgs ;
    else
    {
        std::ostringstream ost;
        int i = 1 ;
        for( char c : r_atcgs)
        {
            ost<<c;
            if( i % weight == 0 || i == (int)r_atcgs.size() )
                ost<<'\n';
            i++ ;
        }
        return ost.str();
    }
}

std::string Aseq::Seq(int weight ) const
{
    if( weight < 1 )
        return atcgs ;
    else
    {
        std::ostringstream ost;
        int i = 1 ;
        for( char c : atcgs)
        {
            ost<<c;
            if( i % weight == 0 || i == (int)atcgs.size() )
                ost<<'\n';
            i++ ;
        }
        return ost.str();
    }
}
// SOAPdenovo2 fasta header
struct SOAP2ContigHead 
{
    unsigned int contigId ;

    int len ;

    float cov;

    int is_tip ;


    void Reset(){ contigId = 0 ; is_tip = 0 ; cov = 0 ; len = 0;  } ;
    void Init( const std::string & line ) 
    {
        //>21 length 64 cvg_0.0_tip_0
        sscanf(line.c_str() 
                ,">%u length %d cvg_%f_tip_%d"
                ,&contigId
                ,&len
                ,&cov 
                ,&is_tip);

    }

    std::string Head() const {
        std::ostringstream ost;
        ost<<'>'<<contigId
            <<" length "<<len
            <<" cvg_"<<cov
            <<"_tip_"<<is_tip;
        return ost.str();
    }
};
// define SOAPdenovo2 format fasta
struct Fasta
{
    typedef SOAP2ContigHead Header;

    FLAGS_INT ;

    ADD_A_FLAG(1,UnSet);
    ADD_A_FLAG(2,Set_head);
    ADD_A_FLAG(3,Set_seq);

    void Reset() 
    {
        head.Reset();
        seq.Reset();
        flags = 0 ;
        Set_UnSet();
    }

    Header head;
    Aseq seq;
    void AddHead(const std::string & line)
    {
        Clean_UnSet();
        Set_Set_head();
        head.Init(line);
    }

    void AddSeq(const std::string & line )
    {
        if( Is_UnSet() || ! Is_Set_head() )
        {
            assert(0);
        }
        Set_Set_seq();
        seq.AddPartSeq(line);
    }
    bool Is_Setted() const { 
        return  (!Is_UnSet() )
            && Is_Set_head() 
            && Is_Set_seq() ;
    }
};

// wrap basic functions to load and parse soapdenovo2 format fasta file.
struct FastaReader
{
    static bool IsHead(const std::string & line)
    {
        return ( ! line.empty()) && line[0] == '>' ;
    }

    static void LoadAllFasta( std::istream & ist , std::vector<Fasta> & buffer )
    {
        std::string line ;
        Fasta fa;
        fa.Reset();
        while ( ! std::getline(ist,line).eof() )
        {
            if( IsHead(line) )
            {
                if(fa.Is_Setted())
                {
                    buffer.push_back(fa);
                }
                fa.Reset();
                fa.AddHead(line);
            }
            else
            {
                fa.AddSeq(line);
            }
        }
        if(fa.Is_Setted())
        {
            buffer.push_back(fa);
        }
        fa.Reset();
    }

    static bool LoadNextFasta(std::istream & ist , Fasta & fa)
    {
        std::string line ;
        fa.Reset();
        while( ! std::getline(ist,line).eof() )
        {
            if( IsHead(line) )
            {
                fa.AddHead(line);
                break ;
            }
        }
        if( ist.eof() || ! fa.Is_Set_head() )
            return false ;

        while( ! std::getline(ist,line).eof() )
        {
            if( IsHead(line) )
            {
                break ;
            }
            else
            {
                fa.AddSeq(line);
            }
        }
        // Put the head line back into istream
        if ( ! ist.eof() )
        {
            ist.rdbuf()->sputbackc('\n');
            for( auto  i = line.rbegin() ; i!= line.rend() ; i++ )
            {
                ist.rdbuf()->sputbackc(*i);
            }
        }

        return  fa.Is_Setted();
    }
};
//
// Struct to wrap all global variables and functions
//
struct AppConfig
{
    BGIQD::MISC::FileNames fNames;

    BGIQD::LOG::logger loger;

    BGIQD::stLFR::ScaffInfoHelper scaff_helper;

    typedef SOAP2ContigHead Header;
    typedef Fasta ContigFasta;
    typedef FastaReader Reader;

    std::map<unsigned int , ContigFasta> contigs;

    std::set<unsigned int> used;

    void Init( const std::string & prefix )
    {
        fNames.Init(prefix);
        loger.Init("ScaffInfo2Seqs");
    }

    void LoadScaffInfos()
    {
        auto  in = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(fNames.scaff_infos());
        if( in == NULL )
            FATAL("failed to open xxx.scaff_infos to read");

        scaff_helper.LoadAllScaff(*in);
        delete in ;
    }

    void LoadAllFasta()
    {
        auto in = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(fNames.contig());
        if( in == NULL )
            FATAL("failed to open xxx.contig to read");

        Reader reader;
        ContigFasta tmp ;
        while( reader.LoadNextFasta(*in , tmp) )
        {
            contigs[tmp.head.contigId] = tmp ;
        }
        delete in ;
    }


    void PrintScaffSeqs()
    {
        auto out = BGIQD::FILES::FileWriterFactory
            ::GenerateWriterFromFileName(fNames.scaff_seqs()) ;

        auto get_atcg = [&] ( const BGIQD::stLFR::ContigDetail & detail
                , Scaff2AGPItem & s2a ) -> std::string
        {
            std::string str = contigs.at(detail.contig_id).seq.atcgs ;
            if( ! detail.orientation )
                str = seqCompleteReverse(str);
            if( detail.gap_size > 0 )
            {
                s2a.AddSeq(
                        contig_name_cache.Id(detail.contig_id)
                        , 1
                        , detail.contig_len 
                        , (detail.orientation ? '+' : '-' )
                        ) ;

                int gap_size = detail.gap_size > min_n ?  detail.gap_size : min_n ;
                str += std::string(gap_size,'N');
                s2a.AddN(gap_size);
            }
            else if ( detail.gap_size < 0 )
            {
                if(  (int)str.length() +  (int)std::abs(detail.gap_size) > 0 )
                {
                    str = str.substr(0,str.length() + detail.gap_size);
                    s2a.AddSeq(
                            contig_name_cache.Id(detail.contig_id)
                            , 1
                            , detail.contig_len + detail.gap_size
                            , (detail.orientation ? '+' : '-' )
                            ) ;
                }
                str += std::string(min_c,'N');
                s2a.AddN(min_c);
            }
            else
            {
                // gap size = 0
                s2a.AddSeq(
                        contig_name_cache.Id(detail.contig_id)
                        , 1
                        , detail.contig_len 
                        , (detail.orientation ? '+' : '-' )
                        ) ;
                ;
            }
            return str ;
        };

        for( const auto & pair : scaff_helper.all_scaff)
        {
            Scaff2AGPItem s2a;
            s2a.InitName("scaffold_" + std::to_string(pair.first));
            (*out)<<">scaffold_"<<pair.first<<'\n';
            std::string str ;
            for( const auto & i : pair.second.a_scaff) 
            {
                used.insert(i.contig_id);
                str+= get_atcg(i,s2a);
            }
            agp_cache.data.insert(agp_cache.data.end() ,s2a.Items().begin() ,s2a.Items().end());
            (*out)<<blockSeq(str,100);
        }
        for( const auto & pair : contigs )
        {
            if( used.find( pair.first ) == used.end () )
            {
                Scaff2AGPItem tmp ;
                tmp.InitName(contig_name_cache.Id(pair.second.head.contigId));
                tmp.AddSeq(contig_name_cache.Id(pair.second.head.contigId),
                        1,
                        pair.second.head.len,
                        '+');
                if( contig_name_cache.HasId(pair.second.head.contigId) )
                    (*out)<<'>'<<contig_name_cache.Id(pair.second.head.contigId) <<'\n';
                else
                    (*out)<<pair.second.head.Head()<<'\n';
                (*out)<<pair.second.seq.Seq(100);
                agp_cache.data.push_back(*tmp.Items().begin());
            }
        }
        delete out ;
    }
    int min_n ;
    int min_c ;

    BGIQD::AGP::AGPFile agp_cache;
    void PrintAGP()
    {
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fNames.scaff_agp());
        if( out == NULL )
            FATAL( " failed to open xxx.scaff_agp to write !!! ");
        agp_cache.Print(*out);
        delete out ;
    }

    BGIQD::MISC::IdStringCache contig_name_cache;
    void LoadFakeSOAPNameMap()
    {
        auto fake = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName( "fakesoap.name2index.map.txt" );
        if( fake == NULL )
            INFO_RETURN(" no fakesoap.name2index.map.txt file !!! ");
        delete fake ;
        contig_name_cache.preload = true ;
        contig_name_cache.LoadStringIdCache("fakesoap.name2index.map.txt");
    }
} config;


int main(int argc , char ** argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string, prefix ,"prefix of files. Input xxx.scaff_infos; Output xxx.scaff_seqs");
    DEFINE_ARG_OPTIONAL(int, min_n,"min N size for gap in scaffold if not filled ","11");
    DEFINE_ARG_OPTIONAL(int, min_c,"min N size for contig overlap cut ","11");
    END_PARSE_ARGS;
    config.min_c = min_c.to_int() ;
    config.min_n = min_n.to_int();
    config.Init( prefix.to_string() );
    BGIQD::LOG::timer t(config.loger,"ScaffInfo2Seqs");

    config.LoadFakeSOAPNameMap();
    config.LoadScaffInfos();
    config.LoadAllFasta();
    config.PrintScaffSeqs();
    config.PrintAGP();
    return 0;
}
