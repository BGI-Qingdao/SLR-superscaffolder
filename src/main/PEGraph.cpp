#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/files/file_writer.h"
#include "common/files/file_reader.h"
#include "common/string/stringtools.h"
#include "common/error/Error.h"
#include "common/stl/mapHelper.h"
#include "common/freq/freq.h"

#include "biocommon/sam_bam/sam_parser.h"

#include "algorithm/distribution/distribution.h"

#include "stLFR/contigPEGraph.h"

#include "soap2/fileName.h"

struct AppConfig
{
    typedef BGIQD::DISTRIBUTION::IntervalDistribution<int> Dist;
    typedef BGIQD::DISTRIBUTION::IntervalPercent<int> Percent;
    BGIQD::LOG::logger loger;
    BGIQD::stLFR::ContigPEGraph pe_graph;
    BGIQD::SOAP2::FileNames fName;

    std::map<unsigned int , int > contigLen_cache ;
    std::map<unsigned int , std::map<unsigned int , std::vector<int> > > pe_cache ;
    std::map<unsigned int , std::map<unsigned int , std::vector<std::vector<int>> > > pe_gap_cache ;

    int dist_bin;
    int insert_size ;
    int max_is ;
    bool ptest;
    Dist dist;
    Percent percert_all;
    void Init(const std::string & prefix)
    {
        fName.Init(prefix);
        BGIQD::LOG::logfilter::singleton().get("PEGraph",BGIQD::LOG::loglevel::INFO, loger);
        dist.Init(dist_bin,0,1000);
    }

    void LoadSeeds()
    {
        BGIQD::LOG::timer t(loger,"LoadSeeds");
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.seeds()) ;
        if( in == NULL )
            FATAL( "open .seeds file to read failed !!! " );

        std::string line ;
        while( in && !std::getline(*in, line).eof() )
        {
            auto items = BGIQD::STRING::split(line,"\t");
            contigLen_cache[std::stoul(items[0])] = std::stoul(items[1]);
            contigLen_cache[std::stoul(items[0])+1] = std::stoul(items[1]);
            pe_graph.AddNode(std::stoul(items[0]) , std::stoul(items[1]));
        }
        delete in ;
    };

    void SavePECahce(unsigned int a1 , int p1 ,int p2 ,
            unsigned int a2 , int p3 , int p4 )
    {
        if(a1< a2)
        {
            std::vector<int> tmp;
            tmp.push_back(p1);
            tmp.push_back(p2);
            tmp.push_back(p3);
            tmp.push_back(p4);
            pe_gap_cache[a1][a2].push_back(tmp);
        }
        else
        {
            std::vector<int> tmp;
            tmp.push_back(p3);
            tmp.push_back(p4);
            tmp.push_back(p1);
            tmp.push_back(p2);
            pe_gap_cache[a1][a2].push_back(tmp);
        }
    }
    //
    //
    //    +-
    //        pos 1bp of r1            pos 1bp of r2
    //    |-----+                    |---------+
    //          | -------insert-size---------  |
    //          V  r1                       r2 V
    //          |----->|                |<-----|
    //    |------------>>| ---gap--> |------------------------>>|
    //     r1_match_seed                r2_match_seed
    //    |------------>>| ---gap--> |------------------------>>|
    //        Pcontig                   Econtig
    //          |--------|           |---------|
    //             Pleft                 ELeft
    //
    //   Pleft = Pcontig.Len - (pos 1bp of r1)
    //   Eleft = (pos 1bp of r2)
    //   gap = insert_size - Pleft - ELeft
    //
    //   ++
    //        pos 1bp of r1            pos 1bp of r2
    //    |-----+                              +----------------|
    //          | -------insert-size---------  |
    //          V  r1                       r2 V
    //          |----->|                |<-----|
    //    |------------>>| ---gap--> |<<------------------------|
    //     r1_match_seed                r2_match_seed
    //    |------------>>| ---gap--> |------------------------>>|
    //        Pcontig                    Econtig
    //          |--------|           |---------|
    //             Pleft                 ELeft
    //
    //   Need think : r1 may be at P contig :
    //                      r1
    //                   |----->|
    //                  |-------------------->>| r1_match_seed
    //                  |----  gap-->|
    //                         r2
    //                       |<-----|
    //    |<<------------------------| r2_match_seed
    //                   |--insert--|
    //   
    //
    //   Pleft = Pcontig.Len - (pos 1bp of r1)
    //   Eleft = Econtig.Len - (pos 1bp of r2)
    //   gap = insert_size - Pleft - ELeft
    //   --
    //        pos 1bp of r1            pos 1bp of r2
    //          +--------|           |---------+
    //          | -------insert-size---------  |
    //          V  r1                       r2 V
    //          |----->|                |<-----|
    //    |<<------------| <--gap--- |------------------------>>|
    //     r1_match_seed                r2_match_seed
    //    |------------>>| ---gap--> |------------------------>>|
    //        Pcontig                    Econtig
    //          |--------|           |---------|
    //             Pleft                 ELeft
    //
    //   Pleft = (pos 1bp of r1)
    //   Eleft = (pos 1bp of r2)
    //   gap = insert_size - Pleft - ELeft
    //
    //
    //   -+
    //        pos 1bp of r1            pos 1bp of r2
    //          +--------|                     +----------------|
    //          | -------insert-size---------  |
    //          V  r1                       r2 V
    //          |----->|                |<-----|
    //    |<<------------| <--gap--- |<<------------------------|
    //     r1_match_seed                r2_match_seed
    //    |------------>>| ---gap--> |------------------------>>|
    //        Pcontig                    Econtig
    //          |--------|           |---------|
    //             Pleft                 ELeft
    //
    //   Pleft = (pos 1bp of r1)
    //   Eleft = Econtig.Len - (pos 1bp of r2)
    //   gap = insert_size - Pleft - ELeft
    //

    void LoadPECahce()
    {
        BGIQD::LOG::timer t(loger,"LoadPECahce");
        BGIQD::FREQ::Freq<int> ISFreq;
        std::string line ;
        {
            auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.read2contig());
            if ( in == NULL )
                FATAL(" open xxx.read2contig to read failed !!! ");
            long total_pair = 0 ;
            long total_pair_pe_same = 0 ;
            long insert_size_sum = 0;
            long insert_size_num = 0;
            // round 1 . calc mean value of insert size.
            unsigned int  prev = -1 ;

            std::vector<int> pos;

            bool isP_wait = false ;
            while( ! std::getline( *in , line ).eof() )
            {
                auto items = BGIQD::STRING::split(line,"\t");

                if( items[5] == "P"  )
                {
                    pos.clear();
                    isP_wait = false ;
                    if( items[6] == "Y" )
                    {
                        total_pair ++ ;
                        prev = std::stoul(items[1]);
                        int first_1bp = std::stoul(items[2]);
                        pos.push_back(first_1bp);
                        int last_1bp = 0 ;
                        if( items[3] == "+" )
                            last_1bp = first_1bp + 100 -1 ;
                        else
                            last_1bp = first_1bp - 100 +1 ;
                        pos.push_back(last_1bp);
                        isP_wait = true ;
                    }
                    else
                    {
                        prev = -1 ;
                        isP_wait = false ;
                    }

                }
                else
                {
                    if( items[6] == "Y" && std::stoul(items[1]) == prev && isP_wait )
                    {
                        total_pair_pe_same ++ ;
                        insert_size_num ++ ;
                        int first_1bp = std::stoul(items[2]);
                        pos.push_back(first_1bp);
                        int last_1bp = 0 ;
                        if( items[3] == "+" )
                            last_1bp = first_1bp + 100 -1 ;
                        else
                            last_1bp = first_1bp - 100 +1 ;
                        pos.push_back(last_1bp);

                        assert(pos.size() == 4 );
                        std::sort(pos.begin() ,pos.end());
                        int smallest = pos[0];
                        int biggest = pos[3] ;
                        int IS = biggest - smallest +1 ;
                        insert_size_sum += IS ;
                        dist.Count(IS);
                        ISFreq.Touch(IS);
                    }
                    isP_wait = false ;
                }
            }

            delete in ;
            percert_all = dist.CalcPercent();
            insert_size = insert_size_sum / insert_size_num ;
            loger<<BGIQD::LOG::lstart() << " average insert_size "<<insert_size<<BGIQD::LOG::lend();
            loger<<BGIQD::LOG::lstart() << " total_pair "<<total_pair<<BGIQD::LOG::lend();
            loger<<BGIQD::LOG::lstart() << " total_pair_pe_same "<<total_pair_pe_same<<BGIQD::LOG::lend();
            loger<<BGIQD::LOG::lstart() << " insert size freq is \n"<<ISFreq.ToString()<<BGIQD::LOG::lend();
            loger<<BGIQD::LOG::lstart() << " distribution of IS is \n"<<percert_all.ToString()<<BGIQD::LOG::lend();
        }

        {
            auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.read2contig());
            if ( in == NULL )
                FATAL(" open xxx.read2contig to read failed !!! ");
            long total_pair_pe_diff= 0 ;
            bool isP_wait = false ;

            int Pleft , Pright ;
            unsigned int  Pcontig , Pcontig_1 ;

            while( ! std::getline( *in , line ).eof() )
            {
                //0    1         2     3   4   5   6   7
                //2\t  140459\t  71\t  -\t 1\t E\t Y\t 0
                auto items = BGIQD::STRING::split(line,"\t");
                if( items[5] == "P"  )
                {
                    isP_wait = false ;
                    if( items[6] == "Y" )
                    {
                        int first_1bp = std::stoul(items[2]);
                        unsigned int  match_contig = std::stoul(items[1]);
                        if( contigLen_cache.find( match_contig ) == contigLen_cache.end() )
                            continue ;
                        isP_wait = true ;
                        Pleft = - first_1bp ;
                        Pright = contigLen_cache[match_contig] - first_1bp ;
                        Pcontig = match_contig ;
                        Pcontig_1  =  match_contig + 1 ;
                        if ( items[3] == "-" )
                        {
                            std::swap(Pcontig , Pcontig_1);
                            Pleft = -(contigLen_cache[match_contig] - first_1bp) ;
                            Pright = first_1bp;
                        }
                    }
                }
                else
                {
                    if( ! isP_wait )
                        continue ;
                    isP_wait = false ;
                    if( items[6] == "N" )
                        continue ;
                    unsigned int match_contig = std::stoul(items[1]);
                    if( contigLen_cache.find( match_contig ) == contigLen_cache.end() )
                        continue ;
                    if( match_contig == Pcontig || match_contig == Pcontig - 1  )
                    {
                        continue ;
                    }

                    total_pair_pe_diff ++ ;
                    int first_1bp_of_E = std::stoul(items[2]) ;
                    unsigned int Econtig = 0 ;
                    unsigned int Econtig_1 = 0 ;
                    int Eleft = 0 , Eright  = 0 ;

                    Eleft = insert_size -  first_1bp_of_E ;
                    Eright = insert_size + (contigLen_cache[match_contig] - first_1bp_of_E ) ;;
                    Econtig = match_contig ;
                    Econtig_1  = match_contig + 1 ;

                    if ( items[3] == "+" )
                    {
                        Eleft =  insert_size -  (contigLen_cache[match_contig] - first_1bp_of_E ) ;
                        Eright = insert_size + first_1bp_of_E ;
                        std::swap(Econtig , Econtig_1);
                    }

                    SavePECahce( Pcontig, Pleft , Pright , Econtig, Eleft, Eright );

                    SavePECahce( Econtig_1, Eleft, Eright ,  Pcontig_1 , Pleft , Pright );

                    if( Pleft < Eleft ) 
                    {
                        int gap = Eleft - Pright ;
                        if( ( insert_size - gap ) < max_is )
                        {
                            pe_cache[Pcontig][Econtig].push_back(gap);
                            pe_cache[Econtig_1][Pcontig_1].push_back(gap);
                        }
                    }
                    else
                    {
                        //int gap = Pleft - Eright ;
                        //pe_cache[Econtig][Pcontig].push_back(gap);
                        //pe_cache[Pcontig_1][Econtig_1].push_back(gap);
                    }
                }
            }
            delete in ;
            loger<<BGIQD::LOG::lstart() << " total_pair_pe_diff "<<total_pair_pe_diff<<BGIQD::LOG::lend();
        }
    }

    void BuildPEGraph()
    {

        for( const auto & pair : pe_gap_cache )
        {
            unsigned int PcontigId = pair.first ;
            for( const auto & pair1 : pair.second )
            {
                std::vector<int> pos;
                if( pair1.second.size() < 10 ) 
                    continue ;
                std::vector<std::tuple<float,int,bool>> p2e;
                unsigned int EcontigId = pair1.first ;
                // P -> E
                for( int i = - 100 ; i < 900 ; i ++ )
                {
                    int incrR2 = contigLen_cache[PcontigId] + i ;
                    Dist tmp ;
                    tmp.Init(dist_bin,0,1000);
                    for( const auto & ape : pair1.second) 
                    {
                        pos.push_back(ape[0]);
                        pos.push_back(ape[1]);
                        pos.push_back(ape[2]+incrR2);
                        pos.push_back(ape[3]+incrR2);
                        std::sort(pos.begin() , pos.end());
                        int IS = pos[3]-pos[0] ;
                        tmp.Count(IS);
                    }
                    auto pet = tmp.CalcPercent();
                    auto keys = pet.ValidKeys();
                    auto base = percert_all.GetSubPercent(keys);
                    float sd ;
                    if( base.SD( pet , sd ) < 0.20 )
                    {
                        p2e.push_back(std::make_tuple(sd,i,true));
                    }
                }

                // R2 R1
                for( int i = - 100 ; i < 900 ; i ++ )
                {
                    int incrR1 = contigLen_cache[PcontigId] + i ;
                    Dist tmp ;
                    tmp.Init(dist_bin,0,1000);
                    for( const auto & ape : pair1.second) 
                    {
                        pos.push_back(ape[2]);
                        pos.push_back(ape[3]);
                        pos.push_back(ape[0]+incrR1);
                        pos.push_back(ape[1]+incrR1);
                        std::sort(pos.begin() , pos.end());
                        int IS = pos[3]-pos[0] ;
                        tmp.Count(IS);
                    }
                    auto pet = tmp.CalcPercent();
                    auto keys = pet.ValidKeys();
                    auto base = percert_all.GetSubPercent(keys);
                    float sd ;
                    if( base.SD( pet , sd ) > 0.75 )
                    {
                        p2e.push_back(std::make_tuple(sd,i,false));
                    }
                }
                std::sort(p2e.rbegin() , p2e.rend() );
                float sd ; int gap ; bool re;
                std::tie(sd,gap,re) = p2e[0] ;
                if( re )
                    std::cout<< PcontigId << "\t"<<EcontigId<<"\t"<<gap<<'\n';
                else
                    std::cout<< EcontigId << "\t"<<PcontigId<<"\t"<<gap<<'\n';
            }
        }

        for( const auto & pair : pe_cache )
        {
            unsigned int PcontigId = pair.first ;
            for( const auto & pair1 : pair.second )
            {
                unsigned int EcontigId = pair1.first ;
                if(ptest)
                {
                    std::string pe = std::to_string(PcontigId) + "\t"+std::to_string(EcontigId);
                    std::cout<<pe;
                }
                int total_len = 0 ;
                for( auto left : pair1.second )
                {
                    total_len += left ;
                    if(ptest)
                    {
                        std::cout<<'\t'<<left;
                    }
                }
                if(ptest)
                {
                    std::cout<<'\n';
                }
                int average_len = total_len / (int) pair1.second.size() ;
                pe_graph.AddEdge(PcontigId ,EcontigId , average_len , pair1.second.size());
            }
        }
    }

    void PrintContigPEGraph()
    {
        BGIQD::LOG::timer t(loger,"PrintContigPEGraph");
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.pe_graph());
        if( out == NULL )
            FATAL( " failed to open xxx.pe_graph for write ");
        pe_graph.PrintAsDOT(*out);
        delete out;
    }

}config;



int main(int argc , char ** argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string, prefix ,"prefix of files.");
    DEFINE_ARG_OPTIONAL(bool, test_data ,"print test data","no");
    DEFINE_ARG_OPTIONAL(int, insert_size ,"insert_size ","500");
    DEFINE_ARG_OPTIONAL(int, max_is ,"max valid insert_size","1000");
    DEFINE_ARG_OPTIONAL(int, dist_bin ,"bin size of insert_size distribution" ,"100");
    END_PARSE_ARGS;

    config.Init(prefix.to_string());
    config.ptest = test_data.to_bool();
    config.insert_size = insert_size.to_int();
    config.max_is= max_is.to_int();
    config.dist_bin = dist_bin.to_int();
    BGIQD::LOG::timer t(config.loger,"PEGraph");
    config.LoadSeeds();
    config.LoadPECahce();
    config.BuildPEGraph();
    config.PrintContigPEGraph();
    return 0 ;
}
