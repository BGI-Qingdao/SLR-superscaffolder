
#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/multithread/MultiThread.h"
#include "common/files/file_writer.h"

#include "soap2/contigGraph.h"
#include "soap2/contigType.h"
#include "soap2/fileName.h"
#include "soap2/contigFasta.h"

#include "stLFR/LineGroup.h"

#include <iostream>

/**************************************************
 *
 * Tools to merge super contig into original data.
 *
 * All contig infomations for SOAPdenovos are :
 *      xxx.contig 
 *      xxx.Arc
 *      xxx.contigIndex
 *      xxx.updated.edge
 *
 * This tool try to merge supercontig into original
 * info and re-generate above files.
 * Since this tool is design to called recursion ,
 * the result file will be in below format :
 *      xxx.contig_round_n
 *      xxx.Arc_round_n
 *      xxx.contigIndex_round_n
 *      xxx.updated.edge_round_n
 *
 * ***********************************************/


struct AppConfig
{
    BGIQD::SOAP2::FileNames fnames;
    BGIQD::LOG::logger loger;
    BGIQD::SOAP2::ContigTypeDetecter detector;
    BGIQD::SOAP2::GraphEA graph_ea; // original graph
    BGIQD::SOAP2::GraphEA new_graph_ea; // for new graph from super_contig
    BGIQD::SOAP2::GraphEA linear_graph_ea; // for new graph from linear step

    BGIQD::stLFR::ContigRoadFills fills;
    BGIQD::SOAP2::ContigFastAMap contig_fasta_map;

    int K;
    int round ;

    void Init(const std::string & prefix , int rd, int k )
    {
        K = k ;
        round = rd;
        // init logger
        BGIQD::LOG::logfilter::singleton().get("MergeContig",BGIQD::LOG::loglevel::INFO , loger);
        // init files name
        fnames.Init(prefix);
        // load original data 
        graph_ea.LoadEdge(fnames.updatedEdge(round),K);
        graph_ea.LoadArc(fnames.Arc(round));
        // prepare for detect contig type 
        detector.Init(K);
        for(unsigned int i = 0 ; i < graph_ea.contigTotalNum; i++ )
        {
            detector.AddContigInfo(graph_ea.edge_array[i].length , graph_ea.edge_array[i].cov);
        }
        detector.GlobalAnalysis();
        // load super contig 
        fills.LoadContigRoadFills(fnames.contigroadfill(round));
        // load contigs fasta
        contig_fasta_map.Init(K) ;
        contig_fasta_map.LoadContig(fnames.contig(round));
    }

    void DeleteUniqueContigInSuperContig()
    {
        AllocNewGraph() ;

        unsigned int new_contig_id = graph_ea.contigTotalNum ;
        unsigned int new_contig_id_new = 0 ;
        unsigned int new_arc_id = 0 ;
        for( const auto & fill: fills.fills )
        {
            auto & head = graph_ea.edge_array[fill[0]];
            auto & head_bal = graph_ea.edge_array[head.bal_id];
            auto & tail = graph_ea.edge_array[*fill.rbegin()];
            auto & tail_bal = graph_ea.edge_array[tail.bal_id];

            //delete old head & tail
            head.SetDelete();
            head_bal.SetDelete();
            tail.SetDelete();
            tail_bal.SetDelete();

            //make new contig
            new_contig_id ++ ;
            auto & new_contig = new_graph_ea.edge_array[new_contig_id_new++];
            new_contig.id = new_contig_id  ;
            new_contig.bal_id = new_contig_id +1   ;
            new_contig.arc = tail.arc ;

            new_contig_id ++ ;
            auto & new_contig_bal = new_graph_ea.edge_array[new_contig_id_new++]; 
            new_contig_bal.id = new_contig_id ;
            new_contig_bal.bal_id = new_contig_id - 1 ;
            new_contig_bal.arc = head_bal.arc ;

            // let old data point to new contig
            volatile BGIQD::SOAP2::Arc * to_me = head_bal.arc ;
            while( to_me != NULL )
            {
                unsigned int id = to_me->to ;
                auto & bal = graph_ea.edge_array[graph_ea.edge_array[id].bal_id];
                auto & curr_arc = new_graph_ea.arc_array[new_arc_id++];
                curr_arc.to = new_contig.id ;
                curr_arc.cov = to_me->cov ;
                curr_arc.next = bal.arc;
                bal.arc = &curr_arc ;
                to_me = to_me->next ;
            }

            // let old data point to new contig
            volatile BGIQD::SOAP2::Arc * me_to = tail.arc ;
            while( me_to != NULL )
            {
                unsigned int id = me_to->to ;
                auto & bal = graph_ea.edge_array[graph_ea.edge_array[id].bal_id];
                auto & curr_arc = new_graph_ea.arc_array[new_arc_id++];
                curr_arc.to = new_contig_bal.id ;
                curr_arc.cov = me_to->cov ;
                curr_arc.next = bal.arc;
                bal.arc = &curr_arc ;

                me_to = me_to->next ;
            }

            // delete contig within fill .
            for ( size_t i = 1 ; i < fill.size() -1 ; i++ )
            {
                auto & curr = graph_ea.edge_array[fill[i]];
                if( detector.ContigType(curr.length , curr.cov ) == BGIQD::SOAP2::ContigTypeDetecter::Type::Unique )
                {
                    curr.SetDelete();
                    if( curr.bal_id != curr.id )
                    {
                        graph_ea.edge_array[curr.bal_id].SetDelete() ;
                    }
                }
            }
        }
        // [ 0 - contigTotalNum ]
        new_graph_ea.contigTotalNum = new_contig_id ;
        // [ 0 - arcNum ]
        new_graph_ea.arcNum= new_arc_id -1 ;
    }

    void Linear() 
    {

    }

    void ReGenerate()
    {

    }

    private:

        void AllocNewGraph()
        {
            int new_arc_num = 0 ;
            for( const auto & fill : fills.fills)
            {
                const auto & head = graph_ea.edge_array[fill[0]];
                const auto & head_bal = graph_ea.edge_array[head.bal_id];
                const auto & tail = graph_ea.edge_array[*fill.rbegin()];
                const auto & tail_bal = graph_ea.edge_array[tail.bal_id];
                new_arc_num += head.ArcNum() + head_bal.ArcNum() + tail.ArcNum() + tail_bal.ArcNum() ;
            }
            int new_contig_num = fills.fills.size() * 2 ;

            new_graph_ea.edge_array =( BGIQD::SOAP2::Edge* ) calloc ( sizeof(BGIQD::SOAP2::Edge), new_contig_num );
            new_graph_ea.arc_array = (BGIQD::SOAP2::Arc *) calloc ( sizeof(BGIQD::SOAP2::Arc) , new_arc_num );

            new_graph_ea.contigTotalNum = 0 ;
            new_graph_ea.arcNum = 0 ;
        }

        void ReGenerate_Arc()
        {
            auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fnames.Arc(round+1));
            for( unsigned int i = 1 ; i <= graph_ea.contigTotalNum ; i++ )
            {
                const auto & curr = graph_ea.edge_array[i] ;
                (*out)<<curr.id<<"\t";
                volatile BGIQD::SOAP2::Arc * next = curr.arc ;
                while(next)
                {
                    (*out)<<next->to<<"\t"<<(int)next->cov<<"\t";
                    next = next->next ;
                }
            }
            for( unsigned int i = 0 ; i <= new_graph_ea.contigTotalNum ; i++ )
            {
                const auto & curr = new_graph_ea.edge_array[i] ;
                (*out)<<curr.id<<"\t";
                volatile BGIQD::SOAP2::Arc * next = curr.arc ;
                while(next)
                {
                    (*out)<<next->to<<"\t"<<(int)next->cov<<"\t";
                    next = next->next ;
                }
            }
            delete out;
        }

        void ReGenerate_updatedEdge()
        {

        }
}config;


int main(int argc , char **argv)
{

    BGIQD::LOG::timer t(config.loger,"MergeContig");

    START_PARSE_ARGS
        DEFINE_ARG_DETAIL(std::string , prefix, 'o',false,"prefix");
        DEFINE_ARG_DETAIL(int , kvalue, 'K',false,"K value");
        DEFINE_ARG_DETAIL(int , t_num, 't',true,"thread num . default[8]");
        DEFINE_ARG_DETAIL(int , round, 'r',true,"round num. default[0]");
    END_PARSE_ARGS

    if(! t_num.setted )
    {
        t_num.setted = true;
        t_num.d.i = 8;
    }
    if( !round.setted )
    {
        round.setted = true ;
        round.d.i = 0;
    }

    config.Init(prefix.to_string() , round.to_int() , kvalue.to_int());
    config.DeleteUniqueContigInSuperContig();
    config.Linear();
    config.ReGenerate();
    return 0 ;
}
