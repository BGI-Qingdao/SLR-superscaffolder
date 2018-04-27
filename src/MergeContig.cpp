#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/multithread/MultiThread.h"
#include "common/files/file_writer.h"
#include "common/freq/freq.h"

#include "biocommon/seq/tool_func.h"

#include "soap2/contigGraph.h"
#include "soap2/contigType.h"
#include "soap2/fileName.h"
#include "soap2/contigFasta.h"

#include "stLFR/LineGroup.h"

#include <iostream>
#include <string.h>
/**************************************************
 *
 * Tools to merge super contig into original data.
 *
 * All contig infomations for SOAPdenovos are :
 *      xxx.contig 
 *      xxx.Arc
 *      xxx.updated.edge
 *
 * This tool try to merge supercontig into original
 * info and re-generate above files.
 * Since this tool is design to called recursion ,
 * the result file will be in below format :
 *      xxx.contig_round_n
 *      xxx.Arc_round_n
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
    BGIQD::stLFR::ContigRoadFills linear_fill;
    BGIQD::SOAP2::ContigFastAMap contig_fasta_map;

    std::map< unsigned int , unsigned int > id_map;

    int K;
    int round ;

    enum ContigStatus 
    {
        UNKNOW = 0 ,
        ORIGINAL = 1 ,
        SUPER = 2 ,
        LINEAR = 3 ,
    };

    std::pair< ContigStatus , BGIQD::SOAP2::Edge & > GetEdge( unsigned int id )
    {
        if( id <= graph_ea.contigTotalNum )
        {
            return std::make_pair( ORIGINAL , std::ref(graph_ea.edge_array[id]) );
        }
        else if ( id <= graph_ea.contigTotalNum + new_graph_ea.contigTotalNum ) 
        {
            return std::make_pair( SUPER, std::ref(new_graph_ea.edge_array[id-graph_ea.contigTotalNum - 1]) );
        }
        else if ( id <= graph_ea.contigTotalNum + new_graph_ea.contigTotalNum + linear_graph_ea.contigTotalNum )
        {
            return std::make_pair( LINEAR, std::ref(linear_graph_ea.edge_array[id-graph_ea.contigTotalNum - new_graph_ea.contigTotalNum - 1]) );
        }
        static BGIQD::SOAP2::Edge error ;
        assert(0);
        return std::make_pair( UNKNOW , std::ref(error));
    }

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
        for( unsigned int i = 0 ; i < graph_ea.contigTotalNum; i++ )
        {
            detector.AddContigInfo(graph_ea.edge_array[i].length , graph_ea.edge_array[i].cov);
        }
        detector.GlobalAnalysis();
        // load super contig 
        fills.LoadContigRoadFills(fnames.contigroadfill(round));
        // load contigs fasta
        contig_fasta_map.Init(K) ;
        contig_fasta_map.LoadContig(fnames.contig(round));
        loger<<BGIQD::LOG::lstart()<<"Load contig "<<graph_ea.contigTotalNum<<" with base "<<contig_fasta_map.contigs.size()<<BGIQD::LOG::lend();
        contig_fasta_map.buildCompeleReverse();
        // init extra graph
        new_graph_ea.contigTotalNum = 0;
        linear_graph_ea.contigTotalNum = 0;
    }

    void DeleteUniqueContigInSuperContig()
    {
        AllocNewGraph_super() ;

        unsigned int new_contig_id = graph_ea.contigTotalNum ;
        unsigned int new_contig_id_new = 0 ;
        unsigned int new_arc_id = 0 ;
        int del_count = 0;

        new_graph_ea.contigTotalNum = fills.fills.size() * 2 ;

        for ( const auto & fill: fills.fills )
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
            del_count += 2 ;

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
                auto ret_to = GetEdge(id) ;
                auto ret_bal_to = GetEdge(ret_to.second.bal_id);
                //auto & bal = graph_ea.edge_array[graph_ea.edge_array[id].bal_id];
                auto & bal = ret_bal_to.second ;
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
                auto ret_to = GetEdge(id) ;
                auto ret_bal_to = GetEdge(ret_to.second.bal_id);
                //auto & bal = graph_ea.edge_array[graph_ea.edge_array[id].bal_id];
                auto & bal = ret_bal_to.second ;
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
                    del_count ++ ;
                    if( curr.bal_id != curr.id )
                    {
                        graph_ea.edge_array[curr.bal_id].SetDelete() ;
                    }
                }
            }
        }
        // [ 0 - contigTotalNum )
        assert( new_graph_ea.contigTotalNum == new_contig_id_new) ;
        // [ 0 - arcNum ]
        new_graph_ea.arcNum= new_arc_id -1 ;
        loger<<BGIQD::LOG::lstart()<<"DeleteUniqueContigInSuperContig merge "<<fills.fills.size()<<" super contig and delete "<<del_count<<" old contig "<<BGIQD::LOG::lend();
    }

    void GenerateNewContigSeq()
    {
        GenerateNewContigSeq_super();
        GenerateNewContigSeq_linear();
    }

    void ResetId()
    {
        // Reset original
        unsigned int new_id = 1 ;
        for( unsigned int i = 1 ; i <= graph_ea.contigTotalNum ; i++ )
        {
            const auto & curr = graph_ea.edge_array[i] ;
            if ( curr.IsDelete() )
                continue;
            id_map[curr.id] = new_id ++ ;
        }

        // Reset super
        for ( unsigned int i = 0 ; i < new_graph_ea.contigTotalNum ; i++ )
        {
            const auto & curr = new_graph_ea.edge_array[i] ;
            if ( curr.IsDelete() )
                continue;

            id_map[curr.id] = new_id ++ ;
            if( contig_fasta_map.contigs.at(curr.id).IsParlindorme() )
            {
                // merge 2 to 1 
                id_map[curr.id+1]=new_id;
                i ++ ;
            }
            else
            {
                id_map[curr.id+1]=new_id ++ ;
                i ++ ;
            }
        }

        // Reset linear
        for( unsigned int i = 0 ; i < linear_graph_ea.contigTotalNum ; i++ )
        {
            const auto & curr = linear_graph_ea.edge_array[i] ;
            if ( curr.IsDelete() )
                continue;

            id_map[curr.id] = new_id ++ ;
            if( contig_fasta_map.contigs.at(curr.id).IsParlindorme() )
            {
                // merge 2 to 1 
                id_map[curr.id+1]=new_id;
                i ++ ;
            }
            else
            {
                id_map[curr.id+1]=new_id ++ ;
                i ++ ;
            }
        }
    }

    void Linear()
    {
        // detect linear node  func
        auto detect_node_linear = [&]( BGIQD::SOAP2::Edge  & node )
        {
            if ( node.bal_id == node.id )
                return false ;
            BGIQD::SOAP2::Arc * arc ;
            // down
            arc = node.arc ;
            int count_to = 0 , count_from = 0;
            while( arc != NULL )
            {
                auto ret   = GetEdge(arc->to) ;
                arc = arc->next;
                if ( ret.first !=  ContigStatus::UNKNOW && ret.second.IsDelete() )
                {
                    continue;
                }
                count_to ++ ;
            }
            // up
            arc = GetEdge(node.bal_id).second.arc ;
            while( arc != NULL )
            {
                auto ret  = GetEdge(arc->to) ;
                arc = arc->next;
                if ( ret.first !=  ContigStatus::UNKNOW &&  ret.second.IsDelete() )
                {
                    continue;
                }
                count_from ++ ;
            }
            return ( count_to == 1 && count_from == 1 );
        };

        // detect linear  for all node
        for ( unsigned int i = 1 ; i <= graph_ea.contigTotalNum +new_graph_ea.contigTotalNum; i++ )
        {
            auto & curr = GetEdge(i).second;
            if( curr.bal_id != curr.id ) 
            {
                i ++ ;
            }
            else
            {
                curr.SetPalindrome();
                continue ;
            }
            if( curr.IsDelete() )
                continue ;
            if( detect_node_linear( curr ) )
            {
                curr.SetLinear();
                auto & bal= GetEdge(curr.bal_id).second;
                bal.SetLinear() ;
            }
        }
        linear_del = 0 ;
        // linear contig
        for( unsigned int i = 1 ; i <= graph_ea.contigTotalNum +new_graph_ea.contigTotalNum ; i++ )
        {
            LinearANode(i) ;
        }

        loger<<BGIQD::LOG::lstart()<<"Linear merge "<<linear_fill.fills.size()<<"linear contig and delete "<<linear_del<<" old contig "<<BGIQD::LOG::lend();
        //loger<<BGIQD::LOG::lstart()<<"linear del freq"<<"\n"<<l_del.ToString()<< BGIQD::LOG::lend();
        //loger<<BGIQD::LOG::lstart()<<"linear result freq"<<"\n"<<l_now.ToString()<< BGIQD::LOG::lend();
        // alloc graph
        AllocNewGraph_linear();
        // rebuild graph
        DeleteLinearContig();
    }


    void ReGenerate()
    {
        ReGenerate_Arc();
        ReGenerate_updatedEdge();
        ReGenerate_contig();
    }
    void PrintSuperOnly()
    {

    }
    void PrintSuperUsedOnly()
    {

    }
    private:
    int linear_del ;
    //BGIQD::FREQ::Freq<int> l_del;
    //BGIQD::FREQ::Freq<int> l_now;
    void LinearANode(unsigned int i)
    {
        auto & curr = GetEdge(i).second;
        if( curr.IsDelete() || curr.IsPalindrome() || ! curr.IsLinear() || curr.IsMarked() )
        {
            return ;
        }
        // get linear path
        auto & bal = GetEdge(curr.bal_id).second;
        std::vector<unsigned int> path_to ;
        std::vector<unsigned int> path_from ;
        GetLinearFromNode( curr , path_to ) ;
        GetLinearFromNode( bal , path_from );
        // mark used node
        int total_len = 0;
        if( path_to.size() ==  0 || path_from.size() == 0 )
        {
            return ;
        }
        else
        {
            linear_del ++ ;
            curr.SetMarked();
            bal.SetMarked() ;
            //total_len += curr.length ;
           // l_del.Touch(curr.length);
        }
        for( auto next : path_to )
        {
            linear_del ++ ;
            auto next_1 = GetEdge( next );
            auto next_2 = GetEdge( next_1.second.bal_id );
            next_1.second.SetMarked() ;
            next_2.second.SetMarked() ;
            //total_len += next_1.second.length ;
            //l_del.Touch(next_1.second.length);
        }
        std::stack<unsigned int> path_from_bal ;
        for( auto next : path_from )
        {
            linear_del ++ ;
            auto next_1 = GetEdge( next );
            auto next_2 = GetEdge( next_1.second.bal_id );
            path_from_bal.push(next_1.second.bal_id);
            next_1.second.SetMarked() ;
            next_2.second.SetMarked() ;
            //total_len += next_1.second.length ;
            //l_del.Touch(next_1.second.length);
        }
        std::vector<unsigned int> total ;
        while( ! path_from_bal.empty() )
        {
            total.push_back(path_from_bal.top());
            path_from_bal.pop();
        }
        total.push_back(curr.id ) ;
        for( auto next : path_to )
        {
            total.push_back(next) ; 
        }
        //l_now.Touch(total_len-K*(total.size()-1));
        linear_fill.fills.push_back(total) ;
    }

    void DeleteLinearContig() 
    {
        linear_graph_ea.contigTotalNum = linear_fill.fills.size() * 2 ;

        unsigned int new_contig_id = graph_ea.contigTotalNum + new_graph_ea.contigTotalNum ;
        unsigned int new_contig_id_new= 0 ;
        unsigned int new_arc_id = 0 ;

        for ( const auto & fill: linear_fill.fills )
        {
            auto & head = GetEdge(fill[0]).second ; 
            auto & head_bal = GetEdge(head.bal_id).second;
            auto & tail = GetEdge(*fill.rbegin()).second;
            auto & tail_bal = GetEdge(tail.bal_id).second;

            //delete old head & tail
            head.SetDelete();
            head_bal.SetDelete();
            tail.SetDelete();
            tail_bal.SetDelete();

            //make new contig
            new_contig_id ++ ;
            auto & new_contig = linear_graph_ea.edge_array[new_contig_id_new++];
            new_contig.id = new_contig_id  ;
            new_contig.bal_id = new_contig_id +1   ;
            new_contig.arc = tail.arc ;

            new_contig_id ++ ;
            auto & new_contig_bal = linear_graph_ea.edge_array[new_contig_id_new++]; 
            new_contig_bal.id = new_contig_id ;
            new_contig_bal.bal_id = new_contig_id - 1 ;
            new_contig_bal.arc = head_bal.arc ;

            // let old data point to new contig
            volatile BGIQD::SOAP2::Arc * to_me = head_bal.arc ;
            while( to_me != NULL )
            {
                unsigned int id = to_me->to ;
                auto ret_to = GetEdge(id) ;
                auto ret_bal_to = GetEdge(ret_to.second.bal_id);
                //auto & bal = graph_ea.edge_array[graph_ea.edge_array[id].bal_id];
                auto & bal = ret_bal_to.second ;
                auto & curr_arc = linear_graph_ea.arc_array[new_arc_id++];
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
                auto ret_to = GetEdge(id) ;
                auto ret_bal_to = GetEdge(ret_to.second.bal_id);
                //auto & bal = graph_ea.edge_array[graph_ea.edge_array[id].bal_id];
                auto & bal = ret_bal_to.second ;
                auto & curr_arc = linear_graph_ea.arc_array[new_arc_id++] ;
                curr_arc.to = new_contig_bal.id ;
                curr_arc.cov = me_to->cov ;
                curr_arc.next = bal.arc;
                bal.arc = &curr_arc ;

                me_to = me_to->next ;
            }

            // delete old contig .
            for ( size_t i = 1 ; i < fill.size() -1 ; i++ )
            {
                auto & curr = GetEdge(fill[i]).second;
                curr.SetDelete();
                if( curr.bal_id != curr.id )
                {
                    GetEdge(curr.bal_id).second.SetDelete() ;
                }
            }
        }
        // [ 0 - contigTotalNum )
        assert(linear_graph_ea.contigTotalNum == new_contig_id_new) ;
        // [ 0 - arcNum ]
        linear_graph_ea.arcNum= new_arc_id -1 ;
    }

    void GetLinearFromNode(const BGIQD::SOAP2::Edge & edge , std::vector<unsigned int> & path)
    {
        auto get_valid_arc = [&](const BGIQD::SOAP2::Edge & edge)
        {
            BGIQD::SOAP2::Arc* arc = edge.arc ;
            while( arc != NULL )
            {
                auto ret = GetEdge(arc->to) ;
                if ( ret.first !=  ContigStatus::UNKNOW && ! ret.second.IsDelete() )
                {
                    return arc ;
                }
                arc = arc->next;
            }
            return (BGIQD::SOAP2::Arc* )NULL ;
        };
        BGIQD::SOAP2::Arc * ret = get_valid_arc(edge) ;
        BGIQD::SOAP2::Edge * next ;
        std::set<unsigned int> uniques ;
        uniques.insert(edge.id);
        next = &GetEdge( ret->to ).second ;
        while( next && next->IsLinear() && uniques.find( next->id) == uniques.end() )
        {
            uniques.insert( next->id );
            path.push_back( ret->to ) ;
            ret = get_valid_arc(*next);
            next = &GetEdge(ret->to).second;
        }
        return ;
    }

    void GenerateNewContigSeq_super()
    {
        for( unsigned int i = 0 ; i < new_graph_ea.contigTotalNum ; i=i+2 )
        {
            auto & curr = new_graph_ea.edge_array[i] ;
            const auto & line = fills.fills[i/2] ;
            auto ret = contig_fasta_map.MergeContig(line);
            assert( ret.id == curr.id );
            curr.length = ret.length + K ;
            unsigned int head = line[0];
            unsigned int tail = line[line.size()-1];
            memcpy ( curr.from , GetEdge(head).second.from , sizeof(BGIQD::SOAP2::Kmer));
            memcpy ( curr.to, GetEdge(tail).second.to, sizeof(BGIQD::SOAP2::Kmer));
            curr.cov = ret.cov * 10 ;
            ret.MarkBase() ;
            if( BGIQD::SEQ::isSeqPalindrome( ret.K + ret.linear))
                ret.MarkParlindorme();
            else
            {
                contig_fasta_map.contigs[ret.id + 1] = ret.ReverseCompelete() ;
            }
            contig_fasta_map.contigs[ret.id] = ret;
            // bal_id
            auto & bal= new_graph_ea.edge_array[i+1] ;
            bal.length = curr.length ;
            bal.cov = curr.cov;
            unsigned int head_bal = GetEdge(tail).second.bal_id;
            unsigned int tail_bal = GetEdge(head).second.bal_id;
            memcpy ( bal.from , GetEdge(head_bal).second.from , sizeof(BGIQD::SOAP2::Kmer));
            memcpy ( bal.to, GetEdge(tail_bal).second.to, sizeof(BGIQD::SOAP2::Kmer));
        }
    }

    void GenerateNewContigSeq_linear()
    {
        for( unsigned int i = 0 ; i < linear_graph_ea.contigTotalNum ; i=i+2 )
        {
            auto & curr = linear_graph_ea.edge_array[i] ;
            const auto & line = linear_fill.fills[i/2] ;
            auto ret = contig_fasta_map.MergeContig(line);
            assert( ret.id == curr.id );
            curr.length = ret.length + K ;
            unsigned int head = line[0];
            unsigned int tail = line[line.size()-1];
            memcpy ( curr.from , GetEdge(head).second.from , sizeof(BGIQD::SOAP2::Kmer));
            memcpy ( curr.to, GetEdge(tail).second.to, sizeof(BGIQD::SOAP2::Kmer));
            curr.cov = ret.cov * 10 ;
            ret.MarkBase() ;
            if( BGIQD::SEQ::isSeqPalindrome( ret.K + ret.linear))
                ret.MarkParlindorme();
            contig_fasta_map.contigs[ret.id] = ret;
            // bal_id
            auto & bal= linear_graph_ea.edge_array[i+1] ;
            bal.length = curr.length ;
            bal.cov = curr.cov;
            unsigned int head_bal = GetEdge(tail).second.bal_id;
            unsigned int tail_bal = GetEdge(head).second.bal_id;
            memcpy ( bal.from , GetEdge(head_bal).second.from , sizeof(BGIQD::SOAP2::Kmer));
            memcpy ( bal.to, GetEdge(tail_bal).second.to, sizeof(BGIQD::SOAP2::Kmer));
        }
    }

    unsigned int NewId(unsigned int old)
    {
        return id_map.at(old);
    }

    void AllocNewGraph_super()
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

    void AllocNewGraph_linear()
    {
        int new_arc_num = 0 ;
        for( const auto & fill : linear_fill.fills)
        {
            const auto & head = GetEdge(fill[0]).second ; 
            const auto & head_bal = GetEdge(head.bal_id).second;
            const auto & tail = GetEdge(*fill.rbegin()).second;
            const auto & tail_bal = GetEdge(tail.bal_id).second;
            new_arc_num += head.ArcNum() + head_bal.ArcNum() + tail.ArcNum() + tail_bal.ArcNum() ;
        }
        int new_contig_num = linear_fill.fills.size() * 2 ;

        linear_graph_ea.edge_array =( BGIQD::SOAP2::Edge* ) calloc ( sizeof(BGIQD::SOAP2::Edge), new_contig_num );
        linear_graph_ea.arc_array = (BGIQD::SOAP2::Arc *) calloc ( sizeof(BGIQD::SOAP2::Arc) , new_arc_num );

        linear_graph_ea.contigTotalNum = 0 ;
        linear_graph_ea.arcNum = 0 ;
    }

    void ReGenerate_Arc()
    {
        auto print = [&](std::ostream & ost , const BGIQD::SOAP2::Edge & curr)
        {
            if ( curr.IsDelete() )
                return ;
            ost<<NewId(curr.id)<<'\t';
            volatile BGIQD::SOAP2::Arc * next = curr.arc ;
            while(next)
            {
                if(! GetEdge(next->to).second.IsDelete() )
                { 
                    ost<<NewId(next->to)<<"\t"<<(int)next->cov<<"\t";
                }
                next = next->next ;
            }
            ost<<std::endl;
        };
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fnames.Arc(round+1));
        for( unsigned int i = 1 ; i <= graph_ea.contigTotalNum ; i++ )
        {
            const auto & curr = graph_ea.edge_array[i] ;
            print(*out,curr);
        }
        bool base = true ;
        for( unsigned int i = 0 ; i < new_graph_ea.contigTotalNum ; i++ )
        {
            const auto & curr = new_graph_ea.edge_array[i] ;
            print( *out , curr);
            if( base )
            {
                // pass the bal_id that not exsist !
                if( contig_fasta_map.contigs.at(curr.id).IsParlindorme() )
                    i ++ ;
                else
                    base = false ;
            }
            else
            {
                base = true ;
            }
        }
        for( unsigned int i = 0 ; i < linear_graph_ea.contigTotalNum ; i++ )
        {
            const auto & curr = linear_graph_ea.edge_array[i] ;
            print( *out , curr);
            if( base )
            {
                // pass the bal_id that not exsist !
                if( contig_fasta_map.contigs.at(curr.id).IsParlindorme() )
                    i ++ ;
                else
                    base = false ;
            }
            else
            {
                base = true ;
            }
        }
        delete out;
    }

    void ReGenerate_updatedEdge()
    {
        auto print = []( std::ostream & out, const BGIQD::SOAP2::Edge & curr )
        {
            (out)<<">length "<<curr.length
                <<','<<(int)(curr.bal_id-curr.id)
                <<','<<(int)(curr.cov)
                <<','<<std::hex<<curr.from[0]
                <<" "<<std::hex<<curr.from[1]
#if K127mer
                <<" "<<std::hex<<curr.from[2]
                <<" "<<std::hex<<curr.from[3]
#endif
                <<','<<std::hex<<curr.to[0]
                <<" "<<std::hex<<curr.to[1]
#if K127mer
                <<" "<<std::hex<<curr.to[2]
                <<" "<<std::hex<<curr.to[3]
#endif
                <<std::dec<<std::endl;
        };
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fnames.updatedEdge(round+1));
        unsigned int line_num = 1 ;

        for( unsigned int i = 1 ; i <= graph_ea.contigTotalNum ; i++ )
        {
            const auto & curr = graph_ea.edge_array[i] ;
            if ( curr.IsDelete() )
                continue;
            assert( line_num = NewId(curr.id)) ;
            print(*out , curr );
            line_num++;
        }
        bool base = true ;
        for( unsigned int i = 0 ; i < new_graph_ea.contigTotalNum ; i++ )
        {
            const auto & curr = new_graph_ea.edge_array[i] ;
            if ( curr.IsDelete() )
                continue;
            assert( line_num = NewId(curr.id)) ;
            print( *out, curr);
            // pass the bal_id that not exsist !
            if( base )
            {
                // pass the bal_id that not exsist !
                if( contig_fasta_map.contigs.at(curr.id).IsParlindorme() )
                    i ++ ;
                else
                    base = false ;
            }
            else
            {
                base = true ;
            }
        }
        for( unsigned int i = 0 ; i < linear_graph_ea.contigTotalNum ; i++ )
        {
            const auto & curr = linear_graph_ea.edge_array[i] ;
            if ( curr.IsDelete() )
                continue;
            assert( line_num = NewId(curr.id)) ;
            print( *out, curr);
            // pass the bal_id that not exsist !
            if( base )
            {
                // pass the bal_id that not exsist !
                if( contig_fasta_map.contigs.at(curr.id).IsParlindorme() )
                    i ++ ;
                else
                    base = false ;
            }
            else
            {
                base = true ;
            }
        }
        delete out;
    }

    void ReGenerate_contig( )
    {
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fnames.contig(round+1));
        int final_num = 0 ;
        int del_count = 0 ;
        for( unsigned int i = 1 ; i <= graph_ea.contigTotalNum ; i++ )
        {
            const auto & curr = graph_ea.edge_array[i] ;
            // check base
            if( curr.id > curr.bal_id ) 
                continue ;
            // check valid
            if ( curr.IsDelete() || curr.length < 1 )
            {
                del_count ++ ;
                continue;
            }
            const auto & c = contig_fasta_map.contigs.at(curr.id) ;
            (*out)<<c.ToString(NewId(curr.id),false)<<std::endl;
            final_num ++ ;
        }
        int step = final_num;
        loger<<BGIQD::LOG::lstart()<<" print base contig "<<final_num;
        loger<<" detect delete "<<del_count<<BGIQD::LOG::lend();

        for( unsigned int i = 0 ; i < new_graph_ea.contigTotalNum ; i+=2 )
        {
            const auto & curr = new_graph_ea.edge_array[i] ;
            if ( curr.IsDelete() || curr.length < 1)
                continue;
            const auto & c = contig_fasta_map.contigs.at(curr.id) ;
            final_num ++ ;
            if( fills.has_circle[i/2] )
                (*out)<<c.ToString(NewId(curr.id),true)<<std::endl;
            else
                (*out)<<c.ToString(NewId(curr.id),false)<<std::endl;

        }
        loger<<BGIQD::LOG::lstart()<<" print super contig "<<final_num-step<<BGIQD::LOG::lend();
        step = final_num;

        for( unsigned int i = 0 ; i < linear_graph_ea.contigTotalNum ; i+=2 )
        {
            const auto & curr = linear_graph_ea.edge_array[i] ;
            if ( curr.IsDelete() || curr.length < 1)
                continue;
            const auto & c = contig_fasta_map.contigs.at(curr.id) ;
            final_num ++ ;
            (*out)<<c.ToString(NewId(curr.id),false)<<std::endl;
        }
        loger<<BGIQD::LOG::lstart()<<" print linear contig "<<final_num-step<<BGIQD::LOG::lend();
        loger<<BGIQD::LOG::lstart()<<" print final contig "<<final_num<<BGIQD::LOG::lend();
        loger<<BGIQD::LOG::lstart()<<" final fasta "<<contig_fasta_map.contigs.size()<<BGIQD::LOG::lend();

        delete out;
    }

} config;

int main(int argc , char **argv)
{

    BGIQD::LOG::timer t(config.loger,"MergeContig");

    START_PARSE_ARGS
    DEFINE_ARG_DETAIL(std::string , prefix, 'o',false,"prefix");
    DEFINE_ARG_DETAIL(int , kvalue, 'K',false,"K value");
    DEFINE_ARG_DETAIL(int , t_num, 't',true,"thread num . default[8]");
    DEFINE_ARG_DETAIL(int , round, 'r',true,"round num. default[0]");
    DEFINE_ARG_DETAIL(bool ,linear, 'l',true,"call linear ?. default false");
    DEFINE_ARG_DETAIL(bool ,super_only, 's' , true , "only print super contig ? default false");
    DEFINE_ARG_DETAIL(bool ,used_only, 's' , true , "only print super used  contig ? default false");
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
    // init
    config.Init(prefix.to_string() , round.to_int() , kvalue.to_int());
    // graph simplify step 1
    config.DeleteUniqueContigInSuperContig();
    // graph simplify step 2
    if( linear.to_bool() )
    {
        config.Linear();
    }
    // management  date
    config.GenerateNewContigSeq();
    // management  date order
    config.ResetId();
    // print
    config.ReGenerate();

    return 0 ;
}
