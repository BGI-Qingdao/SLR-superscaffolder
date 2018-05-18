#ifndef __STLFR_BARCODEONCONTIG_H__
#define __STLFR_BARCODEONCONTIG_H__

#include <set>
#include <map>
#include <vector>
#include "soap2/contigGraph.h"
#include "stLFR/LineGroup.h"
#include "stLFR/ContigCluster.h"
#include "common/flags/flags.h"
namespace BGIQD {
    namespace stLFR {

        typedef unsigned int ContigId; 
        typedef std::map<int, int> BarcodeOnContig;

        struct GraphEA_withBarcode
        {
            BGIQD::SOAP2::GraphEA graph_ea;
            std::map<ContigId,BarcodeOnContig > barcode_on_contig;

            static int Union(const std::map<int,int> & m1 
                    , const std::map<int,int> & m2 
                    , std::map<int, int> & um ); 

            void LoadBarcodeOnConfig(const std::string & file);
        };

        // --------------------SubP2PGraph --------------------

        //
        // This is a sub-graph of GraphEA
        // 
        // In this sub-graph
        //      * all vertex start from 1 root vertex.
        //      * all vertex end to 1 target vertex.
        // 
        // This sub-graph is used to describe the all paths
        //  between 2 key edge, if we can find the 1 correct
        //  path we can use it as a super contig.
        //
        struct  P2PGraph
        {

            GraphEA_withBarcode * base_graph;

            BarcodeOnContig root_target_union;

            unsigned int target;

            unsigned int root ;

            int final_circled ;

            bool deal_circle;

            int K ;

            struct Edge
            {
                unsigned int id;
                unsigned int barcode_cov;
                float cov;
                int length ;
                std::set<unsigned int> tos;
                //std::set<unsigned int> froms;
            };

            std::map<unsigned int , P2PGraph::Edge> sub_graph;

            void Init( unsigned int from , unsigned int to);

            void AddFromTo( unsigned int from , unsigned int to );

            bool CheckPalindrome() const ;

            struct Circle
            {
                std::vector<Edge> cpath;
                std::set<unsigned int>  csets;
                FLAGS_INT ;
                ADD_A_FLAG( 0, set );
                int  circle_run ;
                void SetCircle( const std::vector<Edge> & path , unsigned int root , float ecov)
                {
                    bool flag = false ;
                    for( size_t  i = 0 ; i < path.size() ; i++ )
                    {
                        if( ! flag && path[i].id != root )
                            continue ;
                        flag = true ;
                        cpath.push_back(path[i]);
                        csets.insert(path[i].id);
                    }
                    int total_len = 0 ;
                    float total_cov = 0.0f;
                    for( const auto & i : cpath )
                    {
                        total_len += i.length ;
                        total_cov += i.length * i.cov ; 
                    }
                    float me_cov = total_cov / total_len; 
                    float me_dup = me_cov / ecov ;
                    int me_more1 = (int)me_dup ;
                    int me_more = me_more1;
                    std::vector<Edge> path_base = cpath ;
                    while( me_more > 1  || ( me_more == 1 && me_dup - me_more1 > 0.5f ) )
                    {
                        cpath.insert(cpath.end(), path_base.begin() , path_base.end() );
                        me_more -- ;
                    }
                    circle_run = cpath.size() / csets.size() ;
                    if( ! csets.empty() )
                        Set_set();
                }

                void Clean()
                {
                    cpath.clear();
                    csets.clear();
                    Clean_set();
                    circle_run = 0 ;
                }
            };

            struct Path
            {
                //std::vector<unsigned int > paths;
                std::vector<Edge> paths;
                std::set<unsigned int > nodes;

                Circle circle ;
                int total_length;
                float cov;
                float barcode_cov ;
                int K ;
                void Init()
                {
                    total_length = 0 ;
                    cov = 0;
                    total_barcode = 0;
                    total_cov = 0;
                    barcode_cov = 0;
                    circle.Clean();
                }

                void AddCircle( const Circle & c)
                {
                    circle  =c ;
                }

                void MergeCircle();

                bool IsPathInCircle(const Circle & c)
                {
                    if( ! c.Is_set() )
                        return false ;
                    for(auto i: paths )
                    {
                        if( c.csets.find( i.id ) != c.csets.end() )
                        {
                            return true ;
                        }
                    }
                    return false;
                }


                //bool AddEdge( int id , int length , float cov , int barcode ) 
                bool AddEdge(const Edge & edge ) 
                {
                    if( nodes.find( edge.id ) == nodes.end() )
                    {
                        Edge tmp = edge ;
                        paths.push_back(tmp);
                        nodes.insert(edge.id);
                        return true;
                    }
                    else
                        return false; //circle detected.
                }

                void CalcCov()
                {
                    total_cov = 0 ; 
                    total_length = 0 ;
                    total_barcode = 0 ;
                    for( const auto & edge : paths )
                    {
                        total_length += edge.length - K  ;
                        total_cov += edge.length * cov ;
                        total_barcode += edge.barcode_cov ;
                    }
                    if( total_length != 0 )
                    {
                        cov = total_cov / total_length ;
                        barcode_cov = total_barcode / total_length ;
                    }
                    else
                    {
                        cov = 255 ;
                        barcode_cov = 100000;
                    }
                }

                bool operator < ( const Path & a) const
                {
                    if( barcode_cov != a.barcode_cov )
                        return barcode_cov < a.barcode_cov ;
                    if( cov != a.cov )
                        return cov < a.cov ;
                    return total_length < a.total_length ;
                }
                private:
                int total_barcode;
                float total_cov;
            };


            int path_num;
            bool IsOK();
            std::vector<unsigned int> final_path;
            float ecov ;
            void GeneratePath();
            int ShortestPath();
            private:
            void findAllPath();
            void CleanAndSavePath();
            void InitEdge( unsigned int id);
            std::vector<Path> allPaths;
            void findAllPath(  unsigned int id  , Path  p , Circle & circle_detected );
            void ScoreAllPath();
            void ScoreAllPathByLength() ;

            //struct SubP2PGraphEdge
            //{
            //    unsigned int id;
            //    int x;
            //};
            //

            //SubP2PGraphEdge root_graph;

            //// after addPath and addMid , call this to fill root_graph
            //void ConstructSubP2pGraph();
            //private:
            //    //void ConstructSubP2pGraph(SubP2PGraph & spg);
        };

//        struct MergedEdge
//        {
//            std::vector<unsigned int>  contigs;
//        };
//
    }// namespace stLFR
}// namespace BGIQD

#endif //__STLFR_BARCODEONCONTIG_H__
