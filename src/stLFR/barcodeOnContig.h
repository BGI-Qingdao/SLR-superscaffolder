#ifndef __STLFR_BARCODEONCONTIG_H__
#define __STLFR_BARCODEONCONTIG_H__

#include <set>
#include <map>
#include <vector>
#include "soap2/contigGraph.h"
#include "stLFR/LineGroup.h"
#include "stLFR/ContigCluster.h"
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
            // used below 2 intareface to construct sub_p2p_graph.
            // AddPath must be called only once !!
            void AddPath( unsigned int to ,const std::vector<std::list<SOAP2::Edge>> & paths );
            // AddMid can be called multi-times!!!
            void AddMid( unsigned int to ,const std::vector<std::list<SOAP2::Edge>> & paths );

            struct Path 
            {
                std::vector<unsigned int > paths;
                int total_length;
                float cov;
                float barcode_cov ;

                void Init()
                {
                    total_length = 0 ;
                    cov = 0;
                    total_barcode = 0;
                    total_cov = 0;
                    barcode_cov = 0;
                }

                void AddEdge( int id , int length , float cov , int barcode ) 
                {
                    paths.push_back(id);
                    total_length += length ;
                    total_cov += length * cov ;
                    total_barcode += barcode ;
                }

                void CalcCov()
                {
                    cov = total_cov / total_length ;
                    barcode_cov = total_barcode / total_length ;
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
            Path final_path;

            void GeneratePath();
            private:
            void findAllPath();
            void CleanAndSavePath();
            void InitEdge( unsigned int id);
            std::vector<Path> allPaths;
            void findAllPath(  unsigned int id  , Path  p);
            void ScoreAllPath();
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
