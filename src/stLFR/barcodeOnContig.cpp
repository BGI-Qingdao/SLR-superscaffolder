#include "stLFR/barcodeOnContig.h"
#include <cassert>
#include <algorithm>
#include "common/files/file_reader.h"
#include <sstream>
namespace BGIQD {
    namespace stLFR {

        // ---------------- GraphEA_withBarcode----------------

        int GraphEA_withBarcode::Union(const std::map<int,int> & m1 
                    , const std::map<int,int> & m2 
                    , std::map<int, int> & um )
        {

            int ret = 0 ;
            for(const auto pair : m2)
            {
                if( m1.find( pair.first ) != m1.end() )
                {
                    int min = pair.second < m1.at(pair.first)
                        ? pair.second
                        : m1.at(pair.first) ;

                    ret += min ;
                    um[pair.first] = min ;
                }
                ret += pair.second;
            }
            return ret;
        }
        // ---------------- P2PGraph -------------------------
        void P2PGraph::Init( unsigned int from , unsigned int to)
        {
            target = to ;
            root= from;
            InitEdge(to);
            InitEdge(from);
            GraphEA_withBarcode::Union(base_graph->barcode_on_contig[from] ,base_graph->barcode_on_contig[to],root_target_union);
        }

        void P2PGraph::InitEdge( unsigned int id)
        {
            if(sub_graph.find(id) == sub_graph.end() )
            {
                auto & node = base_graph->graph_ea.edge_array[id];
                sub_graph[id].id = id;
                sub_graph[id].length =  node.id;
                sub_graph[id].cov = node.cov;
                if( id != root && id != target )
                {
                    BarcodeOnContig tmp;
                    if( base_graph->barcode_on_contig.find(id)  != base_graph->barcode_on_contig.end() )
                        sub_graph[id].barcode_cov = GraphEA_withBarcode::Union(base_graph->barcode_on_contig[id],root_target_union,tmp);
                }
            }
        }

        void P2PGraph::AddPath( unsigned int to, const std::vector<std::list<SOAP2::Edge> > & paths )
        {
            assert( target == to ) ;
            if ( paths.empty() )
                return ;
            assert ( root== (paths[0].rbegin()->id) ) ;

            for( const auto & l : paths )
            {
                unsigned int next = target ;
                for( auto i = l.begin() ; i != l.end() ; i++)
                {
                    if( i->id == to )
                        continue;
                    InitEdge( i->id );
                    sub_graph[i->id].tos.insert(next);
                    next  = i->id ;
                }
            }
        }

        void P2PGraph::AddMid( unsigned int to, const std::vector<std::list<SOAP2::Edge> > & paths )
        {
            //assert(root_graph.target == to) ;
            if ( paths.empty() )
                return ;
            assert(root== (paths[0].rbegin()->id)) ;
            // check if to is valid
            if( sub_graph.find( to ) == sub_graph.end() )
            {
                return ;
            }

            for( const auto & l : paths )
            {
                unsigned int next = to ;
                for( auto i = l.begin() ; i != l.end() ; i++ )
                {
                    if( i->id == to )
                        continue;
                    InitEdge( i->id );
                    sub_graph[i->id].tos.insert(next);
                    next  = i->id ;
                }
            }
        }

        void P2PGraph::findAllPath()
        {
            Path p;
            p.Init();
            path_num = 0 ;
            auto & node = sub_graph[root];
            for( const auto &i : node.tos)
            {
                findAllPath(i, p);
                if( path_num < 0 )
                    return ;
            }
            if( path_num < 0 )
            {
                allPaths.clear();
            }
            else
            {
                path_num = allPaths.size();
            }
        }

        void P2PGraph::findAllPath(  unsigned int id  ,Path p)
        {
            unsigned int curr = id;
            while ( sub_graph[curr].tos.size() == 1 && curr != target )
            {
                if( ! p.AddEdge(curr,sub_graph[curr].length, sub_graph[curr].cov,sub_graph[curr].barcode_cov) )
                {
                    //path_num = -1 ; //circle and can't solved  A-B-C-A-B-C-A ....
                    return ;
                }
                curr = * sub_graph[curr].tos.begin();
            }
            if( path_num < 0 )
                return ;
            if ( curr == target )
            {
                allPaths.push_back(p);
            }
            else
            {
                assert( sub_graph[curr].tos.size() >1 );
                if( ! p.AddEdge(curr,sub_graph[curr].length, sub_graph[curr].cov,sub_graph[curr].barcode_cov) )
                {
                    //path_num = -1 ; // circle can solve A-B-C-A-D-E ...
                    //just delete this path and return so it will try another path
                    return ;
                }
                for( const auto & i : sub_graph[curr].tos )
                {
                    findAllPath(i ,p);
                    if( path_num < 0 )
                        return ;
                }
            }
        }

        void P2PGraph::ScoreAllPath()
        {
            if( path_num > 1 )
            {
                for(auto & i : allPaths )
                {
                    i.CalcCov();
                }
                std::sort( allPaths.rbegin() ,allPaths.rend() );
            }
        }

        void P2PGraph::CleanAndSavePath()
        {
            if( path_num > 0 )
            {
                final_path = * allPaths.begin();
                allPaths.clear();
            }
        }


        void P2PGraph::GeneratePath()
        {
            findAllPath();
            ScoreAllPath();
            CleanAndSavePath();
        }

        ////        void P2PGraph::ConstructSubP2pGraph()
        //       {
        //           ConstructSubP2pGraph(root_graph);
        //       }

        void GraphEA_withBarcode::LoadBarcodeOnConfig(const std::string & file)
        {
            auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file);
            std::string line;
            while(!std::getline(*in,line).eof())
            {
                long readId ;
                unsigned int contigId, pos , barcode;
                char dir;
                std::istringstream ist(line);
                ist>>readId>>contigId>>pos>>dir>>barcode;
                if( barcode_on_contig[contigId].find(barcode)
                        == barcode_on_contig[contigId].end())
                    barcode_on_contig[contigId][barcode] =1 ;
                else
                    barcode_on_contig[contigId][barcode] ++ ;
            }
            delete in;
        }

    }// namespace stLFR
}// namespace BGIQD
