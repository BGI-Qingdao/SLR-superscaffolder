#include "algorithm/graph/GraphBasic.h"
#include "common/flags/flags.h"

namespace BGIQD {
    namespace stLFR {
        struct ContigNode : public BGIQD::GRAPH::IGraphNodeBasic<unsigned int , long >
        {
            int contigLen ;
        };

        struct PEEdge : public BGIQD::GRAPH::IDigraphEdgeBase<unsigned int ,long>
        {
            int count ;
            int len;
            FLAGS_INT 
            ADD_A_FLAG( 0 , R2R );

            std::string AttrString() const
            {
                std::ostringstream ost;
                ost<<" count= "<<count<<" len= "<<len;//<<" flags= "<<flags;
                return ost.str();
            }
            std::string ToString() const
            {
                std::ostringstream ost;
                ost<<from<<"\t->\t"<<to<<" [ "<<AttrString()<<" ]";
                return ost.str();
            }
        };

        struct ContigPEGraph : public BGIQD::GRAPH::ListDigraph<ContigNode , PEEdge >
        {
            typedef BGIQD::GRAPH::ListDigraph<ContigNode , PEEdge > Basic;
            // Must guarantee contigId +1 is the compelete reverse contig if this one.
            void AddNode( unsigned int contigId , int len)
            {
                Node tmp ;
                tmp.contigLen = len ;
                tmp.id = contigId ;
                Basic::AddNode(tmp);
                tmp.id = contigId +1 ;
                Basic::AddNode(tmp);
            }

            void AddEdge( unsigned int from , unsigned int to ,int len, int count ) 
            {
                if( ! ( Basic::HasNode( from ) && Basic::HasNode(to) ) )
                    return ;
                Edge edge;
                edge.from = from ;
                edge.to = to ; 
                edge.count = count ;
                edge.len = len ;
                Basic::AddEdge(edge);
            }
        };
    }
}
