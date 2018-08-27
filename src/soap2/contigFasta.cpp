#include "contigFasta.h"
#include "common/files/file_reader.h"
#include "common/error/Error.h"
#include "biocommon/seq/tool_func.h"
#include <sstream>
#include <cassert>

namespace BGIQD {
    namespace SOAP2 {

        void ContigFastA::Init(const std::string & line , int k) {
            id = -1 ;
            length = 0 ;
            cov = 0;
            tip = 0;
            flag = 0;
            K = "";
            linear = "";
            sscanf(line.c_str() ,">%d length %d cvg_%f_tip_%d",&id, &length , &cov , &tip);
            length -= k ;
        }

        ContigFastA ContigFastA::ReverseCompelete() const 
        {
            ContigFastA ret = *this;
            ret.K = "";
            ret.linear = "";
            ret.flag = 0 ;
            ret.id++;
            std::string base = K+linear;
            ret.AddSeq(BGIQD::SEQ::seqCompleteReverse(base) , K.size());
            return ret ;
        }

        void ContigFastA::AddSeq ( const std::string & line, int k)
        {
            if( IsKSet() )
            {
                linear += line;
            }
            else
            {
                int need_k = k - K.size();
                int k_len = (int)line.size() < need_k ? line.size() : need_k ;
                int left = (int)line.size() > need_k ? ( (int)line.size() - need_k ) : 0 ;
                K += line.substr(0,k_len);
                if( left ) 
                {
                    linear += line.substr(k_len, left);
                }
                if( k_len == need_k )
                {
                    MarkSetK();
                }
            }
        }

        std::string ContigFastA::ToString(int new_id ,const std::string & marker ) const
        {
            std::ostringstream ost;
            //">%d length %d cvg_%f_tip_%d"
            ost<<'>'<<new_id<<" length "<<length+K.size()<<" cvg_"<<cov<<"_tip_"<<tip<<'\t';
            ost<<marker<<"\n";

            int index = 0 ;
            for( size_t i = 0 ; i< K.size(); i++)
            {
                ost<<K[i];
                index ++ ;
                if( index % 100  == 0 )
                {
                    ost<<'\n';
                }
            }

            for( size_t i = 0 ; i< linear.size(); i++)
            {
                ost<<linear[i];
                index ++ ;
                if( index % 100  == 0&& (index != (int)linear.size() -1 ) )
                {
                    ost<<'\n';
                }
            }
            return ost.str();
        }
        std::string ContigFastA::ToString() const 
        {
            std::ostringstream ost;
            //">%d length %d cvg_%f_tip_%d"
            ost<<'>'<<id<<" length "<<length+K.size()<<" cov_"<<cov<<"_tip_"<<tip<<'\n';
            size_t index = 0 ;
            for( size_t i = 0 ; i< K.size(); i++)
            {
                ost<<K[i];
                index ++ ;
                if( index % 100  == 0 )
                {
                    ost<<'\n';
                }
            }

            for( size_t i = 0 ; i< linear.size(); i++)
            {
                ost<<linear[i];
                index ++ ;
                if( index % 100  == 0 && (index != linear.size() -1 ))
                {
                    ost<<'\n';
                }
            }
            return ost.str();
        }

        ContigFastA ContigFastAMap::MergeContig(const std::vector<std::string> & line)
        {
            assert(line.size() > 1);
            unsigned int start = std::stoul(line[0]);
            auto & first = contigs.at(start);
            ContigFastA ret = first;
            first.MarkMerge();
            if(! first.IsBase() )
            {
                contigs.at(start-1).MarkMerge();
            }

            float cov = ret.cov * ret.length + ret.cov * K ;

            for( int i = 1 ; i < (int)line.size() ; i++ )
            {
                unsigned int next_id = std::stoul(line[i]);
                if( contigs.find( next_id ) == contigs.end() )
                {
                    std::cerr<<"WARN : invalid contig id : "<<next_id<<std::endl;
                    continue;
                }
                auto & next = contigs.at(next_id) ;
                ret.linear += next.linear ;
                ret.length += next.length ;
                cov += next.cov * next.length ;
                next.MarkMerge();

                if(! next.IsBase() )
                {
                    contigs.at(next_id-1).MarkMerge();
                }
            }

            ret.id = nextContigNum() ;
            ret.cov = cov / (ret.length + K ) ;
            return ret;
        }

        ContigFastA ContigFastAMap::MergeContig(const std::vector<unsigned  int> & line)
        {
            assert(line.size() > 1);
            unsigned int start = line[0];
            auto & first = contigs.at(start);
            ContigFastA ret = first;
            first.MarkMerge();
            if(! first.IsBase() )
            {
                contigs.at(start-1).MarkMerge();
            }

            float cov = ret.cov * ret.length + ret.cov * K ;

            for( int i = 1 ; i < (int)line.size() ; i++ )
            {
                unsigned int next_id = line[i];
                if( contigs.find( next_id ) == contigs.end() )
                {
                    std::cerr<<"WARN : invalid contig id : "<<next_id<<std::endl;
                    continue;
                }
                auto & next = contigs.at(next_id) ;
                ret.linear += next.linear ;
                ret.length += next.length ;
                cov += next.cov * next.length ;
                next.MarkMerge();

                if(! next.IsBase() )
                {
                    contigs.at(next_id-1).MarkMerge();
                }
            }

            ret.id = nextContigNum() ;
            ret.cov = cov / (ret.length + K ) ;
            return ret;
        }

        void ContigFastAMap::LoadContig(const std::string & file)
        {
            auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file);
            if( in == NULL )
                FATAL( " open prefix.contig for read failed !!! ");
            std::string line ;
            ContigFastA tmp ;
            tmp.length = 0; 
            while(!std::getline(*in,line).eof())
            {
                if( line[0] == '>')
                {
                    if( tmp.length > 0 )
                    {
                        assert(tmp.IsSeqComplete(K));
                        contigs[tmp.id] = tmp;
                        contigs[tmp.id].MarkBase() ;
                        if(tmp.id > maxContig )
                        {
                            maxContig = tmp.id;
                        }
                        tmp.length = 0 ;
                    }
                    tmp.Init(line,K);
                }
                else
                {
                    tmp.AddSeq(line,K);
                }
            }
            if( tmp.length > 0 )
            {
                assert(tmp.IsSeqComplete(K));
                contigs[tmp.id] = tmp;
                contigs[tmp.id].MarkBase() ;
                if(tmp.id > maxContig )
                {
                    maxContig = tmp.id;
                }
                tmp.length = 0 ;
            }
            delete in;
        }

        void ContigFastAMap::buildCompeleReverse()
        {
            for(unsigned int i = 1 ; i<= maxContig ; i++)
            {
                auto itr = contigs.find( i ) ;
                if( itr == contigs.end() || ! itr->second.IsBase())
                    continue;
                const auto & base = itr->second ;
                if( BGIQD::SEQ::isSeqPalindrome(base.K + base.linear) 
                        ||
                    contigs.find(i+1) != contigs.end() )
                    continue;
                contigs[i+1] = base.ReverseCompelete();
                if( i == maxContig )
                    maxContig ++ ;
            }
        }
    }//namespace SOAP2
}//namespace BGIQD
