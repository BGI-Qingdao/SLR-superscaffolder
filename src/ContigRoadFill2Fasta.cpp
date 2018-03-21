#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/string/stringtools.h"
#include <cassert>
#include <sstream>
#include <set>

struct GlobalContig
{
    static std::string seqCompleteReverse(const std::string & line)
    {
        std::string ret ;
        ret.resize(line.size(),'N');
        int index = 0;
        for( auto i = line.rbegin() ; i!= line.rend() ; i++)
        {
            if( *i == 'A')
                ret[index++] = 'T';
            if( *i == 'G')
                ret[index++] = 'C';
            if( *i == 'C')
                ret[index++] = 'G';
            if( *i == 'T')
                ret[index++] = 'A';
        }
        return ret;
    }

    static bool isSeqPalindrome(const std::string & line)
    {
        if( line.size() % 2  == 1 )
            return false;
        int len = line.size();
        for( int i = 0 ; i < len /2 ; i++ )
        {
            if( line[i] == 'A' && line[len-i-1] != 'T' )
                return false ;
            if( line[i] == 'T' && line[len-i-1] != 'A' )
                return false ;
            if( line[i] == 'G' && line[len-i-1] != 'C' )
                return false ;
            if( line[i] == 'C' && line[len-i-1] != 'G' )
                return false ;
        }
        return true;
    }

    struct ContigFastA
    {
        unsigned int    id ;
        int             length ;
        float           cov ;
        int             tip;
        std::string     K;
        std::string     linear;
        int flag ;

        void Init(const std::string & line , int k) {
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
        ContigFastA ReverseCompelete() const 
        {
            ContigFastA ret = *this;
            ret.K = "";
            ret.linear = "";
            ret.flag = 0 ;
            ret.id++;
            std::string base = K+linear;
            ret.AddSeq(seqCompleteReverse(base) , K.size());
            return ret ;
        }
        void MarkMerge() { flag  |= 0x1 ; }
        bool IsMerge() const { return flag & 0x1 ; }
        void MarkSetK() { flag |= 0x2 ; }
        bool IsKSet() const { return flag & 0x2 ; }
        void MarkBase() { flag |= 0x4 ; }
        bool IsBase() const { return flag & 0x4; }
        void AddSeq ( const std::string & line, int k)
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

        std::string ToString() const 
        {
            std::ostringstream ost;
            //">%d length %d cvg_%f_tip_%d"
            ost<<'>'<<id<<" length "<<length+K.size()<<" cov_"<<cov<<"_tip_"<<tip<<'\n';
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
                if( index % 100  == 0 )
                {
                    ost<<'\n';
                }
            }
            return ost.str();
        }

        bool IsSeqComplete(int k) const 
        {
            return ( (int)K.size() == k && (int)linear.size() ==length );
        }
    }; //ConfigFasta

    int K;

    std::map<unsigned int , ContigFastA> contigs;

    BGIQD::LOG::logger log;

    unsigned int maxContig;

    unsigned int nextContigNum()
    {
        maxContig += 2 ;
        return maxContig;
    }
    /*
    unsigned int contigId(unsigned int i )
    {
        if( contigs.find( i ) != contigs.end() )
            return i;
        else
        {
            assert( contigs.find(i-1) != contigs.end() );
            return i -1;
        }
    }*/

    ContigFastA MergeContig(const std::vector<std::string> & line)
    {
        assert(line.size() > 1);
        unsigned int start = std::stoul(line[0]);
        auto & first = contigs[start];
        ContigFastA ret = first;
        first.MarkMerge();
        if(! first.IsBase() )
        {
            contigs[start-1].MarkMerge();
        }
        float cov = ret.cov * ret.length + ret.cov * K ;
        for( int i = 1 ; i < (int)line.size() ; i++ )
        {
            unsigned int next_id = std::stoul(line[i]);
            auto & next = contigs[next_id] ;
            ret.linear += next.linear ;
            ret.length += next.length ;
            cov += next.cov * next.length ;
            next.MarkMerge();

            if(! next.IsBase() )
            {
                contigs[next_id-1].MarkMerge();
            }
        }

        ret.id = nextContigNum() ;
        ret.cov = cov / (ret.length + K ) ;
        return ret;
    }

    void LoadContig()
    {
        BGIQD::LOG::timer t(log,"load contigs");
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(contig_file);
        std::string line ;
        GlobalContig::ContigFastA tmp ;
        tmp.length = 0; 
        while(!std::getline(*in,line).eof())
        {
            if( line[0] == '>')
            {
                if( tmp.length > 0 )
                {
                    assert(tmp.IsSeqComplete(K));
                    contigs[tmp.id] = tmp;
                    if(tmp.id > maxContig )
                    {
                        maxContig = tmp.id;
                    }
                }
                tmp.Init(line,K);
            }
            else
            {
                tmp.AddSeq(line,K);
            }
        }
        delete in;
    }

    std::string contig_file;

    std::string contigroadfill;

    void Init(int k , const std::string & f)
    {
        maxContig = 0;
        K = k;
        contig_file = f;
    }


    void buildCompeleReverse()
    {
        for( auto & i : contigs)
        {
            i.second.MarkBase();
        }

        for(unsigned int i = 1 ; i<= maxContig ; i++)
        {
            auto itr = contigs.find( i ) ;
            if( itr == contigs.end() || ! itr->second.IsBase())
                continue;
            const auto & base = itr->second ;
            if( isSeqPalindrome(base.K + base.linear) )
                continue;
            contigs[i+1] = base.ReverseCompelete();
        }
        maxContig ++ ;
    }

} config;

int main(int argc , char ** argv)
{
    START_PARSE_ARGS
    DEFINE_ARG_DETAIL(int , kvalue, 'K',false,"K value");
    DEFINE_ARG_DETAIL(std::string, prefix, 'o',false,"prefix\n\
                                               need   xxx.contig\n\
                                                xxx.contigroadfill");
    DEFINE_ARG_DETAIL(bool , left, 'l',true,"print left contig ? default no");
    DEFINE_ARG_DETAIL(bool , use, 'u',true,"only print used contig ? default no");
    END_PARSE_ARGS

    config.Init( kvalue.to_int() ,prefix.to_string() + ".contig");
    config.LoadContig();
    config.buildCompeleReverse();

    std::string line;

    while(!std::getline(std::cin,line).eof())
    {
        auto items = BGIQD::STRING::split(line,"\t");
        if( ! use.to_bool() )
        {
            auto ret = config.MergeContig(items);
            std::cout<<ret.ToString()<<std::endl;
        }
        else
        {
            std::set<unsigned int> allused;
            for( const auto i : items )
            {
                unsigned int id =std::stoul(i);
                if( config.contigs[id].IsBase() )
                {
                    allused.insert(id);
                }
                else
                {
                    allused.insert(id-1);
                }
            }
            for( const auto i : allused)
            {
                std::cout<<config.contigs[i].ToString()<<std::endl;
            }
        }
    }
    if( left.to_bool() )
    {
        for( const auto & i : config.contigs )
        {
            if( i.second.IsMerge() || ! i.second.IsBase())
                continue;
            std::cout<<i.second.ToString()<<std::endl;
        }
    }
    return 0;
}
