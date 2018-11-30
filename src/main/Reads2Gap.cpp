#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/error/Error.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/freq/freq.h"

#include "biocommon/fasta/fasta.h"
#include "biocommon/fastq/fastq.h"

#include "design_pattern/Publish_Subscribe.h"
#include "design_pattern/Filter.h"

#include "soap2/contigIndex.h"
#include "soap2/fileName.h"

#include "stLFR/EasySam.h"
#include "stLFR/StringIdCache.h"

#include <fstream>
#include <map>
#include <vector>
#include <set>

typedef BGIQD::FASTA::ScaffSplitGapHead GapHead;
typedef BGIQD::FASTA::Fasta<GapHead> GapItem;
typedef BGIQD::FASTA::FastaReader<GapItem> GapReader;

typedef BGIQD::FASTQ::Fastq<BGIQD::FASTQ::stLFRHeader>  stLFRRead ;
typedef BGIQD::FASTQ::FastqReader<stLFRRead> stLFRReader ;

// For reads data flow .
typedef BGIQD::DesignPattern::ISubscriber<stLFRRead> ISub;
typedef BGIQD::DesignPattern::IPublisher<stLFRRead> IPub;
typedef BGIQD::DesignPattern::Filter<std::string, IPub> Filter;

// For pe_single data flow.
typedef BGIQD::DesignPattern::ISubscriber<BGIQD::EASY_SAM::PE_Single>  PESingleSub;
typedef BGIQD::DesignPattern::IPublisher<BGIQD::EASY_SAM::PE_Single>  PESinglePub;
typedef BGIQD::DesignPattern::Filter<unsigned int , PESinglePub> PESingleFilter;

// For pe_both data flow .
typedef BGIQD::DesignPattern::ISubscriber<BGIQD::EASY_SAM::PEInfo>  PEBothSub;
typedef BGIQD::DesignPattern::IPublisher<BGIQD::EASY_SAM::PEInfo>  PEBothPub;
typedef BGIQD::DesignPattern::Filter<unsigned int , PEBothPub> PEBothFilter;

/////////////////////////////////
//
// Global variable
//
/////////////////////////////////

struct AppConfig;

struct PEGap;

struct TrunkGap;

struct ReadNameRegisterPEBoth;

struct ReadNameRegister;

BGIQD::SOAP2::ContigIndexMap contigIndexs;

Filter read_name_filter;

Filter barcode_filter ;

PESingleFilter pe_single_filter;

PEBothFilter pe_both_filter;

BGIQD::stLFR::IdStringCache read_num_2_str;

PEGap * pe_gaps ;

TrunkGap * trunk_gaps;

BGIQD::SOAP2::FileNames fNames;

BGIQD::LOG::logger loger;

struct GapContigs
{
    BGIQD::SOAP2::ContigIndex c1 ;

    BGIQD::SOAP2::ContigIndex c2 ;

    void Init(unsigned int p , unsigned int n )
    {
        c1 = contigIndexs.GetContigIndex(p);
        c2 = contigIndexs.GetContigIndex(n);
    }

    bool IsContain(const BGIQD::EASY_SAM::PE_Single & match_info)
    {
        if( match_info.match_reverse1 )
        {
            if( match_info.contig1 == c1.contig )
                return match_info.pos_1bp1 <= c1.length && match_info.pos_1bp1>= 100 ;
            else if ( match_info.contig1 == c2.contig )
                return match_info.pos_1bp1 <= c2.length && match_info.pos_1bp1>= 100 ;
        }
        else
        {
            if( match_info.contig1 == c1.contig )
                return match_info.pos_1bp1 +100 <= c1.length && match_info.pos_1bp1>= 1 ;
            else if ( match_info.contig1 == c2.contig )
                return match_info.pos_1bp1 +100 <= c2.length && match_info.pos_1bp1>= 1 ;
        }
        return false ;
    }
};

struct ScaffPackage
{
    int id ;
    std::ofstream * scaff ; 
    std::ofstream * r1; 
    std::ofstream * r2;
    const static int reads_buffer_size = 1024;

    void Init(const std::string & id )
    {
        std::string idd = id ;
        for( auto & x :idd )
        {
            if( std::isalnum(x) )
                continue ;
            else
                x = '_';
        }
        std::string prefix = "tmp_filler_" + (idd);
        // open 3 files
        scaff = new std::ofstream(prefix+".scaff.fa");
        r1 = new std::ofstream(prefix+".r1.fq");
        r2 = new std::ofstream(prefix+".r2.fq");
        index = 0 ;
        r1_flag = true ;
    }

    void AddScaff(const GapItem & item )
    {
        // save scaffold
        (*scaff)<<item.head.Head()<<'\n';
        (*scaff)<<item.seq.Seq(100);
        delete scaff ;
        scaff = NULL ;
    }

    stLFRRead buffer[reads_buffer_size];
    int index ;

    void CleanBufferR1()
    {
        for( int i = 0 ; i < index ; i ++ )
        {
            (*r1)<<buffer[i].head.Head()<<'\n';
            (*r1)<<buffer[i].seq.Seq(10000000);
        }
        index = 0 ;
    }

    void AddR1( const stLFRRead & read )
    {
        buffer[index] = read ;
        index ++ ;
        if ( index >= reads_buffer_size )
            CleanBufferR1() ;
    }

    void CleanBufferR2()
    {
        for( int i = 0 ; i < index ; i ++ )
        {
            (*r2)<<buffer[i].head.Head()<<'\n';
            (*r2)<<buffer[i].seq.Seq(10000000);
        }
        index = 0 ;
    }

    void AddR2( const stLFRRead & read )
    {
        buffer[index] = read ;
        index ++ ;
        if ( index >= 10000 )
            CleanBufferR2() ;
    }

    void AddRead(  const stLFRRead & read )
    {
        if( r1_flag )
            AddR1( read );
        else
            AddR2( read );
    }

    ~ScaffPackage()
    {
        if ( scaff ) delete  scaff ;
        if ( r1 ) delete r1;
        if ( r2 ) delete r2;
    }

    void EndR1()
    {
        CleanBufferR1();
        r1_flag = false ;
    }

    void EndR2()
    {
        CleanBufferR2();
    }
    bool r1_flag ;
};

struct ReadNameRegister : public PESingleSub
{
    ISub* item ;

    virtual void update_msg( const BGIQD::EASY_SAM::PE_Single & info )
    {
        assert(item != NULL );
        if( info.read1 %2 == 0 )
        {
            std::string read_name = read_num_2_str.Id(info.read1-1);
            read_name_filter.watch( read_name , item);
        }
        else
        {
            std::string read_name = read_num_2_str.Id(info.read1);
            read_name_filter.watch( read_name , item);
        }
    }
};

struct ITrunkSub : public ISub
{
    virtual void CheckMatch( const BGIQD::EASY_SAM::PE_Single & info , bool p) = 0 ;
};

struct ReadNameRegisterPEBoth : public PEBothSub 
{
    ISub* item ;
    virtual void update_msg( const BGIQD::EASY_SAM::PEInfo& info )
    {
        unsigned int read = info.read1 < info.read2 ? info.read1 : info.read2 ;
        std::string read_name = read_num_2_str.Id(read);
        //PE share the same read_name
        read_name_filter.watch( read_name , item);
        // TODO : need split make r1 & r2
        dynamic_cast<ITrunkSub*>(item)->CheckMatch(info.PInfo(), true);
        dynamic_cast<ITrunkSub*>(item)->CheckMatch(info.EInfo(), false);
    }
};

struct PEGap :  public GapContigs , public ScaffPackage , public ISub
{
    ReadNameRegister rgt;

    void Init(const GapItem & item )
    {
        GapContigs::Init(item.head.prev_base_contig ,
                item.head.next_base_contig);
        ScaffPackage::Init(item.head.Head());
        rgt.item = this ;
        pe_single_filter.watch(c1.contig ,&rgt);
        pe_single_filter.watch(c2.contig ,&rgt);
    }

    virtual void update_msg( const stLFRRead & msg ) final
    {
        AddRead(msg);
    }
};

struct TrunkGap : public GapContigs , public ScaffPackage , public ITrunkSub
{
    ReadNameRegister rgt;
    ReadNameRegisterPEBoth rgtb;
    void Init(const GapItem & item)
    {
        GapContigs::Init(item.head.prev_base_contig ,
                item.head.next_base_contig);
        ScaffPackage::Init(item.head.Head());
        rgt.item = this ;
        pe_single_filter.watch(c1.contig ,&rgt);
        pe_single_filter.watch(c2.contig ,&rgt);

        rgtb.item = this ;
        pe_both_filter.watch(c1.contig ,&rgtb);
        pe_both_filter.watch(c2.contig ,&rgtb);
    }

    std::set<std::string> barcodes ;

    std::set<std::string> maped_reads_str_r1;

    std::set<std::string> maped_reads_str_r2;

    std::vector<stLFRRead> maped_reads ;

    std::vector<stLFRRead> unmaped_reads ;

    std::string prev_read_name ;

    void CheckMatch( const BGIQD::EASY_SAM::PE_Single & info , bool p)
    {
        if( IsContain(info) )
        {
            unsigned int read = info.read1 % 2 == 1 ? info.read1  : info.read1 -1 ;
            if( p ) 
                maped_reads_str_r1.insert(read_num_2_str.Id(read));
            else
                maped_reads_str_r2.insert(read_num_2_str.Id(read));
        }
    }

    virtual void update_msg( const stLFRRead & msg ) final
    {
        // 2 filter may pass one read into here twice .
        if( msg.head.readName == prev_read_name )
            return ;
        prev_read_name =  msg.head.readName ;
        // PE not barcode keep normal ;
        if( barcodes.find(msg.head.barcode_str) ==  barcodes.end() )
            AddRead(msg);
        // barcode reads save first.
        else
        {
            if( r1_flag )
            {
                if( maped_reads_str_r1.find( msg.head.readName ) == maped_reads_str_r1.end())
                    maped_reads.push_back(msg);
                else
                    unmaped_reads.push_back(msg);
            }
            else
            {
                if( maped_reads_str_r2.find( msg.head.readName ) == maped_reads_str_r2.end())
                    maped_reads.push_back(msg);
                else
                    unmaped_reads.push_back(msg);
            }
        }
    }

    // random fake maped reads & unmaped reads into PE reads
    void EndAllReads()
    {
        for(size_t i = 0 ; i < unmaped_reads.size() ; i++ )
        {
            AddR1(maped_reads[i%maped_reads.size()] );
        }
        CleanBufferR1();
        for(size_t i = 0 ; i < unmaped_reads.size() ; i++ )
        {
            AddR2(unmaped_reads[i] );
        }
        CleanBufferR2();
    }
};

struct AppConfig
{
    void Init(
            const std::string & prefix
            ,const std::string & read1
            ,const std::string & read2
            )
    {
        r1 = read1 ;
        r2 = read2 ;
        fNames.Init(prefix);
        BGIQD::LOG::logfilter::singleton().get("Reads2Gap",BGIQD::LOG::loglevel::INFO, loger);
    }

    std::vector<GapItem> gap_buffer ;

    void LoadGaps()
    {
        BGIQD::LOG::timer t( loger,"LoadGaps");
        auto in  = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fNames.scaff_gap2filler_seqs());
        if( in == NULL )
            FATAL(" failed to open xxx.scaff_gap2filler_seqs for read!!! ");
        GapReader reader ;
        reader.LoadAllFasta(*in , gap_buffer);
        delete in;
    }

    int pe_gap_num ;

    int trunk_gap_num ;

    std::map<std::string , ISub * > gapname_ptr_map;

    void BuildGaps()
    {
        BGIQD::LOG::timer t( loger,"BuildGaps");
        pe_gap_num = 0 ;
        trunk_gap_num = 0 ;
        for( const auto & i : gap_buffer )
        {
            if( i.head.gap_type != GapHead::GapType::TRUNK )
                pe_gap_num ++ ;
            else
                trunk_gap_num ++ ;
        }

        pe_gaps = new PEGap[pe_gap_num] ;
        trunk_gaps = new TrunkGap[trunk_gap_num];
        pe_gap_num = 0 ;
        trunk_gap_num = 0 ;
        for( const auto & i : gap_buffer )
        {
            if( i.head.gap_type != GapHead::GapType::TRUNK )
            {
                (pe_gaps[pe_gap_num]).Init(i);
                (pe_gaps[pe_gap_num]).AddScaff(i);
                pe_gap_num ++ ;
                gapname_ptr_map[i.head.Head()] = &(pe_gaps[pe_gap_num]);
            }
            else
            {
                (trunk_gaps[trunk_gap_num]).Init(i);
                (trunk_gaps[trunk_gap_num]).AddScaff(i);
                trunk_gap_num ++ ;
                gapname_ptr_map[i.head.Head()] = &(trunk_gaps[trunk_gap_num]);
            }
        }

    }

    void LoadPESingles()
    {
        BGIQD::LOG::timer t( loger,"LoadPESingles");
        auto in = BGIQD::FILES::FileReaderFactory::
            GenerateReaderFromFileName(fNames.pe_singles());
        if(in == NULL)
            FATAL(" failed to open xxx.pe_singles for read!!! ");
        std::string line ;
        while( ! std::getline(*in, line).eof() )
        {
            BGIQD::EASY_SAM::PE_Single tmp ;
            tmp.InitFromString(line);
            pe_single_filter.notify_with_msg(tmp.contig1 , tmp);
        }
        delete in ;
    }

    void LoadPEBoths()
    {
        BGIQD::LOG::timer t( loger,"LoadPEBoths");
        auto in = BGIQD::FILES::FileReaderFactory::
            GenerateReaderFromFileName(fNames.pe_boths());
        if(in == NULL)
            FATAL(" failed to open xxx.pe_boths for read!!! ");
        std::string line ;
        while( ! std::getline(*in, line).eof() )
        {
            BGIQD::EASY_SAM::PEInfo tmp ;
            tmp.InitFromString(line);
            assert(tmp.contig1 == tmp.contig2);
            pe_both_filter.notify_with_msg(tmp.contig1 , tmp);
        }
        delete in ;
    }
    void LoadContigIndex()
    {
        BGIQD::LOG::timer t( loger,"LoadContigIndex");
        auto in = BGIQD::FILES::FileReaderFactory::
            GenerateReaderFromFileName(fNames.ContigIndex());
        if(in == NULL)
            FATAL(" failed to open xxx.contigIndex for read!!! ");
        contigIndexs.LoadContigIndexs(*in);
        delete in ;
        contigIndexs.BuildReverseCompleteContigs();
    }

    void LoadRead2Num()
    {
        BGIQD::LOG::timer t( loger,"LoadRead2Num");
        read_num_2_str.preload = true ;
        read_num_2_str.LoadStringIdCache(fNames.readNameList());
    }

    void LoadBarcodeOnGaps()
    {
        BGIQD::LOG::timer t( loger,"LoadBarcodeOnGaps");
        auto in = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(fNames.barcodeOnGaps());
        if( in == NULL )
          FATAL(" failed to open xxx.barcodeOnGaps for read !!! ");
        std::string line ;
        std::string gap_name ;
        while( ! std::getline(*in, line).eof() )
        {
            if(line[0] == '>')
                gap_name = line ;
            else
                barcode_filter.watch( line , gapname_ptr_map[gap_name] );
        }
        delete in ;
    }


    void ParseAReadFile(std::string r)
    {
        BGIQD::LOG::timer t( loger,"ParseAReadFile "+r);
        auto in  = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(r);
        if( in == NULL )
        {
            std::cerr<<"open file "<<r<<" for read error"<<std::endl;
        }
        if( in == NULL )
            FATAL(" failed to open r for read!!! ");

        stLFRReader reader;
        stLFRRead tmp ;
        while ( reader.LoadNextFastq(*in , tmp) )
        {
            read_name_filter.notify_with_msg(tmp.head.readName,tmp);
            barcode_filter.notify_with_msg(tmp.head.barcode_str,tmp);
        }
    }

    std::string r1;
    std::string r2;
    void ParseReads()
    {
        BGIQD::LOG::timer t( loger,"ParseReads");
        // Parse read1 
        ParseAReadFile(r1);
        for( int i = 0 ; i < pe_gap_num ; i++ )
            (pe_gaps[i]).EndR1();
        for( int i = 0 ; i < trunk_gap_num ; i++ )
            (trunk_gaps[i]).EndR1();

        // Parse read2 
        ParseAReadFile(r2);
        for( int i = 0 ; i < pe_gap_num ; i++ )
            (pe_gaps[i]).EndR2();
        for( int i = 0 ; i < trunk_gap_num ; i++ )
            (trunk_gaps[i]).EndR2();
        for( int i = 0 ; i < trunk_gap_num ; i++ )
            (trunk_gaps[i]).EndAllReads();
    }

}config;

int main(int argc , char **argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string , prefix , "prefix ");
        DEFINE_ARG_REQUIRED(std::string , r1, "read  1");
        DEFINE_ARG_REQUIRED(std::string , r2, "read  2");
    END_PARSE_ARGS

    config.Init(prefix.to_string()
            ,r1.to_string()
            ,r2.to_string()
            );
    config.LoadRead2Num();
    config.LoadContigIndex();
    config.LoadGaps();
    config.BuildGaps();
    config.LoadPESingles();
    config.LoadPEBoths();
    config.LoadBarcodeOnGaps();
    config.ParseReads();
}
