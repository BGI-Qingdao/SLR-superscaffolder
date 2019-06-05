#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/error/Error.h"
#include "stLFR/CBB.h"

#include "biocommon/agp/agp.h"
#include <map>
#include <sstream>
#include <fstream>
#include <tuple>

typedef BGIQD::AGP::AGP_Item AGP_Item;
typedef BGIQD::stLFR::ContigBarcodeInfo barcodeOnContigItem ;
const std::string invalid_barcode_name("0_0_0"); 
struct AppConf
{
    //barcode_id ->barcode_seq 
    std::map<int , std::string > b2b;
    void LoadBarcodeIdMap(const std::string & file)
    {
        std::ifstream ifs;
        ifs.open(file);
        if( ! ifs.is_open() )
            FATAL( "can not open barcode_list file to read !");
        std::string line ;
        while(! std::getline( ifs , line ).eof() )
        {
            std::string barcode_name;
            unsigned int barcode_id ;
            std::istringstream ist(line);
            ist>>barcode_name>>barcode_id;
            if(barcode_id != 0 )
                b2b[barcode_id] = barcode_name ;
        }
    }

    std::string GetBarcodeSeq( int barcode_id )
    {
        if( b2b.find( barcode_id ) == b2b.end () )
            return invalid_barcode_name;
        else
            return b2b.at(barcode_id);
    }
    //contig_id ->contig_oid
    std::map<unsigned int , std::string > c2c;
    void LoadContigIdMap(const std::string & file)
    {
        std::ifstream ifs;
        ifs.open(file);
        if( ! ifs.is_open() )
            FATAL( "can not open fake_name file to read !");
        std::string line ;
        while(! std::getline( ifs , line ).eof() )
        {
            std::string contig_name ;
            unsigned int contig_id ;
            std::istringstream ist(line);
            ist>>contig_name>>contig_id;
            c2c[contig_id]= contig_name;
        }
    }

    std::string GetContigName( unsigned int contig_id )
    {
        if( c2c.find(contig_id) == c2c.end() )
            return std::to_string(contig_id);
        else
            return c2c.at(contig_id);
    }
    int invalid_agp ;
    std::map<std::string ,AGP_Item>  contig2scaffold;
    void LoadAGP( const std::string & file )
    {
        std::ifstream ifs;
        ifs.open(file);
        if( ! ifs.is_open() )
            FATAL( "can not open scaff_agp file to read !");
        std::string line ;
        while(! std::getline( ifs , line ).eof() )
        {
            if(!AGP_Item::IsAGPStringValid(line))
            {
                invalid_agp ++ ;
                continue ;
            }
            AGP_Item tmp ;
            tmp.InitFromString(line);
            if(tmp.component_type != AGP_Item::ComponentType::U
                && tmp.component_type !=AGP_Item::ComponentType::N )
            {
                contig2scaffold[tmp.lefta.component_id] = tmp ;
            }
        }
    }
    //         succ? scaffold_id      scaff_pos 
    std::tuple<bool ,std::string, int > GetScaffDetail(const std::string &contig_name , int contig_pos )
    {
        if( contig2scaffold.find(contig_name) == contig2scaffold.end() )
            return std::make_tuple(false ,"",-1);
        const AGP_Item & item = contig2scaffold.at(contig_name);
        int scaff_pos = 0 ;
        if( item.lefta.orientation == "+" )
            scaff_pos = contig_pos - item.lefta.component_beg + item.object_beg ;
        else
            scaff_pos = -( contig_pos - item.lefta.component_beg ) + item.object_end ;
        return std::make_tuple( true ,item.object,  scaff_pos );
    }

    int drop  ;

    void Init()
    {
        drop = 0 ;
        invalid_agp = 0 ;
    }
    void Print( unsigned int contig_id , int contig_pos , int barcode_id)
    {
        std::string barcode_name = GetBarcodeSeq(barcode_id);
        if( barcode_name == invalid_barcode_name ) 
            return ;
        std::string contig_name = GetContigName( contig_id);
        std::string scaff_name ; int scaff_pos ; bool check_succ ;
        std::tie(check_succ, scaff_name,scaff_pos) 
            = GetScaffDetail(contig_name,contig_pos );
        if(check_succ )
            std::cout<<barcode_name<<'\t'
                <<barcode_id<<'\t'
                <<contig_name<<'\t'
                <<contig_id<<'\t'
                <<contig_pos<<'\t'
                <<scaff_name<<'\t'
                <<scaff_pos<<'\n' ;
        else
            drop ++ ;
    }
}config;

int main(int argc , char **argv )
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string , barcode_list , "input barcodeList file ");
    DEFINE_ARG_REQUIRED(std::string , barcodeOnContig , "input barcodeOnContig file ");
    DEFINE_ARG_REQUIRED(std::string , scaff_agp , "input scaff_agp file ");
    DEFINE_ARG_OPTIONAL(std::string , fake_name, "input fakesoap.name2index.map.txt file ","");
    END_PARSE_ARGS;
    config.Init();
    if( fake_name.to_string() != "")
    {
        std::cerr<<" Loading "<< fake_name.to_string()<<" now ... "<<std::endl;
        config.LoadContigIdMap( fake_name.to_string());
    }
    std::cerr<<" Loading "<< barcode_list.to_string()<<" now ... "<<std::endl;
    config.LoadBarcodeIdMap(barcode_list.to_string());
    std::cerr<<" Loading "<< scaff_agp.to_string()<<" now ... "<<std::endl;
    config.LoadAGP(scaff_agp.to_string());
    std::cerr<<" Processing "<< barcodeOnContig.to_string()<<" now ... "<<std::endl;
    std::ifstream barcodeOnContig_file;
    barcodeOnContig_file.open(barcodeOnContig.to_string());
    if( ! barcodeOnContig_file.is_open() )
        FATAL( "can not open barcodeOnContig file to read !");
    std::string line ;
    while(! std::getline( barcodeOnContig_file , line ).eof() )
    {
        barcodeOnContigItem tmp ;
        tmp.InitFromString(line);
        for( const auto & pos_pair : tmp.barcodesOnPos )
        {
            for( unsigned int barcode_id : pos_pair.second )
            {
                config.Print(tmp.contig_id,
                        pos_pair.first,
                        barcode_id);
            }
        }
    }
    std::cerr<<"total number of drop contigs is: "<<config.drop<<std::endl;
    std::cerr<<"total number of invalid agp is: "<<config.invalid_agp<<std::endl;
    std::cerr<<"All Done !!! Exit now ... "<<std::endl;
    return 0 ;
}
