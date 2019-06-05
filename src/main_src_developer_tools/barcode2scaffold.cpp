#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/error/Error.h"
#include "stLFR/CBB.h"

#include "biocommon/agp/agp.h"
#include <map>
#include <vector>
#include <sstream>
#include <fstream>

typedef BGIQD::AGP::AGP_Item AGP_Item;
typedef BGIQD::stLFR::ContigBarcodeInfo barcodeOnContigItem ;
struct AppConf
{
    //barcode_id ->barcode_seq 
    std::map<int , std::string > b2b;
    void LoadBarcodeIdMap(const std::string & file)
    {

    }
    std::string GetBarcodeSeq( int barcode_id )
    {

    }
    //contig_id ->contig_oid
    std::map<unsigned int , std::string > c2c;
    void LoadContigIdMap(const std::string & file)
    {

    }

    std::string GetContigName( unsigned int contig_id )
    {

    }

    std::map<std::string ,AGP_Item>  contig2scaffold;
    void LoadAGP( const std::string & file )
    {

    }
    // scaffold_id      scaff_pos 
    std::pair<std::string, int > GetScaffDetail(const std::string &contig_name , int contig_pos ) 
    {

    }

    void Print( unsigned int contig_id , int contig_pos , int barcode_id)
    {
        std::string barcode_name = GetBarcodeSeq(barcode_id);
        std::string contig_name = GetContigName( contig_id);
        auto scaff_detail = GetScaffDetail(contig_name,contig_pos );
        std::cout<<barcode_name<<'\t'
                 <<barcode_id<<'\t'
                 <<contig_name<<'\t'
                 <<contig_id<<'\t'
                 <<scaff_detail.first<<'\t'
                 <<scaff_detail.second<<'\n' ;
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
    std::cerr<<"All Done !!! Exit now ... "<<std::endl;
    return 0 ;
}
