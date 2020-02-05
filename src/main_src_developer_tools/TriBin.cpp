#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/error/Error.h"
#include "common/freq/freq.h"
#include "common/stl/setHelper.h"
#include "stLFR/CBB.h"

#include "soap2/fileName.h"
#include <set>
#include <vector>
#include <stack>
#include <tuple>
#include <sstream>
#include <cassert>

struct AppConfig {

    BGIQD::SOAP2::FileNames fNames;

    struct BarcodeTypeInfo {
        int A_only;
        int B_only;
        int C_only;
        int AB;
        int AC;
        int BC;
        int ABC;
    };

    // Done
    BarcodeTypeInfo GetBarcodeDetail(int ref_id, int start_index ) const {
        BarcodeTypeInfo ret ;
        const auto & bbs = b2bs.at(ref_id);
        assert(  bbs.find(start_index) != bbs.end() ) ;
        assert(  bbs.find(start_index+1) != bbs.end() ) ;
        assert(  bbs.find(start_index+2) != bbs.end() ) ;

        auto collection2set = [](const BGIQD::Collection::Collection<int> & c ) {
            std::set<int> ret ;
            for( const auto & pair : c.elements ) ret.insert(pair.first) ;
            return ret ;
        };
        const auto & Ac = bbs.at(start_index).collections ;
        std::set<int> A = collection2set(Ac) ;

        const auto & Bc = bbs.at(start_index+1).collections ;
        std::set<int> B = collection2set(Bc) ;

        const auto & Cc = bbs.at(start_index+2).collections ;
        std::set<int> C = collection2set(Cc) ;

        auto AB = BGIQD::STL::set_common(A,B);
        auto AC = BGIQD::STL::set_common(A,C);
        auto BC = BGIQD::STL::set_common(B,C);
        auto ABC =  BGIQD::STL::set_common(AB ,BC ) ;
        auto A_only = BGIQD::STL::set_diff_in_s1(A,B);
        A_only = BGIQD::STL::set_diff_in_s1(A_only,C);
        auto B_only = BGIQD::STL::set_diff_in_s1(B,A);
        A_only = BGIQD::STL::set_diff_in_s1(B_only,C);
        auto C_only = BGIQD::STL::set_diff_in_s1(C,B);
        C_only = BGIQD::STL::set_diff_in_s1(C_only,A);

        ret.A_only = A_only.size();
        ret.B_only = B_only.size();
        ret.C_only = C_only.size();
        ret.AB     = AB.size();
        ret.AC     = AC.size();
        ret.BC     = BC.size();
        ret.ABC    = ABC.size();
        return ret ;
    }
    std::map<int , std::map<int,BGIQD::stLFR::BarcodeOnBin> > b2bs;
    // Done
    void LoadBarcodeOnBin()
    {
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(
                fNames.BarcodeOnBin());
        if(! in )
            FATAL( "failed to open xxx.barcodeOnBin to read !");
        std::string line;
        auto add_data = [&] ( const std::string & line)
        {
            BGIQD::stLFR::BarcodeOnBin b2b ;
            b2b.InitFromString(line);
            b2bs[b2b.contigId][b2b.binId] = b2b ;
        };
        BGIQD::FILES::FileReaderFactory::EachLine(*in,add_data);
        assert(b2bs.size()==1);
        delete in;
    }
    void PrintTribinInfos() {
        std::cout<<"ref_id,start_pos,n1,n2,n3,n4,n5,n6.n7\n";
        for( const auto & pair : b2bs ) { 
            int ref_id = pair.first ;
            const auto & bbs = pair.second ;
            for ( int i = 0 ; i < ( (int)bbs.size() )- 2 ; i++ ) {
                assert(  bbs.find(i) != bbs.end() ) ;
                auto info = GetBarcodeDetail(ref_id,i);
                const auto & b2b = bbs.at(i) ;
                std::cout<<ref_id<<','<<b2b.start
                    <<','<<info.A_only
                    <<','<<info.B_only
                    <<','<<info.C_only
                    <<','<<info.AB
                    <<','<<info.BC
                    <<','<<info.AC
                    <<','<<info.ABC<<'\n';
            }
        }
    }
}config;

int main(int argc , char **argv )
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string , prefix , "prefix ");
    END_PARSE_ARGS;
    config.fNames.Init(prefix.to_string());
    config.LoadBarcodeOnBin();
    config.PrintTribinInfos();
    return 0 ;
}
