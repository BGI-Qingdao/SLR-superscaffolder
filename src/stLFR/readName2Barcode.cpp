#include "readName2Barcode.h"
#include "stringtools.h"

namespace BGIQD{
namespace stLFR{
std::string readName2Barcode( const std::string & read ) 
{
    //read looks like : CL200034461L1C015R086_gg2186#531_921_518
    return BGIQD::STRING::split(read,"#").at(1);
}

}
}
