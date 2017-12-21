#ifndef __CONTIG_BARCODE_H__
#define __CONTIG_BARCODE_H__

#include <vector>
#include <map>

namespace BGIQD {
namespace JOB01 {

typedef std::vector<int> barcodeList;

// 1 means contig pos , 2 means contig id.
typedef std::tuple<int,int> contigPosInfo;
typedef std::vector<contigPosInfo> contigList;

// In barcodeList , means contig pos , 2 means barcode ids
typedef std::tuple< int , barcodeList> barcodePosInfo;
typedef std::vector<barcodePosInfo> barcodePosList;

typedef std::map<int , barcodeList> refBarcodeInfo;

typedef std::map<int , contigList> refContigInfo;

typedef std::map<int,barcodePosList> contigBarcodeInfo;

typedef std::map<std::string , int> barcodeNum;

struct BarcodeNum
{
    BarcodeNum() : next(1) {}
    int barcode2num(const std::string & str);
    void save(const std::string & file) const;
    barcodeNum data;
    private:
        int next ;
};

void initLog();

void loadRefBarcodeInfo(const std::string & file , refBarcodeInfo & data);

void loadRefContigInfo( const std::string & file , refContigInfo & data);

void generateConrigBarcodeInfo( const refBarcodeInfo & i_b ,
                                const refContigInfo & i_c,
                                contigBarcodeInfo & data );

void printContigBarcodeInfo( const contigBarcodeInfo & data , const std::string & file );

}//JOB01
}//BGIQD
#endif //__CONTIG_BARCODE_H__
