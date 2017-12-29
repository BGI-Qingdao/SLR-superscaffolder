#ifndef __CONTIG_BARCODE_H__
#define __CONTIG_BARCODE_H__

#include <vector>
#include <map>
#include <set>

namespace BGIQD {
namespace JOB01 {

void initLog(const std::string & module);

struct BarcodeNum
{
    typedef std::map<std::string , int> barcodeNum;
    BarcodeNum() : next(1) {}
    int barcode2num(const std::string & str);
    void save(const std::string & file) const;
    barcodeNum data;
    private:
        int next ;
};

/*****************************************************************************
 *   file name      :
 *                      barcode_onref
 *   file format    :
 *                  1 |123_456|...   <== here may have duplicate.
 *
 * **************************************************************************/
typedef std::set<std::string> barcodeSet;
typedef std::map<int , barcodeSet> refBarcodeUniqueInfo;
void loadRefBarcodeUniqueInfo(const std::string & file , refBarcodeUniqueInfo & data);
void saveRefBarcodeUniqueInfo(const std::string & file , const refBarcodeUniqueInfo & data);

/*****************************************************************************
 *   file name      :
 *                      xxx.sam
 *                      use contig align to ref
 *
 * **************************************************************************/
// 1 means contig pos , 2 means contig id.
typedef std::tuple<int,int> contigPosInfo;
typedef std::vector<contigPosInfo> contigList;
typedef std::map<int , contigList> refContigInfo;

void loadRefContigInfo( const std::string & file , refContigInfo & data);

/*****************************************************************************
 *   file name      :
 *                      barcode_onref
 *   file format    :
 *                  1   0|123_456|...
 *
 * **************************************************************************/
// In barcodeList , means contig pos , 2 means barcode ids
typedef std::vector<int> barcodeList;
typedef std::map<int , barcodeList> refBarcodeInfo;

void loadRefBarcodeInfo(const std::string & file , refBarcodeInfo & data);

/*****************************************************************************
 *   file name      :
 *                      barcode_oncontig
 *   file format    :
 *                  1   0:123|456|...   1:1234|32423 ...
 *          contigId    pos:barcodeList
 *
 * **************************************************************************/
typedef std::tuple< int , barcodeList> barcodePosInfo;
typedef std::vector<barcodePosInfo> barcodePosList;
typedef std::map<int,barcodePosList> contigBarcodeInfo;
void generateConrigBarcodeInfo( const refBarcodeInfo & i_b ,
                                const refContigInfo & i_c,
                                contigBarcodeInfo & data );

void printContigBarcodeInfo( const contigBarcodeInfo & data , const std::string & file );

void loadContigBarcodeInfo( const std::string & file,  contigBarcodeInfo & data );



/*****************************************************************************
 *   file name      :
 *                      barcode_onbin
 *   file format    :
 *                  1:2 33:1
 *          contigId:bin    barcode:num ...
 *
 * **************************************************************************/

//typedef std::map<int , std::map<int,barcodeList> > binBarcodeInfo;
typedef std::map<int, std::map<int, std::map< int, int > > > binBarcodeInfo;

void generateBinBarcodeInfo(const contigBarcodeInfo & data  , int binSize, binBarcodeInfo & d);
void saveBinBarcodeInfo(const std::string & file ,const  binBarcodeInfo &data);

//void loadBinBarcodeInto

}//JOB01
}//BGIQD
#endif //__CONTIG_BARCODE_H__
