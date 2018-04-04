#ifndef __SOAP2_KMER_H__
#define __SOAP2_KMER_H__
#include <stdint.h>
namespace BGIQD {
    namespace SOAP2 {
#if K127mer
        typedef uint64_t Kmer[4];
#else
        typedef uint64_t Kmer[2];
#endif

    }
}
#endif //__SOAP2_KMER_H__
