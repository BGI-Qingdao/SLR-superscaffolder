#ifndef __SOAP2_KMER_H__
#define __SOAP2_KMER_H__
#include <stdint.h>
namespace BGIQD {
    namespace SOAP2 {
#if K127mer
        struct Kmer {
            uint64_t khh;
            uint64_t khl;
            uint64_t klh;
            uint64_t kll;
        };
#else
        struct Kmer {
            uint64_t kh;
            uint64_t kl;
        };
#endif
    }
}
#endif //__SOAP2_KMER_H__
