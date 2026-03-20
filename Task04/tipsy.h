#include "boost/endian/arithmetic.hpp"
#include "boost/endian/conversion.hpp"

namespace tipsy {
    using namespace boost::endian;

    // Header of a Tipsy file
    struct header {
        double dTime;
        std::uint32_t nBodies;
        std::uint32_t nDim;
        std::uint32_t nSph;
        std::uint32_t nDark;
        std::uint32_t nStar;
        std::uint32_t nPad;
    };

    inline void swap(header &hdr) {
        using boost::endian::big_to_native_inplace;
        
        big_to_native_inplace(hdr.dTime);
        big_to_native_inplace(hdr.nBodies);
        big_to_native_inplace(hdr.nDim);
        big_to_native_inplace(hdr.nSph);
        big_to_native_inplace(hdr.nDark); 
        big_to_native_inplace(hdr.nStar);
        big_to_native_inplace(hdr.nPad);
    }

    // Dark matter particle
    struct dark {
        float mass;
        float pos[3];
        float vel[3];
        float eps;
        float phi;
    };

    inline void swap(dark &d) {
        using boost::endian::big_to_native_inplace;
        // Swap all of the fields in "d"
        big_to_native_inplace(d.mass);
        for(int i=0; i<3; ++i) big_to_native_inplace(d.pos[i]);
        for(int i=0; i<3; ++i) big_to_native_inplace(d.vel[i]);
        big_to_native_inplace(d.eps);
        big_to_native_inplace(d.phi);

    }

}
