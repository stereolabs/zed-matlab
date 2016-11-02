/** 
 * types.h
 * types and traits
 */

#ifndef types_h__
#define types_h__

#include <mex.h>
// Uncomment below #define if your Matlab is prior to R2008a
//#ifndef mwSize
//#define mwSize size_t
//#endif

/*
 * Bi-directional mapping between C/C++ native types and Matlab class id
 */

// C/C++ native type -> matlab class id 

template<typename T>
struct cm_traits {
};

template<>
struct cm_traits<double> {
    static const mxClassID CID = mxDOUBLE_CLASS;
};

template<>
struct cm_traits<const double> {
    static const mxClassID CID = mxDOUBLE_CLASS;
};

template<>
struct cm_traits<float> {
    static const mxClassID CID = mxSINGLE_CLASS;
};

template<>
struct cm_traits<const float> {
    static const mxClassID CID = mxSINGLE_CLASS;
};

template<>
struct cm_traits<uint64_T> {
    static const mxClassID CID = mxUINT64_CLASS;
};

template<>
struct cm_traits<const uint64_T> {
    static const mxClassID CID = mxUINT64_CLASS;
};

template<>
struct cm_traits<uint32_T> {
    static const mxClassID CID = mxUINT32_CLASS;
};

template<>
struct cm_traits<const uint32_T> {
    static const mxClassID CID = mxUINT32_CLASS;
};

template<>
struct cm_traits<int32_T> {
    static const mxClassID CID = mxUINT32_CLASS;
};

template<>
struct cm_traits<const int32_T> {
    static const mxClassID CID = mxUINT32_CLASS;
};

template<>
struct cm_traits<uint16_T> {
    static const mxClassID CID = mxUINT16_CLASS;
};

template<>
struct cm_traits<const uint16_T> {
    static const mxClassID CID = mxUINT16_CLASS;
};

template<>
struct cm_traits<uint8_T> {
    static const mxClassID CID = mxUINT8_CLASS;
};

template<>
struct cm_traits<const uint8_T> {
    static const mxClassID CID = mxUINT8_CLASS;
};

template<>
struct cm_traits<wchar_t> {
    static const mxClassID CID = mxCHAR_CLASS;
};

template<>
struct cm_traits<const wchar_t> {
    static const mxClassID CID = mxCHAR_CLASS;
};

template<>
struct cm_traits<char> {
    static const mxClassID CID = mxINT8_CLASS;
};

template<>
struct cm_traits<const char> {
    static const mxClassID CID = mxINT8_CLASS;
};

// matlab class id -> C/C++ native type

template<mxClassID T>
struct mc_traits {
};

template<>
struct mc_traits<mxDOUBLE_CLASS> {
    typedef double CT;
};

template<>
struct mc_traits<mxSINGLE_CLASS> {
    typedef float CT;
};

template<>
struct mc_traits<mxUINT64_CLASS> {
    typedef uint64_T CT;
};

template<>
struct mc_traits<mxUINT32_CLASS> {
    typedef uint32_T CT;
};

template<>
struct mc_traits<mxINT32_CLASS> {
    typedef int32_T CT;
};

template<>
struct mc_traits<mxINT16_CLASS> {
    typedef int16_T CT;
};

template<>
struct mc_traits<mxUINT16_CLASS> {
    typedef uint16_T CT;
};

template<>
struct mc_traits<mxUINT8_CLASS> {
    typedef uint8_T CT;
};

template<>
struct mc_traits<mxINT8_CLASS> {
    typedef char CT;
};

template<>
struct mc_traits<mxLOGICAL_CLASS> {
    typedef uint8_T CT;
};

template<>
struct mc_traits<mxCHAR_CLASS> {
    typedef wchar_t CT;
};


//#ifdef HAS_OPENCV

//#include <cxtypes.h>
/*
 * Bi-directional mapping between openCV and Matlab class id
 */

// openCV -> matlab class id 

template<int>
struct cvm_traits {
};

template<>
struct cvm_traits<IPL_DEPTH_64F> {
    static const mxClassID CID = mxDOUBLE_CLASS;
};

template<>
struct cvm_traits<IPL_DEPTH_32F> {
    static const mxClassID CID = mxSINGLE_CLASS;
};

template<>
struct cvm_traits<IPL_DEPTH_8U> {
    static const mxClassID CID = mxUINT8_CLASS;
};

template<>
struct cvm_traits<CV_64FC1> {
    static const mxClassID CID = mxDOUBLE_CLASS;
};

template<>
struct cvm_traits<CV_32FC1> {
    static const mxClassID CID = mxSINGLE_CLASS;
};

template<>
struct cvm_traits<CV_32SC1> {
    static const mxClassID CID = mxINT32_CLASS;
};

template<>
struct cvm_traits<CV_16SC1> {
    static const mxClassID CID = mxINT16_CLASS;
};

template<>
struct cvm_traits<CV_16UC1> {
    static const mxClassID CID = mxUINT16_CLASS;
};

template<>
struct cvm_traits<CV_8SC1> {
    static const mxClassID CID = mxINT8_CLASS;
};

template<>
struct cvm_traits<CV_8UC1> {
    static const mxClassID CID = mxUINT8_CLASS;
};


// matlab class id-> openCV 

template<mxClassID>
struct mcv_traits {
};

template<>
struct mcv_traits<mxDOUBLE_CLASS> {
    static const int CV_DEPTH = IPL_DEPTH_64F;
    static const int CV_TYPE = CV_64FC1;
};

template<>
struct mcv_traits<mxSINGLE_CLASS> {
    static const int CV_DEPTH = IPL_DEPTH_32F;
    static const int CV_TYPE = CV_32FC1;
};

template<>
struct mcv_traits<mxINT32_CLASS> {
    static const int CV_DEPTH = IPL_DEPTH_32S;
    static const int CV_TYPE = CV_32SC1;
};

template<>
struct mcv_traits<mxINT16_CLASS> {
    static const int CV_DEPTH = IPL_DEPTH_16S;
    static const int CV_TYPE = CV_16SC1;
};

template<>
struct mcv_traits<mxUINT16_CLASS> {
    static const int CV_DEPTH = IPL_DEPTH_16U;
    static const int CV_TYPE = CV_16UC1;
};

template<>
struct mcv_traits<mxINT8_CLASS> {
    static const int CV_DEPTH = IPL_DEPTH_8S;
    static const int CV_TYPE = CV_8SC1;
};

template<>
struct mcv_traits<mxUINT8_CLASS> {
    static const int CV_DEPTH = IPL_DEPTH_8U;
    static const int CV_TYPE = CV_8UC1;
};

template<>
struct mcv_traits<mxLOGICAL_CLASS> {
    static const int CV_DEPTH = IPL_DEPTH_8U;
    static const int CV_TYPE = CV_8UC1;
};

//#endif // HAS_OPENCV

#endif // types_h__