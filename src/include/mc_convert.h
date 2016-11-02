/** 
 * mc_convert.h
 * converter between Matlab and C/C++
 *
 */


#ifndef mc_convert_h__
#define mc_convert_h__



#include "types.h"
#include <iterator>
#include "mex.h"

/** 
 * Convert values in range [beg,end) to a new mxArray.
 * The Matlab will free returned mxArray automatically.
 * @param beg iterator pointing to the beginning
 * @param end iterator pointing to the end (past-end)
 * @return pointer to the created mxArray r. r is N-by-1 where N = std::diff(beg,end); 
 */
template<typename iter_t>
mxArray* values_to_new_mxArr(iter_t beg, iter_t end) {
    using namespace std;
    typedef typename iterator_traits<iter_t>::value_type T;
    typedef typename iterator_traits<iter_t>::pointer ptr;

    const int ndims = 1;
    mwSize dims[ndims];
    dims[0] = distance(beg, end);
    mxClassID id;
    id = cm_traits<T>::CID;

    mxArray* p = mxCreateNumericArray(1, dims, id, mxREAL);
    ptr pp = (ptr) mxGetData(p);
    copy(beg, end, pp);
    return p;
}

/** 
 * Copy values in mxArray mat(:) to the range beginning with beg.
 * The caller should ensure the range is large enough for mat.
 * @param mat pointer to the mxArray to be copied
 * @param beg iterator pointing to the beginning
 */
template<typename iter_t>
void mxArr_to_values(const mxArray* mat, iter_t beg) {
    using namespace std;

    void* pp = mxGetData(mat);
    mxClassID id = mxGetClassID(mat);
    mwSize n = mxGetNumberOfElements(mat);

    if (mxDOUBLE_CLASS == id) {
        typedef mc_traits<mxDOUBLE_CLASS>::CT T;
        T* pbeg = static_cast<T*> (pp);
        T* pend = pbeg + n;
        copy(pbeg, pend, beg);
    }

}

/** 
 * Convert scalar value *it to a new mxArray.
 * The Matlab will free returned mxArray automatically.
 * @param it iterator pointing to the scalar to be converted
 * @return pointer to the created mxArray r. r is 1-by-1
 */
template<typename iter_t>
mxArray* scalar_to_new_mxArr(iter_t it) {
    // TODO
    return 0;
}

/** 
 * Copy the first value in mxArray mat(:) to *it
 * @param mat pointer to the mxArray to be copied
 * @param it iterator to the copied scalar
 */
template<typename iter_t>
void mxArr_to_scalar(const mxArray* mat, iter_t it) {
    using namespace std;

    void* pp = mxGetData(mat);
    mxClassID id = mxGetClassID(mat);

    if (mxDOUBLE_CLASS == id) {
        typedef mc_traits<mxDOUBLE_CLASS>::CT T;
        T* p = static_cast<T*> (pp);
        *it = *p;
    }
}



#ifdef HAS_OPENCV

#include <cxtypes.h>

/** 
 * Convert values in IplImage *img to a new mxArray.
 * The Matlab will free returned mxArray automatically.
 * @param img pointer to the IplImage to be converted
 * @return pointer to the created mxArray r. r is the same dimensions and data type as *img. 
 */
mxArray* IplImage_to_new_mxArr(const IplImage* img);

/** 
 * Convert values in mxArray mat(:) to a new IplImage.
 * The caller is responsible to free returned pointer to IplImage.
 * @param mat pointer to the mxArray to be converted
 * @return pointer to the created IplImage r. r is the same dimensions and data type as *mat. 
 */
IplImage* mxArr_to_new_IplImage(const mxArray* mat);

/** 
 * Convert values in CvMat *mat to a new mxArray.
 * The Matlab will free returned mxArray automatically.
 * @param mat pointer to the CvMat to be converted
 * @return pointer to the created mxArray r. r is the same dimensions and data type as *mat. 
 */
mxArray* CvMat_to_new_mxArr(const CvMat* mat);
/** 
 * Convert values in mxArray mat(:) to a new CvMat.
 * The caller is responsible to free returned pointer to CvMat.
 * @param mat pointer to the mxArray to be converted
 * @return pointer to the created CvMat r. r is the same dimensions and data type as *mat. 
 */
CvMat* mxArr_to_new_CvMat(const mxArray* mat);

#endif // HAS_OPENCV





/*
// Add support to other data type conversion here...
#ifdef HAS_SOME_LIBARARY

#include "lib_header.h"

mxArray* datatype_to_new_mat (datatype* d);

datatype* mat_to_new_datatype (mxArray* mat);

#endif // HAS_SOME_LIBARARY
 */

/** 
 * Deprecated. Use values_to_new_mxArr instead.
 */
template<typename iter_t>
mxArray* values_to_new_mat(iter_t beg, iter_t end) {
    return values_to_new_mxArr(beg, end);
}

/** 
 * Deprecated. Use mxArr_to_values instead.
 */
template<typename iter_t>
void mat_to_values(const mxArray* mat, iter_t beg) {
    mxArr_to_values(mat, beg);
}

/** 
 * Deprecated. Use scalar_to_new_mxArr instead.
 */
template<typename iter_t>
mxArray* scalar_to_new_mat(iter_t it) {
    return scalar_to_new_mxArr(it);
}

/** 
 * Deprecated. Use mxArr_to_scalar instead.
 */
template<typename iter_t>
void mat_to_scalar(const mxArray* mat, iter_t it) {
    mxArr_to_scalar(mat, it);
}


#ifdef HAS_OPENCV

#include <cxtypes.h>

/** 
 * Deprecated. Use image_to_new_mxArr instead.
 */
mxArray* image_to_new_mat(const IplImage* img);

/** 
 *  Deprecated. Use mxArr_to_new_image instead.
 */
IplImage* mat_to_new_image(const mxArray* mat);

#endif // HAS_OPENCV
#endif // mc_convert_h__