/**
 * iter.h
 * iterator for data conversion
 */

#ifndef iter_h__
#define iter_h__

#include <iterator>

/** 
 * Image Storage type.
 * Generally Maltab is column wise, while C++ row wise
 */
enum eStorType {
    eColWise, eRowWise
};

/// Pixel iterator for 2D image, i.e. mono-color image

template<typename data_t, eStorType st = eRowWise>
class pix_iterator_2d : public std::iterator<std::forward_iterator_tag, data_t> {
public:

    /** 
     * Constructor
     * @param beg Beginning of the image
     * @param width width of image in pixels
     * @param height Height of image in pixels
     * @param pitch Width step of aligned image in BYTES
     */
    pix_iterator_2d(data_t* beg, int width, int height, int pitch)
    : pBeg_(beg), pCurPix_(beg), pCurRowBeg_(beg),
    curCol_(0), curRow_(0), width_(width), height_(height) {
        size_t s = sizeof (data_t);
        if (st == eRowWise) pitch_ = ((unsigned) pitch < width * s) ? (width * s) : (pitch);
        else if (st == eColWise) pitch_ = (pitch < height * s) ? (height * s) : (pitch);
#if 0 
        mexEvalString("pcc = [];");
        mexEvalString("pcr = [];");
#endif // _DEBUG
    }

    /// Constructor

    pix_iterator_2d(data_t* beg, int width, int height)
    : pBeg_(beg), pCurPix_(beg), pCurRowBeg_(beg),
    curCol_(0), curRow_(0), width_(width), height_(height) {
        if (st == eRowWise) pitch_ = width * sizeof (data_t);
        else if (st == eColWise) pitch_ = height * sizeof (data_t);
#if 0 
        mexEvalString("pcc = [];");
        mexEvalString("pcr = [];");
#endif // _DEBUG
    }

    /// reset

    void reset(data_t* beg, int width, int height, int pitch) {
        pBeg_ = beg;
        pCurPix_ = beg;
        pCurRowBeg_ = beg;
        curCol_ = 0;
        curRow_ = 0;
        width_ = width;
        height_ = height;

        size_t s = sizeof (data_t);
        if (st == eRowWise) pitch_ = (pitch < width * s) ? (width * s) : (pitch);
        else if (st == eColWise) pitch_ = (pitch < height * s) ? (height * s) : (pitch);

    }

    /// reset 

    void reset(data_t* beg, int width, int height) {
        pBeg_ = beg;
        pCurPix_ = beg;
        pCurRowBeg_ = beg;
        curCol_ = 0;
        curRow_ = 0;
        width_ = width;
        height_ = height;

        if (st == eRowWise) pitch_ = width * sizeof (data_t);
        else if (st == eColWise) pitch_ = height * sizeof (data_t);
#if 0 
        mexEvalString("pcc = [];");
        mexEvalString("pcr = [];");
#endif // _DEBUG
    }

    /// Destructor

    ~pix_iterator_2d() {
    }
    /// operator !=

    bool operator!=(const pix_iterator_2d<data_t, st>& rhs) {
        return pCurPix_ != rhs.pCurPix_;
    }

    /// pre-increment

    pix_iterator_2d<data_t, st>& operator++() {
        if (curCol_ < width_ - 1) {
            ++curCol_;
            if (st == eRowWise) ++pCurPix_;
            else if (st == eColWise) {
                char* tmp = (char*) pCurPix_;
                tmp += pitch_;
                pCurPix_ = (data_t*) tmp;
            }

#if 0
            mxArray* pcc = MexUtil::createArrayFromScalar(curCol_);
            mexPutVariable("base", "tmp", pcc);
            mexEvalString("pcc(end+1) = tmp;");
#endif
        } else { // next row
            curCol_ = 0;
            ++curRow_;
#if 0
            mexEvalString("pcc = 0;");

            mxArray* pcr = MexUtil::createArrayFromScalar(curRow_);
            mexPutVariable("base", "tmp", pcr);
            mexEvalString("pcr(end+1) = tmp;");
#endif // _DEBUG
            if (st == eRowWise) {
                char* tmp = (char*) pCurRowBeg_;
                tmp += pitch_;
                pCurRowBeg_ = (data_t*) tmp;
                pCurPix_ = pCurRowBeg_;
            } else if (st == eColWise) {
                ++pCurRowBeg_;
                if (curRow_ < height_) pCurPix_ = pCurRowBeg_;
                else pCurPix_ = (data_t*) ((char*) pBeg_ + pitch_ * width_);
            }

        }
        return *this;
    }

    /// dereference

    data_t& operator*() {
        return *pCurPix_;
    }
    /// pass-end. STL semantic

    pix_iterator_2d<data_t, st>& end() {
        curCol_ = 0;
        if (st == eRowWise)
            pCurRowBeg_ = (data_t*) ((char*) pBeg_ + pitch_ * height_);
        else if (st == eColWise)
            pCurRowBeg_ = (data_t*) ((char*) pBeg_ + pitch_ * width_);
        pCurPix_ = pCurRowBeg_;
        return *this;
    }
protected:
private:
    data_t *pCurPix_, *pCurRowBeg_; // ��ǰ���أ��е�ָ��
    data_t* pBeg_;
    int curCol_, curRow_;
    int width_, height_; // in pixels
    int pitch_; // in bytes
};



// mxArray_iter_3d
// iterate among pages
// see MATLAB on-line help:
//   MATLAB->Programming->Multidemensional Arrays->Overview
// for the details about the concept of "row", "column" and "page"

template<typename data_t>
class mxArray_iter_3d
: public std::iterator<std::forward_iterator_tag, data_t> {
public:

    mxArray_iter_3d(data_t* start,
            unsigned int width, unsigned int height, unsigned int npage)
    : cur_page_start_(start), start_(start), elemcount_(0),
    width_(width), height_(height), pagesize_(width*height), npage_(npage),
    it_(start, width, height) {
    }

    //

    mxArray_iter_3d<data_t>& operator++() {
        ++elemcount_;
        ++it_;
        if (!(elemcount_ < pagesize_)) { // next page
            cur_page_start_ += pagesize_;
            it_.reset(cur_page_start_, width_, height_);
            elemcount_ = 0;
        }
        return *this;
    }

    //

    bool operator==(const mxArray_iter_3d<data_t>& rhs) {
        return it_ == rhs.it_;
    }

    //

    bool operator!=(const mxArray_iter_3d<data_t>& rhs) {
        return it_ != rhs.it_;
    }

    // 

    data_t& operator*() {
        return *it_;
    }

    //

    mxArray_iter_3d<data_t>& end() {
        cur_page_start_ = start_ + (npage_ - 1) * pagesize_;
        it_.reset(cur_page_start_, width_, height_);
        it_.end();
        return *this;
    }
protected:
private:
    const unsigned int width_, height_, pagesize_, npage_;
    unsigned int elemcount_;
    data_t *cur_page_start_, *start_;
    pix_iterator_2d<data_t, eColWise> it_;
};

/**
 * Pixel iterator for 3-channel interleaved image
 * Data is stored as b0,g0,r0,...bn,gn,rn and
 * fetched as r0,...,rn,g0,...,gn,b0,...,bn
 */
template<typename data_t>
class pix_iter_rgb
: public std::iterator<std::forward_iterator_tag, data_t> {
public:
    /// Constructor

    pix_iter_rgb(data_t* beg, int width, int height, int step)
    : pBeg_(beg),
    w_(0), h_(0), c_(0), width_(width), height_(height) {
        this->reset_ptr_channel();
        size_t s = sizeof (data_t);
        step_ = (step < width * s) ? (width * s) : (step);
    }
    /// reset

    void reset(data_t* beg, int width, int height, int step) {
        w_ = h_ = c_ = 0;
        width_ = width;
        height_ = height;
        pBeg_ = beg;
        this->reset_ptr_channel();
        size_t s = sizeof (data_t);
        step_ = (step < width * s) ? (width * s) : (step);
    }
    /// Destructor

    ~pix_iter_rgb() {
    }

    /// operator 1=

    bool operator!=(const pix_iter_rgb<data_t>& rhs) {
        return p_ != rhs.p_;
    }
    /// pre-increment

    pix_iter_rgb<data_t>& operator++() {
        if (w_ < width_ - 1) { // next pixel
            ++w_;
            p_ += 3;
        } else {
            if (h_ < height_ - 1) { // next row
                w_ = 0;
                ++h_;
                this->reset_ptr_next_row();
            } else {
                if (c_ < 2) { // next channel
                    ++c_;
                    w_ = h_ = 0;
                    this->reset_ptr_channel();
                } else { // pass-end 
                    this->end();
                }
            }
        }
        return *this;
    }
    /// dereference

    data_t& operator*() {
        return *p_;
    }
    /// pass-end iterator. STL style.

    pix_iter_rgb<data_t>& end() {
        w_ = 0;
        h_ = height_;
        p_ = pRow_ = (data_t*) ((char*) pBeg_ + step_ * height_);
        return *this;
    }
protected:
private:
    data_t *p_, *pRow_;
    data_t *pBeg_; // image origin
    int w_, h_, c_; // current width, height and channel
    int width_, height_; // in pixels
    int step_; // in bytes
    // reset p_ and pRow_ to the origin of current channel

    void reset_ptr_channel() {
        pRow_ = pBeg_ + (2 - c_);
        p_ = pRow_;
    }
    // reset p_ and pRow_ to the beginning of next row

    void reset_ptr_next_row() {
        // NOTE: step in BYTES!!!
        char *tmp = (char*) pRow_;
        tmp += step_;
        pRow_ = (data_t*) tmp;
        p_ = pRow_;
    }
};
#endif // iter_h__