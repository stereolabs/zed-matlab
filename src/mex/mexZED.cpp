
/*********************************
 **     Using ZED with Matlab   **
 *********************************/

// MEX header
#include <mex.h>
#include "matrix.h"

// OpenCV
#include <opencv2/opencv.hpp>

// Matrix format conversion by Sun Peng : http://www.mathworks.com/matlabcentral/fileexchange/20927-c-c++-and-matlab-types-convertor
#include "mc_convert.h"
#include "iter.h"
#include "types.h"

// system header
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// ZED
#include <zed/Camera.hpp>

//CUDA includes
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <npp.h>

#ifdef _DEBUG
#error Please select Release mode for compilation
#endif

// global var.
static sl::zed::Camera *zedCam = NULL;
sl::zed::Mat cloudCPU;

// C++ struct

struct Intra {
    float fx, fy, cx, cy;
};

struct StereoParams {
    float baseline;
    Intra left;
    Intra right;
};

const char *fieldsIntra[] = {"fx", "fy", "cx", "cy"};
const char *fieldsStruct[] = {"baseline", "left", "right"};

// Interop. OpenCV-Matlab matrix type (float)

template<class InputIterator, class OutputIterator>
OutputIterator fctCopy(InputIterator first, InputIterator last, OutputIterator result) {
    while (first != last) {
        *result = *first;
        ++result;
        ++first;
    }
    return result;
}

template<int TYPE>
mxArray* floatmat_to_mat(cv::Mat &matrix) {
    void* pBeg = matrix.data;
    int pitch = matrix.step;
    cv::Size size = matrix.size();
    const mxClassID cid = cvm_traits<TYPE>::CID;
    mxArray* pArrOut = mxCreateNumericMatrix(size.height, size.width, cid, mxREAL);
    void* pBegOut = mxGetData(pArrOut);

    typedef typename mc_traits<cid>::CT T;
    pix_iterator_2d<T, eRowWise> it_src1(static_cast<T*> (pBeg), size.width, size.height, pitch);
    pix_iterator_2d<T, eRowWise> it_src2(static_cast<T*> (pBeg), size.width, size.height, pitch);
    it_src2.end();
    pix_iterator_2d<T, eColWise> it_dest(static_cast<T*> (pBegOut), size.width, size.height);

    fctCopy(it_src1, it_src2, it_dest);

    return pArrOut;
}

// Interop. OpenCV-Matlab matrix type (rgb)

template<int DEPTH>
mxArray* rgbimage_to_mat(cv::Mat &image) {
    const int ndim = image.channels();
    cv::Size size = image.size();
    mwSize dims[3];
    dims[0] = size.height;
    dims[1] = size.width;
    dims[2] = ndim;
    const mxClassID cid = cvm_traits<DEPTH>::CID;
    mxArray* pArrOut = mxCreateNumericArray(ndim, dims, cid, mxREAL);

    void* pBeg = image.data;
    int pitch = image.step;
    void* pBegOut = mxGetData(pArrOut);

    typedef typename mc_traits<cid>::CT T;
    pix_iter_rgb<T> it_src1(static_cast<T*> (pBeg), size.width, size.height, pitch);
    pix_iter_rgb<T> it_src2(static_cast<T*> (pBeg), size.width, size.height, pitch);
    it_src2.end();
    mxArray_iter_3d<T> it_dest(static_cast<T*> (pBegOut), size.width, size.height, ndim);
    fctCopy(it_src1, it_src2, it_dest);
    return pArrOut;
}

// Interop. OpenCV-Matlab matrix type

mxArray* CvMat_to_new_mxArr(cv::Mat &matrix) {
    const int TYPE = matrix.type();
    if (CV_32FC1 == TYPE)
        return floatmat_to_mat<CV_32FC1>(matrix);
    else if (CV_8UC3 == TYPE)
        return rgbimage_to_mat<IPL_DEPTH_8U>(matrix);

    return mxCreateDoubleMatrix(0, 0, mxREAL);
}

sl::zed::MODE exctractMode(char *_str) {
    sl::zed::MODE computeMode = sl::zed::MODE::PERFORMANCE;
    if (!strcmp(_str, "PERFORMANCE"))
        computeMode = sl::zed::MODE::PERFORMANCE;
    else if (!strcmp(_str, "MEDIUM"))
        computeMode = sl::zed::MODE::MEDIUM;
    else if (!strcmp(_str, "QUALITY"))
        computeMode = sl::zed::MODE::QUALITY;
    else
        mexPrintf("unknown mode, 'PERFORMANCE' used\n");
    return computeMode;
}

sl::zed::UNIT exctractUnit(char *_str) {
    sl::zed::UNIT unit = sl::zed::UNIT::METER;
    if (!strcmp(_str, "MILLIMETER"))
        unit = sl::zed::UNIT::MILLIMETER;
    else if (!strcmp(_str, "METER"))
        unit = sl::zed::UNIT::METER;
    else if (!strcmp(_str, "INCH"))
        unit = sl::zed::UNIT::INCH;
    else if (!strcmp(_str, "FOOT"))
        unit = sl::zed::UNIT::FOOT;
    else
        mexPrintf("unknown UNIT -> 'METER' used\n");
    return unit;
}

/* MEX entry function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    // each sub function is referred by a string 'command'
    char command[128];
    mxGetString(prhs[0], command, 128);

    // creation of the ZED, in live or playback mode
    if (!strcmp(command, "create")) {
        // is the first param. is a character this represent a string containing a SVO file
        if (mxIsChar(prhs[1])) {
            char svo_path[128];
            mxGetString(prhs[1], svo_path, 128);
            mexPrintf("opening svo : %s\n", svo_path);
            zedCam = new sl::zed::Camera(svo_path);
        } else { // else it should be a resolution for live
            double *ptr_ = mxGetPr(prhs[1]);
            int reso = ptr_[0];
            mexPrintf("capture live at resolution %dp\n", reso);
            sl::zed::ZEDResolution_mode resolution;

            switch (reso) {
                case 376:
                    resolution = sl::zed::VGA;
                    break;
                case 720:
                    resolution = sl::zed::HD720;
                    break;
                case 1080:
                    resolution = sl::zed::HD1080;
                    break;
                case 1242:
                    resolution = sl::zed::HD2K;
                    break;
                default:
                    mexPrintf("unknown resolution, 1080p used\n");
                    resolution = sl::zed::HD1080;
                    break;
            }
            int fps = 0;
            if (nrhs == 3) {
                ptr_ = mxGetPr(prhs[2]);
                fps = ptr_[0];
                mexPrintf("capture live at %dFPS\n", fps);
            }
            zedCam = new sl::zed::Camera(resolution, fps);
        }
    }

    //initialisation
    if (!strcmp(command, "init")) {
        if (zedCam) {
            sl::zed::InitParams params;
            params.mode = sl::zed::MODE::MEDIUM;
            params.unit = sl::zed::UNIT::METER;
            params.verbose = true;

            // check if we have ONE argument
            if (nrhs == 2) {
                // check if we have a parameter structure
                if (mxIsStruct(prhs[1])) {
                    // for all fields of parameter structure overwrite parameters
                    for (int32_t i = 0; i < mxGetNumberOfFields(prhs[1]); i++) {
                        const char *field_name = mxGetFieldNameByNumber(prhs[1], i);
                        //mexPrintf("field %d : %s ", i, field_name);
                        mxArray *field_val = mxGetFieldByNumber(prhs[1], 0, i);
                        int val = *((double*) mxGetPr(field_val));
                        //mexPrintf(" val %d  \n", val);
                        if (!strcmp(field_name, "mode")) params.mode = static_cast<sl::zed::MODE> (val);
                        if (!strcmp(field_name, "unit")) params.unit = static_cast<sl::zed::UNIT> (val);
                        if (!strcmp(field_name, "coordinate")) params.coordinate = static_cast<sl::zed::COORDINATE_SYSTEM> (val);
                        if (!strcmp(field_name, "verbose")) params.verbose = val;
                        if (!strcmp(field_name, "device")) params.device = val;
                        if (!strcmp(field_name, "minimumDistance")) params.minimumDistance = val;
                        if (!strcmp(field_name, "disableSelfCalib")) params.disableSelfCalib = val;
                        if (!strcmp(field_name, "vflip")) params.vflip = val;
                    }
                }
            }
            sl::zed::ERRCODE err = zedCam->init(params);
            // we return the string associated with the error
            plhs[0] = mxCreateString(sl::zed::errcode2str(err).c_str());
            cloudCPU.allocate_cpu(zedCam->getImageSize().width, zedCam->getImageSize().height, 4, sl::zed::DATA_TYPE::FLOAT);
        } else
            plhs[0] = mxCreateString("error");
    }

    // get image size function
    if (!strcmp(command, "getImageSize")) {
        if (zedCam) {
            double ptr_size[2];
            ptr_size[0] = zedCam->getImageSize().width;
            ptr_size[1] = zedCam->getImageSize().height;
            plhs[0] = mxCreateDoubleMatrix(2, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), ptr_size, 2 * sizeof (double));
        } else
            mexErrMsgTxt("ZED is not initialized");
    }

    // grab function, the SENSING MODE can be given "full"-"raw", by default "full" is used
    if (!strcmp(command, "grab")) {
        if (zedCam) {
            if (nrhs == 2) {
                char mode[64];
                mxGetString(prhs[1], mode, 64);
                if (!strcmp(mode, "FILL"))
                    zedCam->grab(sl::zed::SENSING_MODE::FILL);
                else if (!strcmp(mode, "STANDARD"))
                    zedCam->grab(sl::zed::SENSING_MODE::STANDARD);
                else
                    mexErrMsgTxt("wrong sensing mode");
            } else {
                zedCam->grab();
            }
        } else
            mexErrMsgTxt("ZED is not initialized");
    }

    // retrieve image function, the side must be specified : "left" - "right"
    if (!strcmp(command, "retrieveImage")) {
        if (zedCam) {
            char side[64];
            mxGetString(prhs[1], side, 64);

            cv::Mat image_rgb;
            if (!strcmp(side, "left"))
                cv::cvtColor(slMat2cvMat(zedCam->retrieveImage(sl::zed::SIDE::LEFT)), image_rgb, CV_RGBA2RGB);
            else if (!strcmp(side, "right"))
                cv::cvtColor(slMat2cvMat(zedCam->retrieveImage(sl::zed::SIDE::RIGHT)), image_rgb, CV_RGBA2RGB);
            else
                mexErrMsgTxt("wrong image side");

            plhs[0] = CvMat_to_new_mxArr(image_rgb);
        } else
            mexErrMsgTxt("ZED is not initialized");
    }

    // retrieve normalized measure image, the measure type must be specified : "disparity" - "depth" - "confidence"
    if (!strcmp(command, "normalizeMeasure")) {
        if (zedCam) {
            char measure[64];
            mxGetString(prhs[1], measure, 64);

            cv::Mat image_rgb;
            if (!strcmp(measure, "disparity"))
                cv::cvtColor(slMat2cvMat(zedCam->normalizeMeasure(sl::zed::MEASURE::DISPARITY)), image_rgb, CV_RGBA2RGB);
            else if (!strcmp(measure, "depth"))
                cv::cvtColor(slMat2cvMat(zedCam->normalizeMeasure(sl::zed::MEASURE::DEPTH)), image_rgb, CV_RGBA2RGB);
            else if (!strcmp(measure, "confidence"))
                cv::cvtColor(slMat2cvMat(zedCam->normalizeMeasure(sl::zed::MEASURE::CONFIDENCE)), image_rgb, CV_RGBA2RGB);
            else
                mexErrMsgTxt("wrong measure");
            plhs[0] = CvMat_to_new_mxArr(image_rgb);

        } else
            mexErrMsgTxt("ZED is not initialized");
    }

    // retrieve measure matrix, the measure type must be specified : "disparity" - "depth" - "confidence"
    if (!strcmp(command, "retrieveMeasure")) {
        if (zedCam) {
            char measure[64];
            mxGetString(prhs[1], measure, 64);

            cv::Mat measure_mat;
            if (!strcmp(measure, "disparity"))
                measure_mat = slMat2cvMat(zedCam->retrieveMeasure(sl::zed::MEASURE::DISPARITY));
            else if (!strcmp(measure, "depth"))
                measure_mat = slMat2cvMat(zedCam->retrieveMeasure(sl::zed::MEASURE::DEPTH));
            else if (!strcmp(measure, "confidence"))
                measure_mat = slMat2cvMat(zedCam->retrieveMeasure(sl::zed::MEASURE::CONFIDENCE));
            else if (!strcmp(measure, "XYZ")) {
                sl::zed::Mat cloud = zedCam->retrieveMeasure_gpu(sl::zed::MEASURE::XYZ);
                cudaError err = cudaMemcpy2D((uchar*) cloudCPU.data, cloudCPU.step, (Npp32f*) cloud.data, cloud.step, cloud.getWidthByte(), cloud.height, cudaMemcpyDeviceToHost);
                measure_mat = slMat2cvMat(cloudCPU);
            } else
                mexErrMsgTxt("wrong measure");

            if (strcmp(measure, "XYZ"))
                plhs[0] = CvMat_to_new_mxArr(measure_mat);
            else {
                std::vector<cv::Mat> mat_v;
                cv::split(measure_mat, mat_v);
                plhs[0] = CvMat_to_new_mxArr(mat_v[0]);
                plhs[1] = CvMat_to_new_mxArr(mat_v[1]);
                plhs[2] = CvMat_to_new_mxArr(mat_v[2]);
            }
        } else
            mexErrMsgTxt("ZED is not initialized");
    }

    // get number of frames, return -1 if the ZED is used in live mode
    if (!strcmp(command, "getSVONumberOfFrames")) {
        if (zedCam) {
            double nbFrame = zedCam->getSVONumberOfFrames();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &nbFrame, 1 * sizeof (double));
        } else
            mexErrMsgTxt("ZED is not initialized");
    }

    // set the current frame position
    if (!strcmp(command, "setSVOPosition")) {
        if (zedCam) {
            double *ptr_ = mxGetPr(prhs[1]);
            int frame = ptr_[0];
            zedCam->setSVOPosition(frame);
        } else
            mexErrMsgTxt("ZED is not initialized");
    }

    // set the confidence threshold
    if (!strcmp(command, "setConfidenceThreshold")) {
        if (zedCam) {
            double *ptr_ = mxGetPr(prhs[1]);
            int confidence = ptr_[0];
            zedCam->setConfidenceThreshold(confidence);
        } else
            mexErrMsgTxt("ZED is not initialized");
    }

    // set the depth clamp value
    if (!strcmp(command, "setDepthClampValue")) {
        if (zedCam) {
            double *ptr_ = mxGetPr(prhs[1]);
            int depthClamp = ptr_[0];
            zedCam->setDepthClampValue(depthClamp);
        } else
            mexErrMsgTxt("ZED is not initialized");
    }

    // get the projection matrix
    if (!strcmp(command, "getCameraParameters")) {
        if (zedCam) {
            mxArray *p;
            plhs[0] = mxCreateStructMatrix(1, 1, 3, fieldsStruct);
            // Left camera parameters
            p = mxCreateStructMatrix(1, 1, 4, fieldsIntra);
            mxSetField(p, 0, "fx", mxCreateDoubleScalar(zedCam->getParameters()->LeftCam.fx));
            mxSetField(p, 0, "fy", mxCreateDoubleScalar(zedCam->getParameters()->LeftCam.fy));
            mxSetField(p, 0, "cx", mxCreateDoubleScalar(zedCam->getParameters()->LeftCam.cx));
            mxSetField(p, 0, "cy", mxCreateDoubleScalar(zedCam->getParameters()->LeftCam.cy));
            mxSetField(plhs[0], 0, "left", p);

            // Right camera parameters
            p = mxCreateStructMatrix(1, 1, 4, fieldsIntra);
            mxSetField(p, 0, "fx", mxCreateDoubleScalar(zedCam->getParameters()->RightCam.fx));
            mxSetField(p, 0, "fy", mxCreateDoubleScalar(zedCam->getParameters()->RightCam.fy));
            mxSetField(p, 0, "cx", mxCreateDoubleScalar(zedCam->getParameters()->RightCam.cx));
            mxSetField(p, 0, "cy", mxCreateDoubleScalar(zedCam->getParameters()->RightCam.cy));
            mxSetField(plhs[0], 0, "right", p);

            // baseline
            mxSetField(plhs[0], 0, "baseline", mxCreateDoubleScalar(zedCam->getParameters()->baseline));
        } else
            mexErrMsgTxt("ZED is not initialized");
    }

    // set the depth clamp value
    if (!strcmp(command, "enableTracking")) {
        if (zedCam) {
            Eigen::Matrix4f pathInit;
            pathInit.setIdentity(4, 4);
            zedCam->enableTracking(pathInit);
        } else
            mexErrMsgTxt("ZED is not initialized");
    }


    // set the depth clamp value
    if (!strcmp(command, "enableTracking")) {
        if (zedCam) {
            Eigen::Matrix4f pathInit;
            pathInit.setIdentity(4, 4);
            zedCam->enableTracking(pathInit);
        } else
            mexErrMsgTxt("ZED is not initialized");
    }

    if (!strcmp(command, "getPosition")) {
        if (zedCam) {
            cv::Mat position(4, 4, CV_32FC1);
            Eigen::Matrix4f path;
            zedCam->getPosition(path);

            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    position.at<float>(i, j) = path(i, j);
                }
            }

            plhs[0] = CvMat_to_new_mxArr(position);
        } else
            mexErrMsgTxt("ZED is not initialized");
    }

    // get current timestamp, return the current timestamp of the camera.
    if (!strcmp(command, "getCurrentTimestamp")) {
        if (zedCam) {
            double timestamp = zedCam->getCurrentTimestamp();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &timestamp, 1 * sizeof (double));
        } else
            mexErrMsgTxt("ZED is not initialized");
    }

    // get camera timestamp, return the timestamp in ns when the frames has been extracted from the camera. To be called after a grab.
    if (!strcmp(command, "getCameraTimestamp")) {
        if (zedCam) {
            double timestamp = zedCam->getCameraTimestamp();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &timestamp, 1 * sizeof (double));
        } else
            mexErrMsgTxt("ZED is not initialized");
    }

    // set the depth clamp value
    if (!strcmp(command, "delete")) {
        if (zedCam) {
            cloudCPU.deallocate();
            delete zedCam;
            void mexUnlock(void);
        } else
            mexErrMsgTxt("ZED is not initialized");
    }

}
