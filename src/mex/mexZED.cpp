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
#include <sl/Camera.hpp>

//CUDA includes
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <npp.h>

#ifdef _DEBUG
#error Please select Release mode for compilation
#endif

// global var.
static sl::Camera *zedCam = NULL;
sl::Mat cloudCPU;

// C++ struct

struct Intra {
    float fx, fy, cx, cy;
};

struct StereoParams {
    float baseline;
    Intra left;
    Intra right;
};

const char *fieldsIntra[] = {"cx", "cy", "disto", "d_fov", "fx" , "fy" , "h_fov" , "width" , "height" , "v_fov"};
const char *fieldsStruct[] = {"serial_number", "firmware_version", "baseline", "convergence", "r_x", "r_y", "r_z", "t_x", "t_y", "t_z", "left_cam", "right_cam"};

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

sl::DEPTH_MODE exctractMode(char *_str) {
    sl::DEPTH_MODE computeMode = sl::DEPTH_MODE_PERFORMANCE;
    if (!strcmp(_str, "PERFORMANCE"))
        computeMode = sl::DEPTH_MODE_PERFORMANCE;
    else if (!strcmp(_str, "MEDIUM"))
        computeMode = sl::DEPTH_MODE_MEDIUM;
    else if (!strcmp(_str, "QUALITY"))
        computeMode = sl::DEPTH_MODE_QUALITY;
    else
        mexPrintf("unknown mode, 'PERFORMANCE' used\n");
    return computeMode;
}

sl::UNIT exctractUnit(char *_str) {
    sl::UNIT unit = sl::UNIT_METER;
    if (!strcmp(_str, "MILLIMETER"))
        unit = sl::UNIT_MILLIMETER;
    else if (!strcmp(_str, "METER"))
        unit = sl::UNIT_METER;
    else if (!strcmp(_str, "INCH"))
        unit = sl::UNIT_INCH;
    else if (!strcmp(_str, "FOOT"))
        unit = sl::UNIT_FOOT;
    else
        mexPrintf("unknown UNIT -> 'METER' used\n");
    return unit;
}

void deleteOnFail() {
    if (zedCam) {
        delete zedCam;
        zedCam = NULL;
    }
}

bool checkZED() {
    if (zedCam)
        return true;
    else {
        mexErrMsgTxt("ZED is not initialized");
        deleteOnFail();
    }
    return false;
}

bool checkParams(int params, int required) {
    if (params - 1 != required) {
        std::string error = "Invalid parameter number, " + std::to_string(required) + " required, " + std::to_string(params - 1) + " given.";
        mexErrMsgTxt(error.c_str());
        deleteOnFail();
        return false;
    }
    return true;
}

bool checkParams(int params, int required_1, int required_2) {
    if ((params - 1 != required_1) && (params - 1 != required_2)) {
        std::string error = "Invalid parameter number, " + std::to_string(required_1) + " or " + std::to_string(required_2) + " required, " + std::to_string(params - 1) + " given.";
        mexErrMsgTxt(error.c_str());
        deleteOnFail();
        return false;
    }
    return true;
}

void notAvailable() {
    mexErrMsgTxt("This function isn't yet implemented on the Matlab plugin.");
    deleteOnFail();
}

/* MEX entry function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    // each sub function is referred by a string 'command'
    char command[128];
    mxGetString(prhs[0], command, 128);

    if (!strcmp(command, "create")) {
        if (checkParams(nrhs, 0))
            zedCam = new sl::Camera();
    }

    else if (!strcmp(command, "open")) {
        if (checkZED()) {
            sl::ERROR_CODE err;
            sl::InitParameters initParams;
            initParams.depth_mode = sl::DEPTH_MODE_MEDIUM;
            initParams.coordinate_units = sl::UNIT_METER;
            initParams.sdk_verbose = true;

            // check if we have ONE argument
            if (checkParams(nrhs, 1)) {
                // check if we have a parameter structure
                if (mxIsStruct(prhs[1])) {
                    // for all fields of parameter structure overwrite parameters
                    for (int32_t i = 0; i < mxGetNumberOfFields(prhs[1]); i++) {
                        const char *field_name = mxGetFieldNameByNumber(prhs[1], i);
                        mxArray *field_val = mxGetFieldByNumber(prhs[1], 0, i);
                        int val = 0;
                        char string_val[128];
                        if (mxIsChar(field_val))
                            mxGetString(field_val, string_val, 128);
                        else
                            val = *((double*) mxGetPr(field_val));
                        //mexPrintf(" val %d  \n", val);
                        if (!strcmp(field_name, "depth_mode")) initParams.depth_mode = static_cast<sl::DEPTH_MODE> (val);
                        if (!strcmp(field_name, "coordinate_units")) initParams.coordinate_units = static_cast<sl::UNIT> (val);
                        if (!strcmp(field_name, "coordinate_system")) initParams.coordinate_system = static_cast<sl::COORDINATE_SYSTEM> (val);
                        if (!strcmp(field_name, "sdk_verbose")) initParams.sdk_verbose = val;
                        if (!strcmp(field_name, "sdk_gpu_id")) initParams.sdk_gpu_id = val;
                        if (!strcmp(field_name, "depth_minimum_distance")) initParams.depth_minimum_distance = val;
                        if (!strcmp(field_name, "camera_disable_self_calib")) initParams.camera_disable_self_calib = val;
                        if (!strcmp(field_name, "camera_image_flip")) initParams.camera_image_flip = val;
                        if (!strcmp(field_name, "svo_filename")) initParams.svo_input_filename = string_val;
                        if (!strcmp(field_name, "camera_resolution")) initParams.camera_resolution = static_cast<sl::RESOLUTION> (val);
                        if (!strcmp(field_name, "camera_fps")) initParams.camera_fps = val;
                        if (!strcmp(field_name, "svo_real_time_mode")) initParams.svo_real_time_mode = val;
                        if (!strcmp(field_name, "camera_image_flip")) initParams.camera_image_flip = val;
                        if (!strcmp(field_name, "depth_stabilization")) initParams.depth_stabilization = val;
                        if (!strcmp(field_name, "enable_right_side_measure")) initParams.enable_right_side_measure = val;
                    }
                }
                err = zedCam->open(initParams);
            } else {
                err = zedCam->open();
            }
            // we return the string associated with the error
            plhs[0] = mxCreateString(sl::errorCode2str(err).c_str());
        }
    }


    else if (!strcmp(command, "close")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            zedCam->close();
        }
    }

    else if (!strcmp(command, "isOpened")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            double val = zedCam->isOpened();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
        }
    }

    else if (!strcmp(command, "grab")) {
        if (checkZED()) {
            sl::RuntimeParameters grabParams;
            // check if we have ONE argument
            if (checkParams(nrhs, 1)) {
                // check if we have a parameter structure
                if (mxIsStruct(prhs[1])) {
                    // for all fields of parameter structure overwrite parameters
                    for (int32_t i = 0; i < mxGetNumberOfFields(prhs[1]); i++) {
                        const char *field_name = mxGetFieldNameByNumber(prhs[1], i);
                        mxArray *field_val = mxGetFieldByNumber(prhs[1], 0, i);
                        int val = 0;
                        char string_val[128];
                        if (mxIsChar(field_val))
                            mxGetString(field_val, string_val, 128);
                        else
                            val = *((double*) mxGetPr(field_val));
                        if (!strcmp(field_name, "sensing_mode")) grabParams.sensing_mode = static_cast<sl::SENSING_MODE> (val);
                        if (!strcmp(field_name, "enable_depth")) grabParams.enable_depth = val;
                        if (!strcmp(field_name, "enable_point_cloud")) grabParams.enable_point_cloud = val;
                        if (!strcmp(field_name, "measure3D_reference_frame")) grabParams.measure3D_reference_frame = static_cast<sl::REFERENCE_FRAME> (val);
                    }
                }
                zedCam->grab(grabParams);
            } else {
                zedCam->grab();
            }
        }
    }

    else if (!strcmp(command, "retrieveImage")) {
        if (checkZED()) {
            sl::VIEW view = sl::VIEW_LEFT;
            if (checkParams(nrhs, 1, 3)) {
                double *ptr_ = mxGetPr(prhs[1]);
                int val = ptr_[0];
                if (val < sl::VIEW_LAST)
                    view = static_cast<sl::VIEW>(val);
                else {
                    mexErrMsgTxt("Can't find this image");
                    deleteOnFail();
                    return;
                }
            }
            int width = 0, height = 0;
            if (nrhs == 4) {
                double *ptr_ = mxGetPr(prhs[2]);
                width = ptr_[0];
                ptr_ = mxGetPr(prhs[3]);
                height = ptr_[0];
            }

            cv::Mat image_rgb;
            sl::Mat tmp;
            zedCam->retrieveImage(tmp, view, sl::MEM_CPU, width, height);

            int cv_type = CV_8UC4;

            if (view == sl::VIEW_LEFT_GRAY || view == sl::VIEW_RIGHT_GRAY || view == sl::VIEW_LEFT_UNRECTIFIED_GRAY || view == sl::VIEW_RIGHT_UNRECTIFIED_GRAY)
                cv_type = CV_8UC1;

            cv::Mat cv_tmp = cv::Mat(tmp.getHeight(), tmp.getWidth(), cv_type, tmp.getPtr<sl::uchar1>());

            if (cv_type == CV_8UC4)
                cv::cvtColor(cv_tmp, image_rgb, CV_RGBA2RGB);
            plhs[0] = CvMat_to_new_mxArr(image_rgb);
        }
    }

    else if (!strcmp(command, "retrieveMeasure")) {
        if (checkZED()) {
            sl::MEASURE measure = sl::MEASURE_DEPTH;
            if (checkParams(nrhs, 1, 3)) {
                double *ptr_ = mxGetPr(prhs[1]);
                int val = ptr_[0];
                if (val < sl::VIEW_LAST)
                    measure = static_cast<sl::MEASURE>(val);
                else {
                    mexErrMsgTxt("Can't find this measure");
                    deleteOnFail();
                    return;
                }
            }

            int width = 0, height = 0;
            if (nrhs == 4) {
                double *ptr_ = mxGetPr(prhs[2]);
                width = ptr_[0];
                ptr_ = mxGetPr(prhs[3]);
                height = ptr_[0];
            }

            cv::Mat measure_mat;
            sl::Mat tmp;
            zedCam->retrieveMeasure(tmp, measure, sl::MEM_CPU, width, height);
            int cv_type = CV_32FC1;

            if ((measure >= sl::MEASURE_XYZ &&  measure <= sl::MEASURE_NORMALS) || (measure >= sl::MEASURE_XYZ_RIGHT &&  measure < sl::MEASURE_LAST))
                cv_type = CV_32FC4;

            measure_mat = cv::Mat(tmp.getHeight(), tmp.getWidth(), cv_type, tmp.getPtr<sl::uchar1>());

            if (cv_type == CV_32FC4) {
                std::vector<cv::Mat> mat_v;
                cv::split(measure_mat, mat_v);
                plhs[0] = CvMat_to_new_mxArr(mat_v[0]);
                plhs[1] = CvMat_to_new_mxArr(mat_v[1]);
                plhs[2] = CvMat_to_new_mxArr(mat_v[2]);
                plhs[3] = CvMat_to_new_mxArr(mat_v[3]);
            } else
                plhs[0] = CvMat_to_new_mxArr(measure_mat);
        }
    }

    else if (!strcmp(command, "setConfidenceThreshold")) {
        if (checkZED() && checkParams(nrhs, 1)) {
            double *ptr_ = mxGetPr(prhs[1]);
            int val = ptr_[0];
            zedCam->setConfidenceThreshold(val);
        }
    }

    else if (!strcmp(command, "getConfidenceThreshold")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            double val = zedCam->getConfidenceThreshold();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
        }
    }

    else if (!strcmp(command, "getCUDAContext")) {
        notAvailable();
    }

    else if (!strcmp(command, "getResolution")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            double ptr_size[2];
            ptr_size[0] = zedCam->getResolution().width;
            ptr_size[1] = zedCam->getResolution().height;
            plhs[0] = mxCreateDoubleMatrix(2, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), ptr_size, 2 * sizeof(double));
        }
    }

    else if (!strcmp(command, "setDepthMaxRangeValue")) {
        if (checkZED() && checkParams(nrhs, 1)) {
            double *ptr_ = mxGetPr(prhs[1]);
            int val = ptr_[0];
            zedCam->setDepthMaxRangeValue(val);
        }
    }

    else if (!strcmp(command, "getDepthMaxRangeValue")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            double val = zedCam->getDepthMaxRangeValue();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
        }
    }

    else if (!strcmp(command, "getDepthMinRangeValue")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            double val = zedCam->getDepthMinRangeValue();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
        }
    }

    else if (!strcmp(command, "setSVOPosition")) {
        if (checkZED() && checkParams(nrhs, 1)) {
            double *ptr_ = mxGetPr(prhs[1]);
            int val = ptr_[0];
            zedCam->setSVOPosition(val);
        }
    }

    else if (!strcmp(command, "getSVOPosition")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            double val = zedCam->getSVOPosition();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
        }
    }

    else if (!strcmp(command, "getSVONumberOfFrames")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            double val = zedCam->getSVONumberOfFrames();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
        }
    }

    else if (!strcmp(command, "setCameraSettings")) {
        if (checkZED() && checkParams(nrhs, 2)) {
            char settingName[64];
            mxGetString(prhs[1], settingName, 64);
            double *ptr_ = mxGetPr(prhs[2]);
            int val = ptr_[0];
            bool useDefault = false;
            if (val == -1)
                useDefault = true;

            if (!strcmp(settingName, "brightness"))
                zedCam->setCameraSettings(sl::CAMERA_SETTINGS_BRIGHTNESS, static_cast<sl::CAMERA_SETTINGS>(val), useDefault);
            else if (!strcmp(settingName, "contrast"))
                zedCam->setCameraSettings(sl::CAMERA_SETTINGS_CONTRAST, static_cast<sl::CAMERA_SETTINGS>(val), useDefault);
            else if (!strcmp(settingName, "hue"))
                zedCam->setCameraSettings(sl::CAMERA_SETTINGS_HUE, static_cast<sl::CAMERA_SETTINGS>(val), useDefault);
            else if (!strcmp(settingName, "saturation"))
                zedCam->setCameraSettings(sl::CAMERA_SETTINGS_SATURATION, static_cast<sl::CAMERA_SETTINGS>(val), useDefault);
            else if (!strcmp(settingName, "gain"))
                zedCam->setCameraSettings(sl::CAMERA_SETTINGS_GAIN, static_cast<sl::CAMERA_SETTINGS>(val), useDefault);
            else if (!strcmp(settingName, "exposure"))
                zedCam->setCameraSettings(sl::CAMERA_SETTINGS_EXPOSURE, static_cast<sl::CAMERA_SETTINGS>(val), useDefault);
            else if (!strcmp(settingName, "whitebalance")) {
                zedCam->setCameraSettings(sl::CAMERA_SETTINGS_WHITEBALANCE, static_cast<sl::CAMERA_SETTINGS>(val), useDefault);
            } else {
                mexErrMsgTxt("Unkown CameraSettings");
                deleteOnFail();
            }
        }
    }

    else if (!strcmp(command, "getCameraSettings")) {
        if (checkZED() && checkParams(nrhs, 1)) {
            char settingName[64];
            mxGetString(prhs[1], settingName, 64);
            double val = 0;

            if (!strcmp(settingName, "brightness"))
                val = zedCam->getCameraSettings(sl::CAMERA_SETTINGS_BRIGHTNESS);
            else if (!strcmp(settingName, "contrast"))
                val = zedCam->getCameraSettings(sl::CAMERA_SETTINGS_CONTRAST);
            else if (!strcmp(settingName, "hue"))
                val = zedCam->getCameraSettings(sl::CAMERA_SETTINGS_HUE);
            else if (!strcmp(settingName, "saturation"))
                val = zedCam->getCameraSettings(sl::CAMERA_SETTINGS_SATURATION);
            else if (!strcmp(settingName, "gain"))
                val = zedCam->getCameraSettings(sl::CAMERA_SETTINGS_GAIN);
            else if (!strcmp(settingName, "exposure"))
                val = zedCam->getCameraSettings(sl::CAMERA_SETTINGS_EXPOSURE);
            else if (!strcmp(settingName, "whitebalance")) {
                val = zedCam->getCameraSettings(sl::CAMERA_SETTINGS_WHITEBALANCE);
            } else {
                mexErrMsgTxt("Unkown CameraSettings");
                deleteOnFail();
            }

            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
        }
    }

    else if (!strcmp(command, "getCurrentFPS")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            double val = zedCam->getCurrentFPS();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
        }
    }

    else if (!strcmp(command, "getCameraFPS")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            double val = zedCam->getCameraFPS();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
        }
    }

    else if (!strcmp(command, "setCameraFPS")) {
        if (checkZED() && checkParams(nrhs, 1)) {
            double *ptr_ = mxGetPr(prhs[1]);
            int val = ptr_[0];
            zedCam->setCameraFPS(val);
        }
    }

    else if (!strcmp(command, "getCameraTimestamp")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            double val = zedCam->getCameraTimestamp();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
        }
    }

    else if (!strcmp(command, "getCurrentTimestamp")) {
        if (checkZED()) {
            double val = zedCam->getCurrentTimestamp();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
        }
    }

    else if (!strcmp(command, "getFrameDroppedCount")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            double val = zedCam->getFrameDroppedCount();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
        }
    }

    else if (!strcmp(command, "getCameraInformation")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            sl::CameraInformation camInfo = zedCam->getCameraInformation();
            mxArray *pLeft;
            mxArray *pRight;
            plhs[0] = mxCreateStructMatrix(1, 1, 12, fieldsStruct);
            pLeft = mxCreateStructMatrix(1, 1, 10, fieldsIntra);
            pRight = mxCreateStructMatrix(1, 1, 10, fieldsIntra);
            mxSetField(plhs[0], 0, "serial_number", mxCreateDoubleScalar(camInfo.serial_number));
            mxSetField(plhs[0], 0, "firmware_version", mxCreateDoubleScalar(camInfo.firmware_version));
            mxSetField(plhs[0], 0, "baseline", mxCreateDoubleScalar(camInfo.calibration_parameters.T.x));
            mxSetField(plhs[0], 0, "convergence", mxCreateDoubleScalar(camInfo.calibration_parameters.R.y));
            mxSetField(plhs[0], 0, "r_x", mxCreateDoubleScalar(camInfo.calibration_parameters.R.x));
            mxSetField(plhs[0], 0, "r_y", mxCreateDoubleScalar(camInfo.calibration_parameters.R.y));
            mxSetField(plhs[0], 0, "r_z", mxCreateDoubleScalar(camInfo.calibration_parameters.R.z));
            mxSetField(plhs[0], 0, "t_x", mxCreateDoubleScalar(camInfo.calibration_parameters.T.x));
            mxSetField(plhs[0], 0, "t_y", mxCreateDoubleScalar(camInfo.calibration_parameters.T.y));
            mxSetField(plhs[0], 0, "t_z", mxCreateDoubleScalar(camInfo.calibration_parameters.T.z));

            mxSetField(pLeft, 0, "cx", mxCreateDoubleScalar(camInfo.calibration_parameters.left_cam.cx));
            mxSetField(pLeft, 0, "cy", mxCreateDoubleScalar(camInfo.calibration_parameters.left_cam.cy));
            mxArray* matLeftDisto = mxCreateDoubleMatrix(1, 5, mxREAL);
            memcpy(mxGetPr(matLeftDisto), &(camInfo.calibration_parameters.left_cam.disto), 5 * sizeof(double));
            mxSetField(pLeft, 0, "disto", matLeftDisto);
            mxSetField(pLeft, 0, "d_fov", mxCreateDoubleScalar(camInfo.calibration_parameters.left_cam.d_fov));
            mxSetField(pLeft, 0, "fx", mxCreateDoubleScalar(camInfo.calibration_parameters.left_cam.fx));
            mxSetField(pLeft, 0, "fy", mxCreateDoubleScalar(camInfo.calibration_parameters.left_cam.fy));
            mxSetField(pLeft, 0, "h_fov", mxCreateDoubleScalar(camInfo.calibration_parameters.left_cam.h_fov));
            mxSetField(pLeft, 0, "width", mxCreateDoubleScalar(camInfo.calibration_parameters.left_cam.image_size.width));
            mxSetField(pLeft, 0, "height", mxCreateDoubleScalar(camInfo.calibration_parameters.left_cam.image_size.height));
            mxSetField(pLeft, 0, "v_fov", mxCreateDoubleScalar(camInfo.calibration_parameters.left_cam.v_fov));
            mxSetField(plhs[0], 0, "left_cam", pLeft);

            mxSetField(pRight, 0, "cx", mxCreateDoubleScalar(camInfo.calibration_parameters.right_cam.cx));
            mxSetField(pRight, 0, "cy", mxCreateDoubleScalar(camInfo.calibration_parameters.right_cam.cy));
            mxArray* matRightDisto = mxCreateDoubleMatrix(1, 5, mxREAL);
            memcpy(mxGetPr(matRightDisto), &(camInfo.calibration_parameters.right_cam.disto), 5 * sizeof(double));
            mxSetField(pRight, 0, "disto", matRightDisto);
            mxSetField(pRight, 0, "d_fov", mxCreateDoubleScalar(camInfo.calibration_parameters.right_cam.d_fov));
            mxSetField(pRight, 0, "fx", mxCreateDoubleScalar(camInfo.calibration_parameters.right_cam.fx));
            mxSetField(pRight, 0, "fy", mxCreateDoubleScalar(camInfo.calibration_parameters.right_cam.fy));
            mxSetField(pRight, 0, "h_fov", mxCreateDoubleScalar(camInfo.calibration_parameters.right_cam.h_fov));
            mxSetField(pRight, 0, "width", mxCreateDoubleScalar(camInfo.calibration_parameters.right_cam.image_size.width));
            mxSetField(pRight, 0, "height", mxCreateDoubleScalar(camInfo.calibration_parameters.right_cam.image_size.height));
            mxSetField(pRight, 0, "v_fov", mxCreateDoubleScalar(camInfo.calibration_parameters.right_cam.v_fov));
            mxSetField(plhs[0], 0, "right_cam", pRight);
        }
    }

    else if (!strcmp(command, "getSelfCalibrationState")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            double val = zedCam->getSelfCalibrationState();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
        }
    }

    else if (!strcmp(command, "resetSelfCalibration")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            zedCam->resetSelfCalibration();
        }
    }

    else if (!strcmp(command, "enableTracking")) {
        if (checkZED()) {
            double res = 0;
            sl::TrackingParameters  trackParams;
            sl::Transform tranform;
            // check if we have ONE argument
            if (checkParams(nrhs, 1)) {
                // check if we have a parameter structure
                if (mxIsStruct(prhs[1])) {
                    // for all fields of parameter structure overwrite parameters
                    for (int32_t i = 0; i < mxGetNumberOfFields(prhs[1]); i++) {
                        const char *field_name = mxGetFieldNameByNumber(prhs[1], i);
                        mxArray *field_val = mxGetFieldByNumber(prhs[1], 0, i);
                        int val = 0;
                        char string_val[128];
                        if (mxIsChar(field_val))
                            mxGetString(field_val, string_val, 128);
                        else
                            val = *((double*) mxGetPr(field_val));
                        if (!strcmp(field_name, "area_file_path")) trackParams.area_file_path = string_val;
                        if (!strcmp(field_name, "enable_spatial_memory")) trackParams.enable_spatial_memory = val;
                        if (!strcmp(field_name, "initial_world_transform")) {
                            for (int col = 0; col < mxGetN(field_val); col++)
                                for (int row = 0; row < mxGetM(field_val); row++)
                                    tranform(row, col) = (mxGetPr(field_val))[row + col*mxGetM(field_val)];
                            trackParams.initial_world_transform = tranform;
                        }
                    }
                }
                res = zedCam->enableTracking(trackParams);
            } else {
                res = zedCam->enableTracking();
            }
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &res, 1 * sizeof(double));
        }
    }

    else if (!strcmp(command, "getPosition")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            cv::Mat position(4, 4, CV_32FC1);
            sl::Pose path;
            if (nrhs == 2) {
                double *ptr_ = mxGetPr(prhs[1]);
                int val = ptr_[0];
                zedCam->getPosition(path, static_cast<sl::REFERENCE_FRAME>(val));
            } else
                zedCam->getPosition(path);

            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    position.at<float>(i, j) = path.pose_data(i, j);
                }
            }
            plhs[0] = CvMat_to_new_mxArr(position);
        }
    }

    else if (!strcmp(command, "getAreaExportState")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            double val = zedCam->getAreaExportState();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
        }
    }

    else if (!strcmp(command, "disableTracking")) {
        if (checkZED()) {
            if (checkParams(nrhs, 1) && mxIsChar(prhs[1])) {
                char area_file_path[128];
                mxGetString(prhs[1], area_file_path, 128);
                zedCam->disableTracking(area_file_path);
            } else
                zedCam->disableTracking();
        }
    }

    else if (!strcmp(command, "resetTracking")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            sl::Transform path;
            zedCam->resetTracking(path);
        }
    }

    else if (!strcmp(command, "enableSpatialMapping")) {
        if (checkZED()) {
            double res = 0;
            sl::SpatialMappingParameters  mapParams;

            // check if we have ONE argument
            if (checkParams(nrhs, 1)) {
                // check if we have a parameter structure
                if (mxIsStruct(prhs[1])) {
                    // for all fields of parameter structure overwrite parameters
                    for (int32_t i = 0; i < mxGetNumberOfFields(prhs[1]); i++) {
                        const char *field_name = mxGetFieldNameByNumber(prhs[1], i);
                        mxArray *field_val = mxGetFieldByNumber(prhs[1], 0, i);
                        int val = 0;
                        char string_val[128];
                        if (mxIsChar(field_val))
                            mxGetString(field_val, string_val, 128);
                        else
                            val = *((double*) mxGetPr(field_val));
                        //mexPrintf(" val %d  \n", val);
                        if (!strcmp(field_name, "area_file_path")) mapParams.max_memory_usage = val;
                        if (!strcmp(field_name, "range_meter_max")) mapParams.range_meter.second = val;
                        if (!strcmp(field_name, "range_meter_min")) mapParams.range_meter.first = val;
                        if (!strcmp(field_name, "resolution_meter")) mapParams.resolution_meter = val;
                        if (!strcmp(field_name, "save_texture")) mapParams.save_texture = val;
                    }
                }
                res = zedCam->enableSpatialMapping(mapParams);
            } else {
                res = zedCam->enableSpatialMapping();
            }
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &res, 1 * sizeof(double));
        }
    }

    else if (!strcmp(command, "pauseSpatialMapping")) {
        if (checkZED() && checkParams(nrhs, 1)) {
            double *ptr_ = mxGetPr(prhs[1]);
            int val = ptr_[0];
            zedCam->pauseSpatialMapping(val);
        }
    }

    else if (!strcmp(command, "getSpatialMappingState")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            double val = zedCam->getSpatialMappingState();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
        }
    }

    else if (!strcmp(command, "extractAllMesh")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            notAvailable();
        }
    }

    else if (!strcmp(command, "requestMeshAsync")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            zedCam->requestMeshAsync();
        }
    }

    else if (!strcmp(command, "getMeshRequestStatusAsync")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            double val = zedCam->getMeshRequestStatusAsync();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
        }
    }

    else if (!strcmp(command, "retrieveMeshAsync")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            notAvailable();
        }
    }

    else if (!strcmp(command, "disableSpatialMapping")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            zedCam->disableSpatialMapping();
        }
    }

    else if (!strcmp(command, "enableRecording")) {
        if (checkZED()) {
            if (checkParams(nrhs, 1) || checkParams(nrhs, 2)) {
                if (mxIsChar(prhs[1])) {
                    char svo_path[128];
                    mxGetString(prhs[1], svo_path, 128);
                    if (checkParams(nrhs, 2)) {
                        double *ptr_ = mxGetPr(prhs[2]);
                        int val = ptr_[0];
                        zedCam->enableRecording(svo_path, static_cast<sl::SVO_COMPRESSION_MODE>(val));
                    } else
                        zedCam->enableRecording(svo_path);
                }
            }
        }
    }

    else if (!strcmp(command, "record")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            sl::RecordingState state = zedCam->record();
            plhs[0] = mxCreateStructMatrix(1, 1, 5, fieldsStruct);
            mxSetField(plhs[0], 0, "average_compression_ratio", mxCreateDoubleScalar(state.average_compression_ratio));
            mxSetField(plhs[0], 0, "average_compression_time", mxCreateDoubleScalar(state.average_compression_time));
            mxSetField(plhs[0], 0, "current_compression_ratio", mxCreateDoubleScalar(state.current_compression_ratio));
            mxSetField(plhs[0], 0, "current_compression_time", mxCreateDoubleScalar(state.current_compression_time));
            mxSetField(plhs[0], 0, "status", mxCreateDoubleScalar(state.status));
        }
    }

    else if (!strcmp(command, "disableRecording")) {
        if (checkZED() && checkParams(nrhs, 0))
            zedCam->disableTracking();
    }


    else if (!strcmp(command, "getSDKVersion")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            const char* version = zedCam->getSDKVersion();
            plhs[0] = mxCreateString(version);
        }
    }

    else if (!strcmp(command, "isZEDconnected")) {
        if (checkZED() && checkParams(nrhs, 0)) {
            double val = zedCam->isZEDconnected();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
        }
    } else {
        mexErrMsgTxt("Can't find the specified function");
        deleteOnFail();
    }
}
