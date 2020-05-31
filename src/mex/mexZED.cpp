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

#ifdef _DEBUG
#error Select Release mode for compilation
#endif

// global var.
static sl::Camera *zedCam = nullptr;

// C++ struct
const char* fieldsIntra[] = {"cx", "cy", "disto", "d_fov", "fx" , "fy" , "h_fov" , "width" , "height" , "v_fov"};
const char* fieldsParameters[] = {"serial_number", "firmware_version", "R", "t", "left_cam", "right_cam", "model"};
const char* fieldsIMU[] = {"pose","angular_velocity", "linear_acceleration"};
const char* fieldsObject[] = { "time_stamp", "is_new", "is_tracked", "object_list" };
const char* fieldsObjectData[] = { "id","label", "tracking_state", "position", "bounding_box_2d", "bounding_box_3d"};

// Interop. OpenCV-Matlab matrix type (float)
template<class InputIterator, class OutputIterator>
OutputIterator fctCopy(InputIterator first, InputIterator last, OutputIterator result) {
    while(first != last) {
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
    if(CV_32FC1 == TYPE)
        return floatmat_to_mat<CV_32FC1>(matrix);
    else if(CV_8UC3 == TYPE)
        return rgbimage_to_mat<IPL_DEPTH_8U>(matrix);

    return mxCreateDoubleMatrix(0, 0, mxREAL);
}

inline bool checkZED() {
    if(zedCam)
        return true;
    else 
		mexWarnMsgTxt("ZED is not initialized");    
    return false;
}

inline bool checkParams(int params, int required, bool print_warning = true) {
    if(params - 1 != required) {
		if (print_warning) {
			std::string error = "Invalid parameter number, " + std::to_string(required) + " required, " + std::to_string(params - 1) + " given.";
			mexWarnMsgTxt(error.c_str());
		}
        return false;
    }
    return true;
}

inline bool checkParams(int params, int required_1, int required_2, bool print_warning = true) {
    if((params - 1 != required_1) && (params - 1 != required_2)) {
		if (print_warning) {
			std::string error = "Invalid parameter number, " + std::to_string(required_1) + " or " + std::to_string(required_2) + " required, " + std::to_string(params - 1) + " given.";
			mexWarnMsgTxt(error.c_str());
		}
        return false;
    }
    return true;
}

inline bool checkParamsAtLeast(int params, int number, bool print_warning = true) {
    if (params - 1 < number) {
        if (print_warning) {
            std::string error = "Invalid parameter number, " + std::to_string(number) + " required, " + std::to_string(params - 1) + " given.";
            mexWarnMsgTxt(error.c_str());
        }
        return false;
    }
    return true;
}

void fillCameraParam(mxArray *pArray, sl::CameraParameters &param) {
    mxSetField(pArray, 0, fieldsIntra[0], mxCreateDoubleScalar(param.cx));
    mxSetField(pArray, 0, fieldsIntra[1], mxCreateDoubleScalar(param.cy));
    mxArray* matDisto = mxCreateDoubleMatrix(1, 5, mxREAL);
    memcpy(mxGetPr(matDisto), &(param.disto), 5 * sizeof(double));
    mxSetField(pArray, 0, fieldsIntra[2], matDisto);
    mxSetField(pArray, 0, fieldsIntra[3], mxCreateDoubleScalar(param.d_fov));
    mxSetField(pArray, 0, fieldsIntra[4], mxCreateDoubleScalar(param.fx));
    mxSetField(pArray, 0, fieldsIntra[5], mxCreateDoubleScalar(param.fy));
    mxSetField(pArray, 0, fieldsIntra[6], mxCreateDoubleScalar(param.h_fov));
    mxSetField(pArray, 0, fieldsIntra[7], mxCreateDoubleScalar(param.image_size.width));
    mxSetField(pArray, 0, fieldsIntra[8], mxCreateDoubleScalar(param.image_size.height));
    mxSetField(pArray, 0, fieldsIntra[9], mxCreateDoubleScalar(param.v_fov));
}

const char* point3D[] = { "x", "y", "z"};
const char* point2D[] = { "u", "v"};
//const char* fieldsObjectData[] = { "id", "label", "tracking_state", "position", "bounding_box_2d", "bounding_box_3d" };
void fillObjectData(mxArray* pArray, sl::ObjectData& obj, int idx) {

    mxSetField(pArray, idx, fieldsObjectData[0], mxCreateDoubleScalar(static_cast<double>(obj.id)));
    mxSetField(pArray, idx, fieldsObjectData[1], mxCreateDoubleScalar(static_cast<double>(obj.label)));
    mxSetField(pArray, idx, fieldsObjectData[2], mxCreateDoubleScalar(static_cast<double>(obj.tracking_state)));

    mxArray* position3d = mxCreateStructMatrix(1, 1, 3, point3D);
    mxSetField(position3d, 0, "x", mxCreateDoubleScalar(obj.position.x));
    mxSetField(position3d, 0, "y", mxCreateDoubleScalar(obj.position.y));
    mxSetField(position3d, 0, "z", mxCreateDoubleScalar(obj.position.z));
    mxSetField(pArray, idx, fieldsObjectData[3], position3d);
    
    mxArray* bb2d = mxCreateStructMatrix(1, obj.bounding_box_2d.size(), 2, point2D);
    for (int i = 0; i < obj.bounding_box_2d.size(); i++){
        mxSetField(bb2d, i, "u", mxCreateDoubleScalar(obj.bounding_box_2d[i].x));
        mxSetField(bb2d, i, "v", mxCreateDoubleScalar(obj.bounding_box_2d[i].y));
    }
    mxSetField(pArray, idx, fieldsObjectData[4], bb2d);

    mxArray* bb3d = mxCreateStructMatrix(1, obj.bounding_box.size(), 3, point3D);
    for (int i = 0; i < obj.bounding_box.size(); i++) {
        mxSetField(bb3d, i, "x", mxCreateDoubleScalar(obj.bounding_box[i].x));
        mxSetField(bb3d, i, "y", mxCreateDoubleScalar(obj.bounding_box[i].y));
        mxSetField(bb3d, i, "z", mxCreateDoubleScalar(obj.bounding_box[i].z));
    }
    mxSetField(pArray, idx, fieldsObjectData[5], bb3d);
}

//const char* fieldsObject[] = { "time_stamp", "is_new", "is_tracked", "object_list"};
void fillObjects(mxArray* pArray, sl::Objects& objs) {
    mxSetField(pArray, 0, fieldsObject[0], mxCreateDoubleScalar(static_cast<double>(objs.timestamp.getMilliseconds())));
    mxSetField(pArray, 0, fieldsObject[1], mxCreateDoubleScalar(static_cast<double>(objs.is_new)));
    mxSetField(pArray, 0, fieldsObject[2], mxCreateDoubleScalar(static_cast<double>(objs.is_tracked)));

    mxArray* obj = mxCreateStructMatrix(1, objs.object_list.size(), 6, fieldsObjectData);
    for (int i = 0; i < objs.object_list.size(); i++)
        fillObjectData(obj, objs.object_list[i], i);    
    mxSetField(pArray, 0, fieldsObject[3], obj);
}


/* MEX entry function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	sl::ERROR_CODE err = sl::ERROR_CODE::FAILURE;
    // each sub function is referred by a string 'command'
    char command[128];
    mxGetString(prhs[0], command, 128);

    if(!strcmp(command, "open")) {
        zedCam = new sl::Camera();
        sl::InitParameters initParams;
        // check if we have ONE argument
        if(checkParams(nrhs, 1)) {
            // check if we have a parameter structure
            if(mxIsStruct(prhs[1])) {
                // for all fields of parameter structure overwrite parameters
                for(int32_t i = 0; i < mxGetNumberOfFields(prhs[1]); i++) {
                    const char *field_name = mxGetFieldNameByNumber(prhs[1], i);
                    mxArray *field_val = mxGetFieldByNumber(prhs[1], 0, i);
                    int val = 0;
                    char string_val[128];
                    if(mxIsChar(field_val))
                        mxGetString(field_val, string_val, 128);
                    else
                        val = *((double*) mxGetPr(field_val));
                    //mexPrintf(" val %d  \n", val);
                    if(!strcmp(field_name, "depth_mode")) initParams.depth_mode = static_cast<sl::DEPTH_MODE> (val);
                    if(!strcmp(field_name, "coordinate_units")) initParams.coordinate_units = static_cast<sl::UNIT> (val);
                    if(!strcmp(field_name, "coordinate_system")) initParams.coordinate_system = static_cast<sl::COORDINATE_SYSTEM> (val);
                    if(!strcmp(field_name, "sdk_verbose")) initParams.sdk_verbose = val;
                    if(!strcmp(field_name, "sdk_gpu_id")) initParams.sdk_gpu_id = val;
                    if(!strcmp(field_name, "depth_minimum_distance")) initParams.depth_minimum_distance = val;
                    if(!strcmp(field_name, "depth_maximum_distance")) initParams.depth_maximum_distance = val;
                    if(!strcmp(field_name, "camera_disable_self_calib")) initParams.camera_disable_self_calib = val;
                    if(!strcmp(field_name, "camera_image_flip")) initParams.camera_image_flip = val;                        
                    if(!strcmp(field_name, "svo_input_filename")) initParams.input.setFromSVOFile(string_val);
                    if(!strcmp(field_name, "camera_resolution")) initParams.camera_resolution = static_cast<sl::RESOLUTION> (val);
                    if(!strcmp(field_name, "camera_fps")) initParams.camera_fps = val;
                    if(!strcmp(field_name, "svo_real_time_mode")) initParams.svo_real_time_mode = val;
                    if(!strcmp(field_name, "depth_stabilization")) initParams.depth_stabilization = val;
                    if(!strcmp(field_name, "enable_right_side_measure")) initParams.enable_right_side_measure = val;
                }
            }
        }
        err = zedCam->open(initParams);
        
        // we return the string associated with the error
        plhs[0] = mxCreateString(sl::toString(err).c_str());
    }
    
    else if(!strcmp(command, "close")) {
        if(checkZED() && checkParams(nrhs, 0)) {
            zedCam->close();
            delete zedCam;
            zedCam = nullptr;
            err = sl::ERROR_CODE::SUCCESS;
        }
    }

    else if(!strcmp(command, "isOpened")) {
        if(checkZED() && checkParams(nrhs, 0)) {
            double val = zedCam->isOpened();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
            err = sl::ERROR_CODE::SUCCESS;
        }
    }

    else if(!strcmp(command, "grab")) {
        if(checkZED()) {
            sl::RuntimeParameters grabParams;
			// check if we have a parameter structure
			if ((nrhs > 1) && (mxIsStruct(prhs[1]))) {
				// for all fields of parameter structure overwrite parameters
				for (int32_t i = 0; i < mxGetNumberOfFields(prhs[1]); i++) {
					const char *field_name = mxGetFieldNameByNumber(prhs[1], i);
					mxArray *field_val = mxGetFieldByNumber(prhs[1], 0, i);
					int val = 0;
					char string_val[128];
					if (mxIsChar(field_val))
						mxGetString(field_val, string_val, 128);
					else
						val = *((double*)mxGetPr(field_val));
					if (!strcmp(field_name, "sensing_mode")) grabParams.sensing_mode = static_cast<sl::SENSING_MODE> (val);
					if (!strcmp(field_name, "enable_depth")) grabParams.enable_depth = val;
					if (!strcmp(field_name, "measure3D_reference_frame")) grabParams.measure3D_reference_frame = static_cast<sl::REFERENCE_FRAME> (val);
				}
			}
			err = zedCam->grab(grabParams);
        }
		// we return the string associated with the error
		plhs[0] = mxCreateString(sl::toString(err).c_str());
    }

    else if(!strcmp(command, "retrieveImage")) {
        if(checkZED()) {
            sl::VIEW view = sl::VIEW::LEFT;
            if(checkParams(nrhs, 1, 3)) {
                double *ptr_ = mxGetPr(prhs[1]);
                int val = ptr_[0];
                if(val < static_cast<int>(sl::VIEW::LAST))
                    view = static_cast<sl::VIEW>(val);
                else {
					mexWarnMsgTxt("Unknown VIEW requested");
                    return;
                }
            }
            int width = 0, height = 0;
            if(nrhs == 4) {
                double *ptr_ = mxGetPr(prhs[2]);
                width = ptr_[0];
                ptr_ = mxGetPr(prhs[3]);
                height = ptr_[0];
            }

            cv::Mat image_rgb;
            sl::Mat tmp;
            err = zedCam->retrieveImage(tmp, view, sl::MEM::CPU, sl::Resolution(width, height));

            int cv_type = CV_8UC4;

            if(view == sl::VIEW::LEFT_GRAY || view == sl::VIEW::RIGHT_GRAY || view == sl::VIEW::LEFT_UNRECTIFIED_GRAY || view == sl::VIEW::RIGHT_UNRECTIFIED_GRAY)
                cv_type = CV_8UC1;

            cv::Mat cv_tmp = cv::Mat(tmp.getHeight(), tmp.getWidth(), cv_type, tmp.getPtr<sl::uchar1>());

            if(cv_type == CV_8UC4)
                cv::cvtColor(cv_tmp, image_rgb, cv::COLOR_RGBA2RGB);
            plhs[0] = CvMat_to_new_mxArr(image_rgb);
        }
    }

    else if(!strcmp(command, "retrieveMeasure")) {
        if(checkZED()) {
            sl::MEASURE measure = sl::MEASURE::DEPTH;
            if(checkParams(nrhs, 1, 3)) {
                double *ptr_ = mxGetPr(prhs[1]);
                int val = ptr_[0];
                if(val < static_cast<int>(sl::MEASURE::LAST))
                    measure = static_cast<sl::MEASURE>(val);
                else {
                    mexWarnMsgTxt("Unknown MEASURE requested");
                    return;
                }
            }

            int width = 0, height = 0;
            if(nrhs == 4) {
                double *ptr_ = mxGetPr(prhs[2]);
                width = ptr_[0];
                ptr_ = mxGetPr(prhs[3]);
                height = ptr_[0];
            }
            
            cv::Mat measure_mat;
            sl::Mat tmp;
            err = zedCam->retrieveMeasure(tmp, measure, sl::MEM::CPU, sl::Resolution(width, height));
            int cv_type = CV_32FC1;

            if((measure >= sl::MEASURE::XYZ &&  measure <= sl::MEASURE::NORMALS) || (measure >= sl::MEASURE::XYZ_RIGHT &&  measure < sl::MEASURE::LAST))
                cv_type = CV_32FC4;

            measure_mat = cv::Mat(tmp.getHeight(), tmp.getWidth(), cv_type, tmp.getPtr<sl::uchar1>());

            if(cv_type == CV_32FC4) {
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
            
    else if(!strcmp(command, "setSVOPosition")) {
        if(checkZED() && checkParams(nrhs, 1)) {
            double *ptr_ = mxGetPr(prhs[1]);
            int val = ptr_[0];
            zedCam->setSVOPosition(val);
            err = sl::ERROR_CODE::SUCCESS;
        }
    }

    else if(!strcmp(command, "getSVOPosition")) {
        if(checkZED() && checkParams(nrhs, 0)) {
            double val = zedCam->getSVOPosition();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
            err = sl::ERROR_CODE::SUCCESS;
        }
    }

    else if(!strcmp(command, "getSVONumberOfFrames")) {
        if(checkZED() && checkParams(nrhs, 0)) {
            double val = zedCam->getSVONumberOfFrames();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
            err = sl::ERROR_CODE::SUCCESS;
        }
    }

    else if(!strcmp(command, "setCameraSettings")) {
        if(checkZED() && checkParamsAtLeast(nrhs, 2)) {
            char settingName[64];
            mxGetString(prhs[1], settingName, 64);

            if (!strcmp(settingName, "aec_agc_roi")) {
                double* ptr_ = mxGetPr(prhs[2]);
                sl::Rect roi(ptr_[0], ptr_[1], ptr_[2], ptr_[3]);

                sl::SIDE side = sl::SIDE::BOTH;
                if (nrhs == 4) {
                    ptr_ = mxGetPr(prhs[3]);
                    side = static_cast<sl::SIDE>((int)ptr_[0]);
                }

                bool reset = false;
                if (nrhs == 5) {
                    ptr_ = mxGetPr(prhs[3]);
                    reset = (int)ptr_[0];
                }
                zedCam->setCameraSettings(sl::VIDEO_SETTINGS::AEC_AGC_ROI, roi, side, reset);
            } else {

                double* ptr_ = mxGetPr(prhs[2]);
                int val = ptr_[0];

                if (!strcmp(settingName, "brightness"))
                    zedCam->setCameraSettings(sl::VIDEO_SETTINGS::BRIGHTNESS, val);
                else if (!strcmp(settingName, "contrast"))
                    zedCam->setCameraSettings(sl::VIDEO_SETTINGS::CONTRAST, val);
                else if (!strcmp(settingName, "hue"))
                    zedCam->setCameraSettings(sl::VIDEO_SETTINGS::HUE, val);
                else if (!strcmp(settingName, "saturation"))
                    zedCam->setCameraSettings(sl::VIDEO_SETTINGS::SATURATION, val);
                else if (!strcmp(settingName, "gain"))
                    zedCam->setCameraSettings(sl::VIDEO_SETTINGS::GAIN, val);
                else if (!strcmp(settingName, "exposure"))
                    zedCam->setCameraSettings(sl::VIDEO_SETTINGS::EXPOSURE, val);
                else if (!strcmp(settingName, "aec_agc"))
                    zedCam->setCameraSettings(sl::VIDEO_SETTINGS::AEC_AGC, val);
                else if (!strcmp(settingName, "whitebalance"))
                    zedCam->setCameraSettings(sl::VIDEO_SETTINGS::WHITEBALANCE_TEMPERATURE, val);
                else
                    mexWarnMsgTxt("Unknown VIDEO SETTINGS");
            }
            err = sl::ERROR_CODE::SUCCESS;
        }
    }

    else if(!strcmp(command, "getCameraSettings")) {
        if(checkZED() && checkParams(nrhs, 1)) {
            char settingName[64];
            mxGetString(prhs[1], settingName, 64);
            
            if (!strcmp(settingName, "aec_agc_roi")) {
                sl::Rect roi;
                double val[4];
                zedCam->getCameraSettings(sl::VIDEO_SETTINGS::AEC_AGC_ROI, roi);
                val[0] = roi.x;
                val[1] = roi.y;
                val[2] = roi.width;
                val[3] = roi.height;
                plhs[0] = mxCreateDoubleMatrix(1, 4, mxREAL);
                memcpy(mxGetPr(plhs[0]), &val, 4 * sizeof(double));
            } else {
                double val = 0;
                if (!strcmp(settingName, "brightness"))
                    val = zedCam->getCameraSettings(sl::VIDEO_SETTINGS::BRIGHTNESS);
                else if (!strcmp(settingName, "contrast"))
                    val = zedCam->getCameraSettings(sl::VIDEO_SETTINGS::CONTRAST);
                else if (!strcmp(settingName, "hue"))
                    val = zedCam->getCameraSettings(sl::VIDEO_SETTINGS::HUE);
                else if (!strcmp(settingName, "saturation"))
                    val = zedCam->getCameraSettings(sl::VIDEO_SETTINGS::SATURATION);
                else if (!strcmp(settingName, "gain"))
                    val = zedCam->getCameraSettings(sl::VIDEO_SETTINGS::GAIN);
                else if (!strcmp(settingName, "exposure"))
                    val = zedCam->getCameraSettings(sl::VIDEO_SETTINGS::EXPOSURE);
                else if (!strcmp(settingName, "aec_agc"))
                    val = zedCam->getCameraSettings(sl::VIDEO_SETTINGS::AEC_AGC);
                else if (!strcmp(settingName, "whitebalance"))
                    val = zedCam->getCameraSettings(sl::VIDEO_SETTINGS::WHITEBALANCE_TEMPERATURE);
                else
                    mexWarnMsgTxt("Unknown VIDEO SETTINGS");

                plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
                memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
            }
            err = sl::ERROR_CODE::SUCCESS;
        }
    }

    else if(!strcmp(command, "getCurrentFPS")) {
        if(checkZED() && checkParams(nrhs, 0)) {
            double val = zedCam->getCurrentFPS();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
            err = sl::ERROR_CODE::SUCCESS;
        }
    }    

    else if(!strcmp(command, "getTimestamp")) {
        if(checkZED() && checkParams(nrhs, 1)) {
            double* ptr_ = mxGetPr(prhs[1]);
            int time_ref = ptr_[0];
            double val = zedCam->getTimestamp(static_cast<sl::TIME_REFERENCE>(time_ref)).getMilliseconds();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
            err = sl::ERROR_CODE::SUCCESS;
        }
    }

    else if(!strcmp(command, "getCameraInformation")) {
        if(checkZED() && checkParams(nrhs, 0)) {
            sl::CameraInformation camInfo = zedCam->getCameraInformation();
            mxArray *pLeft;
            mxArray *pRight;
            plhs[0] = mxCreateStructMatrix(1, 1, 7, fieldsParameters);
            pLeft = mxCreateStructMatrix(1, 1, 10, fieldsIntra);
            pRight = mxCreateStructMatrix(1, 1, 10, fieldsIntra);

            mxSetField(plhs[0], 0, fieldsParameters[0], mxCreateDoubleScalar(camInfo.serial_number));
            mxSetField(plhs[0], 0, fieldsParameters[1], mxCreateDoubleScalar(camInfo.camera_firmware_version));

            const mxClassID cid = cvm_traits<CV_32FC1>::CID;
            mxArray* rotation = mxCreateNumericMatrix(1, 3, cid, mxREAL);
            memcpy(mxGetPr(rotation), &(camInfo.calibration_parameters.R), sizeof(camInfo.calibration_parameters.R));
            mxSetField(plhs[0], 0, fieldsParameters[2], rotation);

            mxArray* translation = mxCreateNumericMatrix(1, 3, cid, mxREAL);
            memcpy(mxGetPr(translation), &(camInfo.calibration_parameters.T), sizeof(camInfo.calibration_parameters.T));
            mxSetField(plhs[0], 0, fieldsParameters[3], translation);

            fillCameraParam(pLeft, camInfo.calibration_parameters.left_cam);
            mxSetField(plhs[0], 0, fieldsParameters[4], pLeft);

            fillCameraParam(pRight, camInfo.calibration_parameters.right_cam);
            mxSetField(plhs[0], 0, fieldsParameters[5], pRight);

            mxSetField(plhs[0], 0, fieldsParameters[6], mxCreateDoubleScalar(static_cast<int>(camInfo.camera_model)));
            err = sl::ERROR_CODE::SUCCESS;
        }
    }

    else if(!strcmp(command, "enablePositionalTracking")) {
        if(checkZED()) {
            sl::PositionalTrackingParameters  trackParams;
            // check if we have a parameter structure
            if((nrhs > 1) && (mxIsStruct(prhs[1]))) {
				sl::Transform tranform;
                // for all fields of parameter structure overwrite parameters
                for(int32_t i = 0; i < mxGetNumberOfFields(prhs[1]); i++) {
                    const char *field_name = mxGetFieldNameByNumber(prhs[1], i);
                    mxArray *field_val = mxGetFieldByNumber(prhs[1], 0, i);
                    int val = 0;
                    char string_val[128];
                    if(mxIsChar(field_val))
                        mxGetString(field_val, string_val, 128);
                    else
                        val = *((double*) mxGetPr(field_val));
                    if(!strcmp(field_name, "area_file_path")) trackParams.area_file_path = string_val;
                    if(!strcmp(field_name, "enable_area_memory")) trackParams.enable_area_memory = val;
                    if(!strcmp(field_name, "initial_world_transform")) {
                        for(int col = 0; col < mxGetN(field_val); col++)
                            for(int row = 0; row < mxGetM(field_val); row++)
                                tranform(row, col) = (mxGetPr(field_val))[row + col*mxGetM(field_val)];
                        trackParams.initial_world_transform = tranform;
                    }
                }
            }
			err = zedCam->enablePositionalTracking(trackParams);            
        }
		// we return the string associated with the error
		plhs[0] = mxCreateString(sl::toString(err).c_str());
    }

    else if(!strcmp(command, "getPosition")) {
        if(checkZED() && checkParams(nrhs, 0)) {
            cv::Mat position(4, 4, CV_32FC1);
            sl::Pose path;
            if(nrhs == 2) {
                double *ptr_ = mxGetPr(prhs[1]);
                int val = ptr_[0];
                zedCam->getPosition(path, static_cast<sl::REFERENCE_FRAME>(val));
            } else
                zedCam->getPosition(path);

            for(int i = 0; i < 4; i++) {
                for(int j = 0; j < 4; j++) {
                    position.at<float>(i, j) = path.pose_data(i, j);
                }
            }
            plhs[0] = CvMat_to_new_mxArr(position);
            err = sl::ERROR_CODE::SUCCESS;
        }
    }

    else if(!strcmp(command, "getSensorsData")) {
        if(checkZED() && checkParams(nrhs, 0)) {
            sl::SensorsData sensors;
            sl::TIME_REFERENCE t_ref = sl::TIME_REFERENCE::IMAGE;
            if(nrhs == 2) {
                double *ptr_ = mxGetPr(prhs[1]);
                int val = ptr_[0];
                t_ref = static_cast<sl::TIME_REFERENCE>(val);
            }

            err = zedCam->getSensorsData(sensors, t_ref);

            plhs[0] = mxCreateStructMatrix(1, 1, 3, fieldsIMU);

            cv::Mat position(4, 4, CV_32FC1);
            for(int i = 0; i < 4; i++) {
                for(int j = 0; j < 4; j++) {
                    position.at<float>(i, j) = sensors.imu.pose(i, j);
                }
            }
            mxSetField(plhs[0], 0, fieldsIMU[0], CvMat_to_new_mxArr(position));

            const mxClassID cid = cvm_traits<CV_32FC1>::CID;
            mxArray* angular_v = mxCreateNumericMatrix(1, 3, cid, mxREAL);
            memcpy(mxGetPr(angular_v), &(sensors.imu.angular_velocity), 3 * sizeof(float));
            mxSetField(plhs[0], 0, fieldsIMU[1], angular_v);

            mxArray* linear_a = mxCreateNumericMatrix(1, 3, cid, mxREAL);
            memcpy(mxGetPr(linear_a), &(sensors.imu.linear_acceleration), 3 * sizeof(float));
            mxSetField(plhs[0], 0, fieldsIMU[2], linear_a);
        }
    }
   
    else if(!strcmp(command, "disablePositionalTracking")) {
        if(checkZED()) {
            if((nrhs > 1) && (mxIsChar(prhs[1]))) {
                char area_file_path[128];
                mxGetString(prhs[1], area_file_path, 128);
                zedCam->disablePositionalTracking(area_file_path);
            } else
                zedCam->disablePositionalTracking();
            err = sl::ERROR_CODE::SUCCESS;
        }
    }

    else if(!strcmp(command, "resetTracking")) {
        if(checkZED() && checkParams(nrhs, 0)) {
            sl::Transform path;
            err = zedCam->resetPositionalTracking(path);
        }
    }
    
    else if(!strcmp(command, "enableRecording")) {
        if(checkZED()) {
            if(checkParams(nrhs, 1, false) || checkParams(nrhs, 2, false)) {
                if(mxIsChar(prhs[1])) {
                    char svo_path[128];
                    mxGetString(prhs[1], svo_path, 128);
                    sl::RecordingParameters rec_param;
                    rec_param.video_filename = svo_path;
                    if(checkParams(nrhs, 2, false)) {
                        double *ptr_ = mxGetPr(prhs[2]);
                        rec_param.compression_mode = static_cast<sl::SVO_COMPRESSION_MODE>((int)ptr_[0]);
                    }
                    err = zedCam->enableRecording(rec_param);
                }
            }
        }
		// we return the string associated with the error
		plhs[0] = mxCreateString(sl::toString(err).c_str());
    }
    
    else if(!strcmp(command, "disableRecording")) {
        if (checkZED()) {
            zedCam->disableRecording();
            err = sl::ERROR_CODE::SUCCESS;
        }
    }

    else if (!strcmp(command, "enableObjectDetection")) {
        if (checkZED()) {
            sl::ObjectDetectionParameters param;
            err = zedCam->enableObjectDetection(param);
        }
    }

    else if (!strcmp(command, "retrieveObjects")) {
        if (checkZED()) {
            sl::ObjectDetectionRuntimeParameters rt_param;
            if (checkParams(nrhs, 1, false)) {
                double* ptr_ = mxGetPr(prhs[1]);
                rt_param.detection_confidence_threshold = (int)ptr_[0];
            }

            sl::Objects objs;
            err = zedCam->retrieveObjects(objs, rt_param);
            plhs[0] = mxCreateStructMatrix(1, 1, 4, fieldsObject);
            fillObjects(plhs[0], objs);
        }
    }

    else if (!strcmp(command, "disableObjectDetection")) {
        if (checkZED()) {
            zedCam->disableObjectDetection();
            err = sl::ERROR_CODE::SUCCESS;
        }
    }

    else if(!strcmp(command, "getSDKVersion")) {
        if(checkParams(nrhs, 0)) {
            const char* version = sl::Camera::getSDKVersion();
            plhs[0] = mxCreateString(version);
            err = sl::ERROR_CODE::SUCCESS;
        }
    }

    else if(!strcmp(command, "isZEDconnected")) {
        if(checkParams(nrhs, 0)) {
            double val = sl::Camera::getDeviceList().size();
            plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
            memcpy(mxGetPr(plhs[0]), &val, 1 * sizeof(double));
            err = sl::ERROR_CODE::SUCCESS;
        }
    }

    else if (!strcmp(command, "savePointCloudAs")) {
        if (checkParams(nrhs, 2) || checkParams(nrhs, 3)) {
            double *ptr_ = mxGetPr(prhs[1]);
            int format = ptr_[0];

            char save_path[128];
            if (mxIsChar(prhs[2])) {
                mxGetString(prhs[2], save_path, 128);
            }

            int color = 0;
            if (checkParams(nrhs, 3)) {
                double *ptr1_ = mxGetPr(prhs[3]);
                color = ptr1_[0];
            }
            sl::Mat pc;
            zedCam->retrieveMeasure(pc, sl::MEASURE::XYZBGRA);
            err = pc.write(save_path);
            plhs[0] = mxCreateString(sl::toString(err).c_str());
        }
    }

    else if (!strcmp(command, "saveDepthAs")) {
        if (checkParams(nrhs, 2) || checkParams(nrhs, 3)) {
            double *ptr_ = mxGetPr(prhs[1]);
            int format = ptr_[0];

            char save_path[128];
            if (mxIsChar(prhs[2])) {
                mxGetString(prhs[2], save_path, 128);
            }

            int color = 0;
            if (checkParams(nrhs, 3)) {
                double *ptr1_ = mxGetPr(prhs[3]);
                color = ptr1_[0];
            }

            sl::Mat depth;
            zedCam->retrieveMeasure(depth, sl::MEASURE::DEPTH);
            err = depth.write(save_path);
            plhs[0] = mxCreateString(sl::toString(err).c_str());
        }
    } 
    else {
        std::string undefined_fct("ZED SDK MEX does not contains specified function: " + std::string(command));
		mexWarnMsgTxt(undefined_fct.c_str());
        err = sl::ERROR_CODE::SUCCESS;
    }

    if (err != sl::ERROR_CODE::SUCCESS) {
        std::string error_msg("ZED SDK error when calling: " + std::string(command) + ", Error: "+ std::string(sl::toString(err)));
        mexWarnMsgTxt(error_msg.c_str());
    }
}
