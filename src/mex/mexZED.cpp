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
sl::SpatialMappingParameters sp_p;

const int BUFF_SIZE = 256;

// C++ struct
const char* fieldsParameters[] = {"serial_number", "firmware_version", "R", "t", "left_cam", "right_cam", "model"};

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
    else if (CV_16UC1 == TYPE)
        return floatmat_to_mat<CV_16UC1>(matrix);
    else if((CV_8UC3 == TYPE) || (CV_8UC1 == TYPE))
        return rgbimage_to_mat<IPL_DEPTH_8U>(matrix);
    return mxCreateDoubleMatrix(0, 0, mxREAL);
}

inline cv::Mat cvtMat(sl::Matrix3f &in){
    return cv::Mat(3, 3, CV_32FC1, in.r);
}

inline cv::Mat cvtMat(sl::Matrix4f &in){
    return cv::Mat(4, 4, CV_32FC1, in.m);
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

const char* fieldsIntra[] = {"cx", "cy", "disto", "d_fov", "fx" , "fy" , "h_fov" , "width" , "height" , "v_fov"};
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
const char* fieldsObjectData[] = { "id", "label", "tracking_state", "action_state", "confidence", "position", "velocity", "bounding_box_2d", "bounding_box_3d", "keypoint_2d", "keypoint"};
void fillObjectData(mxArray* pArray, sl::ObjectData& obj, int idx) {

    mxSetField(pArray, idx, fieldsObjectData[0], mxCreateDoubleScalar(static_cast<double>(obj.id)));
    mxSetField(pArray, idx, fieldsObjectData[1], mxCreateDoubleScalar(static_cast<double>(obj.label)));
    mxSetField(pArray, idx, fieldsObjectData[2], mxCreateDoubleScalar(static_cast<double>(obj.tracking_state)));
    mxSetField(pArray, idx, fieldsObjectData[3], mxCreateDoubleScalar(static_cast<double>(obj.action_state)));
    mxSetField(pArray, idx, fieldsObjectData[4], mxCreateDoubleScalar(static_cast<double>(obj.confidence)));

    mxArray* position = mxCreateStructMatrix(1, 1, 3, point3D);
    mxSetField(position, 0, "x", mxCreateDoubleScalar(obj.position.x));
    mxSetField(position, 0, "y", mxCreateDoubleScalar(obj.position.y));
    mxSetField(position, 0, "z", mxCreateDoubleScalar(obj.position.z));
    mxSetField(pArray, idx, fieldsObjectData[5], position);
    
    mxArray* velocity = mxCreateStructMatrix(1, 1, 3, point3D);
    mxSetField(velocity, 0, "x", mxCreateDoubleScalar(obj.velocity.x));
    mxSetField(velocity, 0, "y", mxCreateDoubleScalar(obj.velocity.y));
    mxSetField(velocity, 0, "z", mxCreateDoubleScalar(obj.velocity.z));
    mxSetField(pArray, idx, fieldsObjectData[6], velocity);
    
    mxArray* bounding_box_2d = mxCreateStructMatrix(1, obj.bounding_box_2d.size(), 2, point2D);
    for (int i = 0; i < obj.bounding_box_2d.size(); i++){
        mxSetField(bounding_box_2d, i, "u", mxCreateDoubleScalar(obj.bounding_box_2d[i].x));
        mxSetField(bounding_box_2d, i, "v", mxCreateDoubleScalar(obj.bounding_box_2d[i].y));
    }
    mxSetField(pArray, idx, fieldsObjectData[7], bounding_box_2d);

    mxArray* bounding_box = mxCreateStructMatrix(1, obj.bounding_box.size(), 3, point3D);
    for (int i = 0; i < obj.bounding_box.size(); i++) {
        mxSetField(bounding_box, i, "x", mxCreateDoubleScalar(obj.bounding_box[i].x));
        mxSetField(bounding_box, i, "y", mxCreateDoubleScalar(obj.bounding_box[i].y));
        mxSetField(bounding_box, i, "z", mxCreateDoubleScalar(obj.bounding_box[i].z));
    }
    mxSetField(pArray, idx, fieldsObjectData[8], bounding_box);

    
    mxArray* keypoint_2d = mxCreateStructMatrix(1, obj.keypoint_2d.size(), 2, point2D);
    for (int i = 0; i < obj.keypoint_2d.size(); i++){
        mxSetField(keypoint_2d, i, "u", mxCreateDoubleScalar(obj.keypoint_2d[i].x));
        mxSetField(keypoint_2d, i, "v", mxCreateDoubleScalar(obj.keypoint_2d[i].y));
    }
    mxSetField(pArray, idx, fieldsObjectData[9], keypoint_2d);

    mxArray* keypoint = mxCreateStructMatrix(1, obj.keypoint.size(), 3, point3D);
    for (int i = 0; i < obj.keypoint.size(); i++) {
        mxSetField(keypoint, i, "x", mxCreateDoubleScalar(obj.keypoint[i].x));
        mxSetField(keypoint, i, "y", mxCreateDoubleScalar(obj.keypoint[i].y));
        mxSetField(keypoint, i, "z", mxCreateDoubleScalar(obj.keypoint[i].z));
    }
    mxSetField(pArray, idx, fieldsObjectData[10], keypoint);

}

const char* fieldsObject[] = { "time_stamp", "is_new", "is_tracked", "object_list" };
void fillObjects(mxArray* pArray, sl::Objects& objs) {
    mxSetField(pArray, 0, fieldsObject[0], mxCreateDoubleScalar(static_cast<double>(objs.timestamp.getMilliseconds())));
    mxSetField(pArray, 0, fieldsObject[1], mxCreateDoubleScalar(static_cast<double>(objs.is_new)));
    mxSetField(pArray, 0, fieldsObject[2], mxCreateDoubleScalar(static_cast<double>(objs.is_tracked)));

    mxArray* obj = mxCreateStructMatrix(1, objs.object_list.size(), 11, fieldsObjectData);
    for (int i = 0; i < objs.object_list.size(); i++)
        fillObjectData(obj, objs.object_list[i], i);    
    mxSetField(pArray, 0, fieldsObject[3], obj);
}

const char* fieldsTemp[] = {"IMU", "BAROMETER", "ONBOARD_LEFT", "ONBOARD_RIGHT" };
void fillTemp(mxArray* pArray, sl::SensorsData::TemperatureData& temp) { 
    mxSetField(pArray, 0, fieldsTemp[0], mxCreateDoubleScalar(temp.temperature_map[sl::SensorsData::TemperatureData::SENSOR_LOCATION::IMU]));
    mxSetField(pArray, 0, fieldsTemp[1], mxCreateDoubleScalar(temp.temperature_map[sl::SensorsData::TemperatureData::SENSOR_LOCATION::BAROMETER]));
    mxSetField(pArray, 0, fieldsTemp[2], mxCreateDoubleScalar(temp.temperature_map[sl::SensorsData::TemperatureData::SENSOR_LOCATION::ONBOARD_LEFT]));
    mxSetField(pArray, 0, fieldsTemp[3], mxCreateDoubleScalar(temp.temperature_map[sl::SensorsData::TemperatureData::SENSOR_LOCATION::ONBOARD_RIGHT]));
}

const char* fieldsMag[] = {"is_available", "timestamp", "magnetic_field_calibrated", "magnetic_field_uncalibrated", "effective_rate" };
void fillMag(mxArray* pArray, sl::SensorsData::MagnetometerData& mag) { 
    //is_available
    mxSetField(pArray, 0, fieldsMag[0], mxCreateDoubleScalar(mag.is_available));

    //timestamp
    mxSetField(pArray, 0, fieldsMag[1], mxCreateDoubleScalar(mag.timestamp.getMilliseconds()));

    //magnetic_field_calibrated
    const mxClassID cid = cvm_traits<CV_32FC1>::CID;
    mxArray* magnetic_field_calibrated = mxCreateNumericMatrix(1, 3, cid, mxREAL);
    memcpy(mxGetPr(magnetic_field_calibrated), &(mag.magnetic_field_calibrated), 3 * sizeof(float));
    mxSetField(pArray, 0, fieldsMag[2], magnetic_field_calibrated);

    //magnetic_field_uncalibrated
    mxArray* magnetic_field_uncalibrated = mxCreateNumericMatrix(1, 3, cid, mxREAL);
    memcpy(mxGetPr(magnetic_field_uncalibrated), &(mag.magnetic_field_uncalibrated), 3 * sizeof(float));
    mxSetField(pArray, 0, fieldsMag[3], magnetic_field_uncalibrated);    

    //effective_rate
    mxSetField(pArray, 0, fieldsMag[4], mxCreateDoubleScalar(mag.effective_rate));
}

const char* fieldsBaro[] = {"is_available", "timestamp", "pressure", "relative_altitude", "effective_rate" };
void fillBaro(mxArray* pArray, sl::SensorsData::BarometerData& baro) { 
    //is_available
    mxSetField(pArray, 0, fieldsBaro[0], mxCreateDoubleScalar(baro.is_available));

    //timestamp
    mxSetField(pArray, 0, fieldsBaro[1], mxCreateDoubleScalar(baro.timestamp.getMilliseconds()));

    //pressure
    mxSetField(pArray, 0, fieldsBaro[2], mxCreateDoubleScalar(baro.pressure));
    
    //relative_altitude
    mxSetField(pArray, 0, fieldsBaro[3], mxCreateDoubleScalar(baro.relative_altitude));

    //effective_rate
    mxSetField(pArray, 0, fieldsBaro[4], mxCreateDoubleScalar(baro.effective_rate));
}

const char* fieldsIMU[] = {"is_available", "timestamp", "pose", "pose_covariance", "angular_velocity", "linear_acceleration","angular_velocity_covariance","linear_acceleration_covariance","effective_rate"};
void fillImu(mxArray* pArray, sl::SensorsData::IMUData& imu) { 
    //is_available
    mxSetField(pArray, 0, fieldsIMU[0], mxCreateDoubleScalar(imu.is_available));

    //timestamp
    mxSetField(pArray, 0, fieldsIMU[1], mxCreateDoubleScalar(imu.timestamp.getMilliseconds()));

    //pose
    cv::Mat pose = cvtMat(imu.pose);
    mxSetField(pArray, 0, fieldsIMU[2], CvMat_to_new_mxArr(pose));

    //pose_covariance
    cv::Mat pose_covariance = cvtMat(imu.pose_covariance);
    mxSetField(pArray, 0, fieldsIMU[3], CvMat_to_new_mxArr(pose_covariance));

    //angular_velocity
    const mxClassID cid = cvm_traits<CV_32FC1>::CID;
    mxArray* angular_v = mxCreateNumericMatrix(1, 3, cid, mxREAL);
    memcpy(mxGetPr(angular_v), &(imu.angular_velocity), 3 * sizeof(float));
    mxSetField(pArray, 0, fieldsIMU[4], angular_v);

    //linear_acceleration
    mxArray* linear_a = mxCreateNumericMatrix(1, 3, cid, mxREAL);
    memcpy(mxGetPr(linear_a), &(imu.linear_acceleration), 3 * sizeof(float));
    mxSetField(pArray, 0, fieldsIMU[5], linear_a);

    //angular_velocity_covariance
    cv::Mat angular_velocity_covariance = cvtMat(imu.angular_velocity_covariance);
    mxSetField(pArray, 0, fieldsIMU[6], CvMat_to_new_mxArr(angular_velocity_covariance));

    //linear_acceleration_covariance
    cv::Mat linear_acceleration_covariance = cvtMat(imu.linear_acceleration_covariance);
    mxSetField(pArray, 0, fieldsIMU[7], CvMat_to_new_mxArr(linear_acceleration_covariance));

    //effective_rate
    mxSetField(pArray, 0, fieldsIMU[8], mxCreateDoubleScalar(imu.effective_rate));
}

const char* fieldsSensors[] = {"BarometerData", "IMUData", "MagnetometerData", "TemperatureData"};
void fillSensors(mxArray* pArray, sl::SensorsData& data) {
    auto baro = mxCreateStructMatrix(1, 1, 5, fieldsBaro);
    fillBaro(baro, data.barometer);
    mxSetField(pArray, 0, fieldsSensors[0], baro);
    
    auto imu = mxCreateStructMatrix(1, 1, 9, fieldsIMU);
    fillImu(imu, data.imu);
    mxSetField(pArray, 0, fieldsSensors[1], imu);
    
    auto mag = mxCreateStructMatrix(1, 1, 5, fieldsMag);
    fillMag(mag, data.magnetometer);
    mxSetField(pArray, 0, fieldsSensors[2], mag);
    
    auto temp = mxCreateStructMatrix(1, 1, 4, fieldsTemp);
    fillTemp(temp, data.temperature);
    mxSetField(pArray, 0, fieldsSensors[3], temp);
}

template <typename T>
inline void getValue(std::string ref, std::string curr, mxArray* v_in, T& output) {
   if(ref == curr) {
        if(mxIsNumeric(v_in)) {
            int output_ID = mxGetPr(v_in)[0];
            output = static_cast<T>(output_ID);
            //mexPrintf("Set %s\n", curr.c_str());
        }
        else
            mexPrintf("Can not convert %s\n", curr.c_str());
   }
}

/* MEX entry function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	sl::ERROR_CODE err = sl::ERROR_CODE::FAILURE;
    // each sub function is referred by a string 'command'
    char command[BUFF_SIZE];
    mxGetString(prhs[0], command, BUFF_SIZE);

    if(!strcmp(command, "open")) {
        zedCam = new sl::Camera();
        sl::InitParameters initParams;
        // check if we have ONE argument
        if(checkParams(nrhs, 1)) {
            auto arg = prhs[1];
            // check if we have a parameter structure
            if(mxIsStruct(arg)) {
                // for all fields of parameter structure overwrite parameters
                for(int32_t i = 0; i < mxGetNumberOfFields(arg); i++) {
                    const char *field_name = mxGetFieldNameByNumber(arg, i);
                    mxArray *field_val = mxGetFieldByNumber(arg, 0, i);
                    if(mxIsChar(field_val) && (!strcmp(field_name, "svo_input_filename"))){
                        char string_val[BUFF_SIZE];
                        mxGetString(field_val, string_val, BUFF_SIZE);
                        initParams.input.setFromSVOFile(string_val);
                    }else{
                        getValue(field_name, "depth_mode", field_val, initParams.depth_mode);
                        getValue(field_name, "coordinate_units", field_val, initParams.coordinate_units);
                        getValue(field_name, "coordinate_system", field_val, initParams.coordinate_system);
                        getValue(field_name, "camera_resolution", field_val, initParams.camera_resolution);
                        getValue(field_name, "sdk_verbose", field_val, initParams.sdk_verbose);
                        getValue(field_name, "sdk_gpu_id", field_val, initParams.sdk_gpu_id);
                        getValue(field_name, "depth_minimum_distance", field_val, initParams.depth_minimum_distance);
                        getValue(field_name, "depth_maximum_distance", field_val, initParams.depth_maximum_distance);
                        getValue(field_name, "camera_disable_self_calib", field_val, initParams.camera_disable_self_calib);
                        getValue(field_name, "camera_image_flip", field_val, initParams.camera_image_flip);
                        getValue(field_name, "camera_fps", field_val, initParams.camera_fps);
                        getValue(field_name, "svo_real_time_mode", field_val, initParams.svo_real_time_mode);
                        getValue(field_name, "depth_stabilization", field_val, initParams.depth_stabilization);
                        getValue(field_name, "enable_right_side_measure", field_val, initParams.enable_right_side_measure);  
                    }                  
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
            sl::RuntimeParameters param;
			// check if we have a parameter structure
			if ((nrhs > 1) && (mxIsStruct(prhs[1]))) {
				// for all fields of parameter structure overwrite parameters
				for (int32_t i = 0; i < mxGetNumberOfFields(prhs[1]); i++) {
					const char *field_name = mxGetFieldNameByNumber(prhs[1], i);
					mxArray *field_val = mxGetFieldByNumber(prhs[1], 0, i);
                    
                    getValue(field_name, "sensing_mode", field_val, param.sensing_mode);
                    getValue(field_name, "enable_depth", field_val, param.enable_depth);
                    getValue(field_name, "measure3D_reference_frame", field_val, param.measure3D_reference_frame);
				}
			}
			err = zedCam->grab(param);
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

            sl::Mat tmp;
            err = zedCam->retrieveImage(tmp, view, sl::MEM::CPU, sl::Resolution(width, height));
            int cv_type = (view == sl::VIEW::LEFT_GRAY || view == sl::VIEW::RIGHT_GRAY || view == sl::VIEW::LEFT_UNRECTIFIED_GRAY || view == sl::VIEW::RIGHT_UNRECTIFIED_GRAY) ? CV_8UC1 : CV_8UC4;
            cv::Mat cv_tmp(tmp.getHeight(), tmp.getWidth(), cv_type, tmp.getPtr<sl::uchar1>());
                        
            if(cv_type == CV_8UC4){
                cv::Mat image_rgb;
                cv::cvtColor(cv_tmp, image_rgb, cv::COLOR_RGBA2RGB);
                plhs[0] = CvMat_to_new_mxArr(image_rgb);
            }else{
                cv::Mat image_gray;
                cv::cvtColor(cv_tmp, image_gray, cv::COLOR_GRAY2RGB);
                plhs[0] = CvMat_to_new_mxArr(image_gray);
            }
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
            
            sl::Mat tmp;
            err = zedCam->retrieveMeasure(tmp, measure, sl::MEM::CPU, sl::Resolution(width, height));
            int cv_type = ((measure >= sl::MEASURE::XYZ &&  measure <= sl::MEASURE::NORMALS) || (measure >= sl::MEASURE::XYZ_RIGHT &&  measure < sl::MEASURE::LAST)) ? CV_32FC4 : CV_32FC1;
            cv::Mat measure_mat(tmp.getHeight(), tmp.getWidth(), cv_type, tmp.getPtr<sl::uchar1>());
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
            int time_ref = mxGetPr(prhs[1])[0];
            double val = zedCam->getTimestamp(static_cast<sl::TIME_REFERENCE>(time_ref)).getMilliseconds();
            plhs[0] = mxCreateDoubleScalar(val);
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
            sl::PositionalTrackingParameters  param;
            // check if we have a parameter structure
            if((nrhs > 1) && (mxIsStruct(prhs[1]))) {
				sl::Transform tranform;
                // for all fields of parameter structure overwrite parameters
                for(int32_t i = 0; i < mxGetNumberOfFields(prhs[1]); i++) {
                    const char *field_name = mxGetFieldNameByNumber(prhs[1], i);
                    mxArray *field_val = mxGetFieldByNumber(prhs[1], 0, i);
                    if(mxIsChar(field_val) && (!strcmp(field_name, "area_file_path"))){
                        char string_val[BUFF_SIZE];
                        mxGetString(field_val, string_val, BUFF_SIZE);
                        param.area_file_path.set(string_val);
                    } else {                        
                        getValue(field_name, "enable_area_memory", field_val, param.enable_area_memory);
                        getValue(field_name, "enable_imu_fusion", field_val, param.enable_imu_fusion);
                        getValue(field_name, "enable_pose_smoothing", field_val, param.enable_pose_smoothing);
                        getValue(field_name, "set_as_static", field_val, param.set_as_static);
                        getValue(field_name, "set_floor_as_origin", field_val, param.set_floor_as_origin);
                        if(!strcmp(field_name, "initial_world_transform")) {
                            for(int col = 0; col < mxGetN(field_val); col++)
                                for(int row = 0; row < mxGetM(field_val); row++)
                                    tranform(row, col) = (mxGetPr(field_val))[row + col*mxGetM(field_val)];
                            param.initial_world_transform = tranform;
                        }
                    }
                }
            }
			err = zedCam->enablePositionalTracking(param);            
        }
		// we return the string associated with the error
		plhs[0] = mxCreateString(sl::toString(err).c_str());
    }

    else if(!strcmp(command, "getPosition")) {
        if(checkZED() && checkParams(nrhs, 0)) {
            sl::REFERENCE_FRAME ref_f = sl::REFERENCE_FRAME::WORLD;
            if (nrhs == 2) {
                int val = mxGetPr(prhs[1])[0];
                ref_f = static_cast<sl::REFERENCE_FRAME>(val);
            }

            sl::Pose path;          
            zedCam->getPosition(path,ref_f);
            auto position = cvtMat(path.pose_data);
            plhs[0] = CvMat_to_new_mxArr(position);
            err = sl::ERROR_CODE::SUCCESS;
        }
    }

    else if(!strcmp(command, "getSensorsData")) {
        if(checkZED() && checkParams(nrhs, 1)) {
            sl::SensorsData sensors;
            sl::TIME_REFERENCE t_ref = sl::TIME_REFERENCE::IMAGE;
            if (nrhs == 2) {
                int val = mxGetPr(prhs[1])[0];
                t_ref = static_cast<sl::TIME_REFERENCE>(val);                
            }
            err = zedCam->getSensorsData(sensors, t_ref);
            plhs[0] = mxCreateStructMatrix(1, 1, 4, fieldsSensors);
            fillSensors(plhs[0], sensors);
        }
    }
   
    else if(!strcmp(command, "disablePositionalTracking")) {
        if(checkZED()) {
            if((nrhs > 1) && (mxIsChar(prhs[1]))) {
                char area_file_path[BUFF_SIZE];
                mxGetString(prhs[1], area_file_path, BUFF_SIZE);
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
             if(nrhs == 2) {
                sl::RecordingParameters param;
                if(mxIsStruct(prhs[1])){
                    for(int32_t i = 0; i < mxGetNumberOfFields(prhs[1]); i++) {
                        const char *field_name = mxGetFieldNameByNumber(prhs[1], i);
                        mxArray *field_val = mxGetFieldByNumber(prhs[1], 0, i);
                        if(mxIsChar(field_val) && (!strcmp(field_name, "video_filename"))) {
                            char svo_path[BUFF_SIZE];
                            mxGetString(field_val, svo_path, BUFF_SIZE);
                            param.video_filename.set(svo_path);
                        }else{
                            getValue(field_name, "bitrate", field_val, param.bitrate);
                            getValue(field_name, "compression_mode", field_val, param.compression_mode);
                            getValue(field_name, "target_framerate", field_val, param.target_framerate);
                            getValue(field_name, "transcode_streaming_input", field_val, param.transcode_streaming_input);
                        }
                    }
                }else if(mxIsChar(prhs[1])) {
                    char svo_path[BUFF_SIZE];
                    mxGetString(prhs[1], svo_path, BUFF_SIZE);
                    param.video_filename.set(svo_path);
                }
                err = zedCam->enableRecording(param);
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
            if(nrhs == 2) {
                auto arg = prhs[1];
                // check if we have a parameter structure            
                if(mxIsStruct(arg)) {
                    // for all fields of parameter structure overwrite parameters
                    for(int32_t i = 0; i < mxGetNumberOfFields(arg); i++) {
                        const char *field_name = mxGetFieldNameByNumber(arg, i);
                        mxArray *field_val = mxGetFieldByNumber(arg, 0, i);
                        getValue(field_name, "detection_model", field_val, param.detection_model);
                        getValue(field_name, "enable_body_fitting", field_val, param.enable_body_fitting);
                        getValue(field_name, "enable_mask_output", field_val, param.enable_mask_output);
                        getValue(field_name, "enable_tracking", field_val, param.enable_tracking);
                    }
                }
            }
            err = zedCam->enableObjectDetection(param);
        }
    }

    else if (!strcmp(command, "retrieveObjects")) {
        if (checkZED()) {
            sl::ObjectDetectionRuntimeParameters param;
             if(nrhs == 2) {
                auto arg = prhs[1];
                // check if we have a parameter structure            
                if(mxIsStruct(arg)) {
                    // for all fields of parameter structure overwrite parameters
                    for(int32_t i = 0; i < mxGetNumberOfFields(arg); i++) {
                        const char *field_name = mxGetFieldNameByNumber(arg, i);
                        auto field_val = mxGetFieldByNumber(arg, 0, i);
                        getValue(field_name, "detection_confidence_threshold", field_val, param.detection_confidence_threshold);
                    }
                } else
                    param.detection_confidence_threshold = static_cast<float>(mxGetPr(arg)[0]);                
            }           

            sl::Objects objs;
            err = zedCam->retrieveObjects(objs, param);
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
            plhs[0] = mxCreateDoubleScalar(val);
            err = sl::ERROR_CODE::SUCCESS;
        }
    }

    else if (!strcmp(command, "savePointCloudAs")) {
        if (checkParams(nrhs, 2) || checkParams(nrhs, 3)) {
            double *ptr_ = mxGetPr(prhs[1]);
            int format = ptr_[0];

            char save_path[BUFF_SIZE];
            if (mxIsChar(prhs[2])) {
                mxGetString(prhs[2], save_path, BUFF_SIZE);
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

            char save_path[BUFF_SIZE];
            if (mxIsChar(prhs[2])) {
                mxGetString(prhs[2], save_path, BUFF_SIZE);
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
    else if (!strcmp(command, "enableSpatialMapping")) {
        if (checkZED()) {
            if (nrhs == 2) {
                auto arg = prhs[1];
                // check if we have a parameter structure
                if (mxIsStruct(arg)) {
                    // for all fields of parameter structure overwrite parameters
                    for (int32_t i = 0; i < mxGetNumberOfFields(arg); i++) {
                        const char* field_name = mxGetFieldNameByNumber(arg, i);
                        mxArray* field_val = mxGetFieldByNumber(arg, 0, i);
                        getValue(field_name, "map_type", field_val, sp_p.map_type);
                        getValue(field_name, "max_memory_usage", field_val, sp_p.max_memory_usage);
                        getValue(field_name, "range_meter", field_val, sp_p.range_meter);
                        getValue(field_name, "resolution_meter", field_val, sp_p.resolution_meter);
                    }
                }
            }
            err = zedCam->enableSpatialMapping(sp_p);
        }
    }
    else if (!strcmp(command, "extractWholeSpatialMap")) {
        if (checkZED()) {
            if (sp_p.map_type == sl::SpatialMappingParameters::SPATIAL_MAP_TYPE::MESH) {
                sl::Mesh mesh;
                err = zedCam->extractWholeSpatialMap(mesh);
                if (err == sl::ERROR_CODE::SUCCESS) {
                    cv::Mat vert(mesh.vertices.size(), 3, CV_32FC1);
                    for (int i = 0; i < vert.rows; i++) {
                        auto v = mesh.vertices[i];
                        vert.at<float>(i, 0) = v.x;
                        vert.at<float>(i, 1) = v.y;
                        vert.at<float>(i, 2) = v.z;
                    }
                    cv::Mat faces(mesh.triangles.size(), 3, CV_16UC1);
                    for (int i = 0; i < faces.rows; i++) {
                        auto f = mesh.triangles[i];
                        faces.at<ushort>(i, 0) = f.x;
                        faces.at<ushort>(i, 1) = f.y;
                        faces.at<ushort>(i, 2) = f.z;
                    }
                    plhs[0] = CvMat_to_new_mxArr(vert);
                    plhs[1] = CvMat_to_new_mxArr(faces);
                }
            }
            else{
                sl::FusedPointCloud fpc;
                err = zedCam->extractWholeSpatialMap(fpc);
                if (err == sl::ERROR_CODE::SUCCESS) {
                    cv::Mat vert(fpc.vertices.size(), 3, CV_32FC1);
                    cv::Mat clrs(fpc.vertices.size(), 1, CV_8UC3);
                    for (int i = 0; i < vert.rows; i++) {
                        auto v = fpc.vertices[i];
                        vert.at<float>(i, 0) = v.x;
                        vert.at<float>(i, 1) = v.y;
                        vert.at<float>(i, 2) = v.z;
                        // depack the color
                        uint32_t color_uint = *(uint32_t*)&v.w;
                        auto color_uchar = (uchar*)&color_uint;
                        clrs.at<cv::Vec3b>(i, 0) = cv::Vec3b(color_uchar[0], color_uchar[1], color_uchar[2]);
                    }
                    plhs[0] = CvMat_to_new_mxArr(vert);
                    plhs[1] = CvMat_to_new_mxArr(clrs);
                }
            }
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
