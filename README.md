# Stereolabs ZED - Matlab Integration

This sample shows how to use the ZED SDK functionalities  within Matlab.

## Getting started

- First, download the latest version of the ZED SDK on [stereolabs.com](https://www.stereolabs.com/developers)
- For more information, read the ZED [Documentation](https://www.stereolabs.com/docs/) and [API documentation](https://www.stereolabs.com/docs/api/python/) or our [Community page](https://community.stereolabs.com)

### Prerequisites

- Windows 10 64bits or later, Ubuntu 16.04 or higher.
- Matlab with MEX compiler installed
- [ZED SDK](https://www.stereolabs.com/developers/release/) and its dependency ([CUDA](https://developer.nvidia.com/cuda-downloads))


## Build the program

#### Build for Windows

For detailed installation instructions, check out our [documentation](https://www.stereolabs.com/docs/matlab#using-the-zed-sdk).

#### Build for Linux

Download the sample and execute the following command in a terminal:

    export MATLAB_ROOT=/usr/local/MATLAB/R2012b # Put your actual Matlab path here
    mkdir build
    cd build
    cmake ../src
    make
    make install

## Run the program
In the Matlab directory, open the file `ZED_DepthSensing.m` with Matlab and press run. Press any key to exit the program.

## Features

This sample is split into 6 scripts, each of them shows a specific feature of the ZED SDK in Matlab.

1. `ZED_DepthSensing.m` :

	In this part, we first initialize the ZED Camera and then retrieve both stereo images (Left and Right) as well as the depth data (32-bits float buffer) and the depth image normalized to 8-bits (grayscale). Then, with Matlab, we compute the depth histogram from the depth data values and display it.

	<img src="ZED_Camera.png" alt="Drawing" style="width: 600px;"/>
	
2. `ZED_Tracking.m` :

	This part shows how to get the positional tracking information. First, we initialize the ZED camera and the tracking, then we display the current left image as well as all the positions of the camera.

	<img src="ZED_Tracking.png" alt="Drawing" style="width: 600px;"/>

3. `ZED_PointCloud.m` :

	This part shows how to initialize the ZED camera and retrieve for each frame the current point cloud (X,Y,Z) of the scene. Then we use the mesh function of Matlab to display it.

	<img src="ZED_PointCloud.png" alt="Drawing" style="width: 600px;"/>

4. `ZED_Recording.m` :

	This script shows how to record a video sequence into 'svo' format which can be playback with the ZED SDK.
	
5. `ZED_ObjectDetection.m` :

	This script use the ability of the ZED2 to detect objects, a video window shows the results.
	
6. `ZED_CameraControl.m` :

	This script shows how to control the video settings of the ZED camera.
	
### Limitations
This sample application is not designed to run in real time


### Troubleshooting
On Linux, If you get the following error :

     /usr/local/MATLAB/R2014a/bin/glnxa64/../../sys/os/glnxa64/libstdc++.so.6: version `GLIBCXX_3.4.21` not found

Launch matlab with `LD_PRELOAD`:

    LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.21:/usr/lib/x86_64-linux-gnu/libfreetype.so.6 matlab&

You can also create an alias in a bash terminal (just once) to simplify the launch:

    echo "alias matlab=LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.21:/usr/lib/x86_64-linux-gnu/libfreetype.so.6 matlab&" >> ~/.bashrc

Then just launch the application by typing:

    matlab
    
## Support
If you need assistance go to our Community site at https://community.stereolabs.com/
