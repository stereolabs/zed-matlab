# zed-matlab

**This sample is designed to work with the ZED stereo camera only and requires the ZED SDK. For more information: https://www.stereolabs.com**

It demonstrates how to use most of the ZED SDK functionalities with Matlab.

**Warning :**
 - This sample is not designed to operate in real time

This sample displays the both left and right images of the ZED as well as the normalized depth map.
This sample also retrieve the depth information and then compute thru Matlab the depth histogram.

##Build the program
To get more detailed instructions (especially for Windows) check out our blog.

Open a terminal in zed-matlab directory and execute the following command:

    $ export MATLAB_ROOT=/usr/local/MATLAB/R2012b # Change this with your actual Matlab path
    $ mkdir build
    $ cd build
    $ cmake ../src
    $ make
    $ make install


##Run the program
In the matlab directory, open the file ZED_Camera.m with Matlab and press run.


**Quit :**
Press any key to exit the program.
