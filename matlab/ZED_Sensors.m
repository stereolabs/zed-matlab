clc;
disp('========= ZED SDK PLUGIN =========');
disp('-- Sensors --');
close all;
clear mex; clear functions; clear all;

% initial parameter structure, the same as sl::InitParameters
% values as enum number, defines in : sl/defines.hpp
% or from https://www.stereolabs.com/docs/api/structsl_1_1InitParameters.html

InitParameters.camera_resolution = 0; %HD2K
InitParameters.depth_mode = 0; % No depth required
result = mexZED('open', InitParameters);

if(strcmp(result,'SUCCESS')) % the Camera is open 
	% basic informations
    camInfo = mexZED('getCameraInformation');    
    if(camInfo.model ~= 0) % ZED has no other sensors
        
        % Create Figure and wait for keyboard interruption to quit
        f = figure('name','ZED SDK: Sensors','NumberTitle','off', 'keypressfcn',@(obj,evt) 0);
        ha1 = axes('Position',[0.1,0.6,0.89,0.35]);
        hold on;
        title('Linear Acceleration')
        xlabel('frame')
        ylabel('Acceleration')
        ha2 = axes('Position',[0.1,0.05,0.89,0.35]);
        hold on;
        title('Angular Velocity')
        xlabel('frame')
        ylabel('Velocity')
        sample =0;
        key = 1;
        % loop over frames, till Esc is pressed
        while (key ~= 27)
            
            sensors_data = mexZED('getSensorsData', 1); % ask CURRENT sensors data

            % Baro
            if(sensors_data.BarometerData.is_available)
                baro_pressure = sensors_data.BarometerData.pressure
                baro_rate = sensors_data.BarometerData.effective_rate
            end
            
            % IMU
            if(sensors_data.IMUData.is_available)
                imu_pose = sensors_data.IMUData.pose
                imu_angular_v = sensors_data.IMUData.angular_velocity
                imu_linear_a = sensors_data.IMUData.linear_acceleration
                imu_rate = sensors_data.IMUData.effective_rate
                % Display
                axes(ha1);
                plot(sample, imu_linear_a(1), 'r+')
                plot(sample, imu_linear_a(2), 'b+')
                plot(sample, imu_linear_a(3), 'c+')
                axes(ha2);
                plot(sample, imu_angular_v(1), 'r+')
                plot(sample, imu_angular_v(2), 'b+')
                plot(sample, imu_angular_v(3), 'c+')
            end
            
            % Magnetometer
            if(sensors_data.MagnetometerData.is_available)
                mag_field_uncalib = sensors_data.MagnetometerData.magnetic_field_uncalibrated
                mag_field = sensors_data.MagnetometerData.magnetic_field_calibrated
                mag_rate = sensors_data.MagnetometerData.effective_rate
            end

            % Temperature
            temp_imu = sensors_data.TemperatureData.IMU
            temp_baro = sensors_data.TemperatureData.BAROMETER
            temp_left = sensors_data.TemperatureData.ONBOARD_LEFT
            temp_right = sensors_data.TemperatureData.ONBOARD_RIGHT
            
            sample = sample+1;
            drawnow;  %this checks for interrupts  
            key = uint8(get(f,'CurrentCharacter'));
            if(~length(key))
                key=0;
            end
        end
        close(f)
    end
    mexZED('close')
    disp('========= END =========');
    clear mex;
end
