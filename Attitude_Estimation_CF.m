%% Atitude Estimation Complementary Filter
clc
close all
clear all

% Load Data
load('slow_v4.mat');

input = XS2(:,2:10);
out_q = XS2(:,11:14);
% The nominal average value at Earth's surface (standard gravity)
% Thompson, Ambler, and Barry N. Taylor. "Use of the International System of Units (SI)." (2008).
g = 9.80665;

for i = 1:length(out_q)
    q0 = out_q(i,1);
    q1 = out_q(i,2);
    q2 = out_q(i,3);
    q3 = out_q(i,4);
    
    [phi, theta, psi] = quat_2_eul(q0,q1,q2,q3);
    phi_ref(i,1) = phi;
    theta_ref(i,1)=theta;
    psi_ref(i,1)=psi;
end

    dt      =   1/100; % Sassari Dataset Sampling Rate = 100 Hz
    t       =   0:dt:(length(input) - 1)*dt;
    
    acc_x   =   input(:,1);
    acc_y   =   input(:,2);
    acc_z   =   input(:,3);
    gyro_x  =   input(:,4);
    gyro_y  =   input(:,5);
    gyro_z  =   input(:,6);


%% Calculate Phi and Theta using Accelerometer and Gyroscope Measurements  
    phi_gyro    =   zeros(length(input),1);
    theta_gyro  =   zeros(length(input),1);
    
    % Accelerometer
    for i = 1:length(input)
        % First Equ.
        phi1_acc(i,1)	=	atan2(acc_y(i),sqrt(acc_x(i)^2 + acc_z(i)^2));
        theta1_acc(i,1)	=	atan2(acc_x(i),sqrt(acc_y(i)^2 + acc_z(i)^2));
        
        % Second Equ.
        phi2_acc(i,1)   =   atan2(acc_y(i),acc_z(i));
        theta2_acc(i,1) =   atan2(acc_x(i),(acc_y(i)*sin(phi2_acc(i))+acc_z(i)*cos(phi2_acc(i))));
        % Third Equ.
        theta3_acc(i,1) =   atan2(acc_x(i),g);
    end
    % Gyroscope
     for i = 2:length(input)
        p   =   gyro_x(i);
        q   =   gyro_y(i);
        r   =   gyro_z(i);

        phi     =	phi_gyro(i-1);
        theta   =   theta_gyro(i-1);

        phi_gyro(i,1)   =     phi     +   dt*(p + (sin(phi)*tan(theta)*q) + (cos(phi)*tan(theta)*r));
        theta_gyro(i,1) =     theta   +   dt*((cos(phi)*q) - (sin(phi)*r));
     end
     
    
    figure(1)
    plot(t,rad2deg(phi1_acc),t,rad2deg(phi2_acc),t,rad2deg(phi_ref),LineWidth = 1.2)
    legend('Phi1','Phi2', 'Ref')
    title('Phi Acc')
    RMSE_phi_est_1 = sqrt(mean( (rad2deg(phi1_acc)-rad2deg(phi_ref)).^2 ));
    RMSE_phi_est_2 = sqrt(mean( (rad2deg(phi2_acc)-rad2deg(phi_ref)).^2 ));
    
    figure(2)
    plot(t,rad2deg(theta1_acc),t,rad2deg(theta2_acc),t,rad2deg(theta3_acc),t,-rad2deg(theta_ref),LineWidth = 1.2)
    legend('Theta1','Theta2','Theta3', 'Ref')
    title('Theta Acc')
    RMSE_theta_est1 = sqrt(mean( (rad2deg(theta1_acc)+rad2deg(theta_ref)).^2 ));
    RMSE_theta_est2 = sqrt(mean( (rad2deg(theta2_acc)+rad2deg(theta_ref)).^2 ));
    RMSE_theta_est3 = sqrt(mean( (rad2deg(theta3_acc)+rad2deg(theta_ref)).^2 ));
    
    
    figure(3)
    plot(t,phi_gyro,t,phi_ref)
    legend('Gyro','Ref')
    title('Phi Gyro')
    RMSE_phi_gyro_ = sqrt(mean( (rad2deg(phi_gyro)-rad2deg(phi_ref)).^2 ));
    
    figure(4)
    plot(t,theta_gyro,t,+theta_ref)
    legend('Gyro','Ref')
    title('Theta Gyro')
    RMSE_theta_gyro_ = sqrt(mean( (rad2deg(theta_gyro)-rad2deg(theta_ref)).^2 ));
    
    %% Complementary Filter
    phi_CF1   = zeros(length(input),1);
    phi_CF2   = zeros(length(input),1);
    theta_CF1 = zeros(length(input),1);
    theta_CF2 = zeros(length(input),1);
    theta_CF3 = zeros(length(input),1);
    alpha = 0.3;
    for i = 2: length (input)
  
        phi_CF1(i,1)	= ((1-alpha) * phi_gyro(i))     + (alpha * phi1_acc(i));
        phi_CF2(i,1)	= ((1-alpha) * phi_gyro(i))     + (alpha * phi2_acc(i));
        
        theta_CF1(i,1)	= ((1-alpha) * theta_gyro(i))	+ (alpha * theta1_acc(i));
        theta_CF2(i,1)	= ((1-alpha) * theta_gyro(i))	+ (alpha * theta2_acc(i));
        theta_CF3(i,1)	= ((1-alpha) * theta_gyro(i))	+ (alpha * theta3_acc(i));
    end
    
    figure(5)
    plot(t,rad2deg(phi_CF1),t,rad2deg(phi_CF2),t,rad2deg(phi_ref))
    legend('CF1','CF2','Ref')
     title('Phi CF')
    RMSE_phi_CF1 = sqrt(mean( (rad2deg(rad2deg(phi_CF1))-rad2deg(phi_ref)).^2 ))
    RMSE_phi_CF2 = sqrt(mean( (rad2deg(rad2deg(phi_CF2))-rad2deg(phi_ref)).^2 ))
    
    figure(6)
    plot(t,rad2deg(theta_CF1),t,rad2deg(theta_CF2),t,rad2deg(theta_CF3),t,-rad2deg(theta_ref))
    legend('CF1','CF2','CF3','Ref')
    title('Theta CF')
    RMSE_theta_CF1 = sqrt(mean( (rad2deg(theta_CF1)+rad2deg(theta_ref)).^2 ))
    RMSE_theta_CF2 = sqrt(mean( (rad2deg(theta_CF2)+rad2deg(theta_ref)).^2 ))
    RMSE_theta_CF3 = sqrt(mean( (rad2deg(theta_CF3)+rad2deg(theta_ref)).^2 ))
    
    