%% SAMUEL R PILON
% AE 426 - Project 3 
% For personal use of SAMUEL PILON only. Not to be distributed
% Collaborators: Samuel Horine, Ashley Tirado, Sage Webster, Connor
% Bramhall, Saaim Rahman
clc
close
clear all; 

%% SIMULATION
%-------------------------------------------------------------------------
%{
1) Initial conditions: the satellite is in the equatorial orbit of 
400 km altitude with the given attitude configuration and always 
pointing at the Earth center.
a.)Perturbation about â1 axis:
b.)Perturbation about â1 axis (change quaternion): 
   
%---------------------------

2) Perturbation about â1 axis with non-zero initial spin about â2
Initial conditions: the satellite is in the equatorial orbit of 400 
km altitude with the given attitude configuration and angular velocity 
of the (orbital rate + spin rate).

%---------------------------

3) Perturbation for all three axes with zero initial spin.
Initial conditions: the satellite is in the polar orbit of 
400 km altitude with the arbitrary attitude configuration 
(explain your initial attitude) and angular velocity of the orbital rate.
%}
%-------------------------------------------------------------------------
%% SOLVER OPTIONS
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
w1 = 0.01; % Initial precession rate
w2 = 0.01; % Initial nutation rate
w3 = 0.04; % Initial spin rate
q1 = 0.0; % Initial q1 (x)
q2 = 0.0; % Initial q2 (y) 
q3 = 0.707; % Initial q3 (z)
q4 = 0.707; % Initial q4 (axis)
t= 0:10:1500; % T-span

% Solve system of ODEs using ode45
[t,x] = ode45(@OdeFun,[t],[w1,w2,w3,q1,q2,q3,q4],options);

%% PLOT THE RESULTS
figure
tiledlayout(3,3)
sgtitle('Simulation Results for Perturbation for all three axes with initial spin (Angled Case)', 'FontSize', 16, 'FontWeight', 'bold')
nexttile
plot(t,x(:,1)*180/pi); 
grid minor
title('\boldmath{$\omega_1$} (Precession Rate)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time (seconds)');
ylabel('Angular Velocity (deg/s)');

nexttile
plot(t,x(:,2)*180/pi,'b'); 
grid minor
title('\boldmath{$\omega_2$} (Nutation Rate)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time (seconds)');
ylabel('Angular Velocity (deg/s)');

nexttile
plot(t,x(:,3)*180/pi,'r'); 
grid minor
title('\boldmath{$\omega_3$} (Spin Rate)', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('Time (seconds)');
ylabel('Angular Velocity (deg/s)');

nexttile
plot(t,x(:,4)); 
grid minor
title('First Component of Quaternion');
xlabel('Time (seconds)');

nexttile
plot(t,x(:,5),'b'); 
grid minor
title('Second Component of Quaternion');
xlabel('Time (seconds)');

nexttile
plot(t,x(:,6),'r'); 
grid minor
title('Third Component of Quaternion');
xlabel('Time (seconds)');

nexttile
plot(t,x(:,7),'r'); 
grid minor
title('Fourth Component of Quaternion');
xlabel('Time (seconds)');

%% ODE 45 FUNCTION 
function quateqn=OdeFun(t,x)
    A=500;   
    C=1000;      
    spin=0; %Spin Rate rad/s
    capOmega = 0.00113; %Omega rad/s (Calculate Rate by altitude and orbital period) 
%_______________________________
    % ODE KEY 
    %x(1) = omegadot_1
    %x(2) = omegadot_2
    %x(3) = omegadot_3
    %x(4) = quaternion1
    %x(5) = quaternion2
    %x(6) = quaternion3
    %x(7) = quaternion4
%_______________________________
     
    quateqn=[(-spin*x(2))+(1-C/A)*(x(2)*x(3))*-12*capOmega^2*(x(4)*x(5)-x(6)*x(7))*((x(6)*x(4)+x(5)*x(7)));... %w1dot             
        (spin*x(1))-(1-C/A)*(x(1)*x(3))-6*capOmega^2*(x(6)*x(4)+x(5)*x(7))*(1-2*x(5)^2-2*x(6)^2);... %w2dot
        0;... %w3dot
        (x(5)*(x(3)-spin+capOmega)-x(6)*x(2)+x(7)*x(1))/2;... %time rate of change for quaternion 1
        (x(6)*x(1)+x(7)*x(2)-x(4)*(x(3)-spin+capOmega))/2;... %time rate of change for quaternion 2
        (x(7)*(x(3)-spin-capOmega)+x(4)*x(2)-x(5)*x(1))/2;...  %time rate of change for quaternion 3
        ((-x(4)*x(1)-x(5)*x(2)-x(6)*(x(3)-spin-capOmega))/2)... %time rate of change for quaternion 4
        ];
    
end