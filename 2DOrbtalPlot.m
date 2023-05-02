clc;
clear all;
G = 6.673e-11; %Gravitational constant
M = 5.98e24; %mass of earth in (kg)
ra = 400000; %orbit distance in (m)
r = 6.37e6 + ra; %total radius of orbit in (m)
m = 1.5; % mass satelite (kg)
a = (G*M)/(r^2); % check for acceleration
v_orb = sqrt((G*M)/r); % orbital velocity (m/s)
T = sqrt(((4*(pi^2))*r^3)/(G*M)); % period (s)
%lets graph it for the period in steps of 150 to reduce the computation
steps = T/35;
simt = -steps;
for i = 1:1:36
  simt = simt+steps;
  t(i) = simt;
  v_o(i) = sqrt((G*M)/r);% array of velocity corresponding to time in steps of 150 .. should not change and be 1x34 - but ERROR
  rn(i) = 6.37e6 + ra; % array of radius corresponding to time in steps of 150 .. should not change and be 1x34 - but ERROR
  %disp(t2) % goes from 1 to 5101 in steps of 150
  
  
  end
  
figure
plot(t, v_o); xlabel('Simulation Time', 'FontSize', 12);
ylabel('Velocity', 'FontSize', 12);
grid;  %%% expecting a graph of straight line %%%
figure
plot(t, rn); xlabel('Simulation Time', 'FontSize', 12);
ylabel('Radius', 'FontSize', 12);
grid;  %%% expecting a graph of straight line %%%
%%% take the points of radius and time convert to polar coordinates and plot %%%
%%% expecting a circle as satellite moves around the orbit of constant radius in a sice of time of the period%%%
%[theta,rho] = cart2pol(t,rn);
%th = atan2(rn,t);
th = 2*pi * t / T;
rho = sqrt((t.^2)+(rn.^2));
%disp(th);
figure
polar(th,rho)
title('Orbit')