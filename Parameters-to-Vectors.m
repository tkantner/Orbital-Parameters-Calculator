%Housekeeping
%clear all;
clc;
close all;

%Constants
a = 3.1184E7/1000; %km
e = .4658; 
i = 62.52; %deg
OMEGA = 73.74; %deg
omega = 22.08; %deg
theta = 353.6; %deg
mu = 398600; %km^3/s^2
h = sqrt(a.*mu.*(1-e.^2)); %km^2/s

%e = 0.971; 
%i = 162.79; %deg
%OMEGA = 183.27; %deg
%omega = 358.36; %deg
%theta = 190; %deg
%mu = 398600; %km^3/s^2
%h = 11850.418; %km^2/s

%Create Matrices
A = zeros(3,1); %r transformation
A(1) = cosd(theta);
A(2) = sind(theta);

B = zeros(3,1); %v transformation
B(1) = -sind(theta);
B(2) = e+cosd(theta);

%Q
Q = zeros(3,3);
Q(1,1) = -sind(OMEGA)*cosd(i)*sind(omega)+cosd(OMEGA)*cosd(omega);
Q(1,2) = -sind(OMEGA)*cosd(i)*cosd(omega)-cosd(OMEGA)*sind(omega);
Q(1,3) = sind(OMEGA)*sind(i);
Q(2,1) = cosd(OMEGA)*cosd(i)*sind(omega)+sind(OMEGA)*cosd(omega);
Q(2,2) = cosd(OMEGA)*cosd(i)*cosd(omega)-sind(OMEGA)*sind(omega);
Q(2,3) = -cosd(OMEGA)*sind(i);
Q(3,1) = sind(i)*sind(omega);
Q(3,2) = sind(i)*cosd(omega);
Q(3,3) = cosd(i);

%Position perifocal
rpf = (h.^2./mu).*(1./(1+(e.*cosd(theta))))*A;

%Velocity perifocal
vpf = (mu./h)*B;

%Position ECI
reci = Q*rpf; %ANS, km

%Velocity ECI
veci = Q*vpf; %ANS, km/s

%Print Results

fprintf('Rx = %.4f km\n', reci(1));
fprintf('Ry = %.4f km\n', reci(2));
fprintf('Rz = %.4f km\n', reci(3));
fprintf('Vx = %.4f km/s\n', veci(1));
fprintf('Vy = %.4f km/s\n', veci(2));
fprintf('Vz = %.4f km/s\n', veci(3));
