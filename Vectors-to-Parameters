
%Housekeeping
clear all;
clc;
close all;

%Constants
r = [2500, 16000, 4000]; %km
v = [-3, -1, 5]; %km/s
mu = 389600; %km^3/s^2

%Calculations

%Radial Velocity
vr = dot(r,v)/norm(r);

%Angular Velocity
h = cross(r,v); %km^2/s

%Inclination
i = acosd(h(3)/norm(h)); %ANS

N = cross([0, 0, 1], h);

%Big Omega
if N(2) >= 0    %ANS
    OMEGA = acosd(N(1)/norm(N));
else
    OMEGA = 360 - acosd(N(1)/norm(N));
end

%Eccentricity
e = (1/mu)*((norm(v).^2-(mu./norm(r))).*r -(norm(r).*vr.*v)); %ANS

%Semi-Major axis
a = (norm(h).^2/mu)*(1./(1-norm(e).^2));  %ANS

%Little Omega
if e(3) >= 0  %ANS
    omega = acosd(dot(N./norm(N),e./norm(e)));
else
    omega = 360 - acosd(dot(N./norm(N),e./norm(e)));
end

%True Anomaly
if vr >= 0 %ANS
    theta = acosd(dot(e./norm(e),r./norm(r)));
else
    theta = 360 - acosd(dot(e./norm(e),r./norm(r)));
end

%Print Answers
fprintf('a = %.3f km\n', a);  %semi-major
fprintf('h = %.3f km^2/s\n', norm(h));  %ang momentum
fprintf('e = %.3f\n', norm(e));  %eccentricity
fprintf('i = %.2f degrees\n', i); %inclination
fprintf('(Big) Omega = %.2f degrees\n', OMEGA); %RAAN
fprintf('(Little) omega = %.2f degrees\n', omega); %arg of periapsis
fprintf('Theta = %.2f degrees\n', theta);  %true anomaly
