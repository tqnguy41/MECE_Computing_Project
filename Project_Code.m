%project MECE 
%Toan Nguyen

%housekeepin
clear all
close all
clc

%known variable
ax = 0;
ay = 0;
bx = 2*pi;
by = 2*pi;
y = 0;
%known function
f_functtion = @(x) sin*(pi*(x-ay)/(bx-ax))*cos(pi/2*(2*(y-ay)/by-ay)+1);
gb_function = @(x) (bx-x)^2*cos(pi*x)/bx;
fb_function = @(x) x*(bx-x)^2;

%matrix size
i = 2*pi;

%step size
n = i/5;
m  = i/5;

%total count
sizes = i/n;

%matrix
u_maxtrix = zeros(sizes,sizes);

counter = 1;
%initlize the first knows
for x = n:n:i
    %since step size is all the same
    y = x;
    %top boundary
    u_maxtrix(1,counter) = fb_function(x);
    %bottom bounday
    u_maxtrix(sizes,counter) = gb_function(x);
    %left boundary
    if (x ~= n) && (x ~=i )
        u_maxtrix(sizes-counter+1,1) = gb_function(ax) + (y-ay)/(by-ay)*(fb_function(ax)-gb_function(ax));
    end
    counter = counter+1;
end

%for loop to get the get the values





%for loop
