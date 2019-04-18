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
%known function
f_functtion = @(x) sin*(pi*(x-ay)/(bx-ax))*cos(pi/2*(2*(y-ay)/by-ay)+1);
gb_function = @(x) (bx-x)^2*cos(pi*x)/bx;
fb_function = @(x) x*(bx-x)^2;

%matrix size
i = 2*pi;

%step size
n = i/4;
m  = i/4;

%total count
counts = i/n;

%matrix
u_maxtrix = zeros(counts,counts);

counter = 1;
for x = n:n:i
    u_maxtrix(1,counter) = fb_function(x);
    u_maxtrix(counts,counter) = gb_function(x);
    u_maxtrix(counts,counter) = gb_function(x);
    counter = counter+1;
end
%creating known values

u_maxtrix
%%
c = 1;
%get the knowns
for x = n:n:i
    for y = m:m:i
        c = 1+c;
    end
end



%for loop
