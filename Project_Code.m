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
%reset varible to ensure values
x = 0;
y = 0;
%using temp to test logic
u_maxtrix_temp = u_maxtrix;
%need to insert ghost node
u_maxtrix_ghost_node = 0;

[r c] = size(u_maxtrix_temp);
for k = 1:1:r
    for j = 1:1:c
        %top boundary
        if k == 1
            
        %bottom boundary
        elseif k == r
            
        %left
        elseif j == 1
            
        %right
        elseif j == c
            
        else
             u_maxtrix_temp(k,j) = u_maxtrix_temp(k,j-1) + u_maxtrix_temp(k,j+1) + u_maxtrix_temp(k-1,j) + u_maxtrix_temp(k+1,j)
        end
    end
end



%for loop
