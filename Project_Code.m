%project MECE 
%Toan Nguyen

%% Method 1
%Gauss-Seidel Method
%housekeepin

%housekeepin
clear all
close all
clc

%assumption
error_min = 10^-4;
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
%start node
node_start = 5;
%ending node
node_end = 6;
%increment
node_increment = 1;

count_temp = 1;
matrix_holder = cell(1,3);
matrix_holder_after = cell(1,3);

%initlize the first knows
for nodes_n = node_start:node_increment:node_end
    %step size
    n = i/nodes_n;
    m  = i/nodes_n;
    %total count
    sizes = round(i/n);
    %matrix
    u_maxtrix = zeros(sizes,sizes);
    counter = 1;
    for t = 1:1:sizes
        %since step size is all the same
        x = t*n;
        y = x*n;
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
    matrix_holder{count_temp} = u_maxtrix;
    %for loop to get the get the values
    %reset varible to ensure values
    x = 0;
    y = 0;
    interative = 0;
    errors = 10;
    %using temp to test logic
    u_maxtrix1 = u_maxtrix;
    u_maxtrix_temp = u_maxtrix;
    while error_min <= errors
        %need to insert ghost node
        u_maxtrix_ghost_node = 0;

        [r c] = size(u_maxtrix_temp);
        for k = 2:1:(r-1)
            for j = 2:1:(c-1)
                u_maxtrix_temp(k,j) = 1/4 *(u_maxtrix_temp(k-1,j)+u_maxtrix(k+1,j)+u_maxtrix_temp(k,j-1)+u_maxtrix(k+1,j+1));
            end
        end

        errors = sqrt(sum(sum((u_maxtrix_temp-u_maxtrix).^2)));
        u_maxtrix = u_maxtrix_temp;
        interative = interative +1;
        
        
    end
    matrix_holder_after{count_temp} = u_maxtrix;
    count_temp = count_temp + 1;
    fprintf('We went thorugh this many round for Gauss-Seidel Method %f\n', interative)
end

%% Method 2

%Successive Over-Relaxation

%housekeepin
clear all
close all
clc

%assumption
error_min = 10^-4;

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
%start node
node_start = 5;
%ending node
node_end = 5;
%increment
node_increment = 1;

count_temp = 1;
matrix_holder = cell(1,3);
matrix_holder_after = cell(1,3);
length_omega = length(1:0.01:2);
counter_omega = 1;
omega_tracker = zeros(1,length_omega);
%initlize the first knows
for nodes_n = node_start:node_increment:node_end
    for omega_w = 1:0.01:2
        %step size
        n = i/nodes_n;
        m  = i/nodes_n;
        %total count
        sizes = round(i/n);
        %matrix
        u_maxtrix = zeros(sizes,sizes);
        counter = 1;
        for t = 1:1:sizes
            %since step size is all the same
            x = t*n;
            y = x*n;
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
        matrix_holder{count_temp} = u_maxtrix;
        %for loop to get the get the values
        %reset varible to ensure values
        x = 0;
        y = 0;
        interative = 0;
        errors = 10;
        %using temp to test logic
        u_maxtrix1 = u_maxtrix;
        u_maxtrix_temp = u_maxtrix;
        while error_min <= errors
            %need to insert ghost node
            u_maxtrix_ghost_node = 0;

            [r c] = size(u_maxtrix_temp);
            for k = 2:1:(r-1)
                for j = 2:1:(c-1)
                    u_maxtrix_temp(k,j) = (1-omega_w)*u_maxtrix(k,j)+(omega_w/4)*(u_maxtrix_temp(k-1,j)+u_maxtrix(k+1,j)+u_maxtrix_temp(k,j-1)+u_maxtrix(k+1,j+1));
%                     u_maxtrix_temp(k,j) = 1/4 *(u_maxtrix_temp(k-1,j)+u_maxtrix(k+1,j)+u_maxtrix_temp(k,j-1)+u_maxtrix(k+1,j+1));
                end
            end

            errors = sqrt(sum(sum((u_maxtrix_temp-u_maxtrix).^2)));
            u_maxtrix = u_maxtrix_temp;
            interative = interative +1;
           
        end
        omega_tracker(counter_omega) = interative;
        counter_omega = counter_omega + 1;
        fprintf('We went thorugh this many round for Successive Over-Relaxation Method %f\n', interative)
    end
    matrix_holder_after{count_temp} = u_maxtrix;
    count_temp = count_temp + 1;
end






