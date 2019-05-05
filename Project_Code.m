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
f_function = @(x,y) sin(pi*(x-ay)/(bx-ax))*cos(pi/2*(2*(y-ay)/by-ay)+1);
gb_function = @(x) (bx-x)^2*cos(pi*x)/bx;
fb_function = @(x) x*(bx-x)^2;

%matrix size
i = 2*pi;
%start node
node_start = 5;
%ending node
node_end = 50;
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
    %pre-locate matrix
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
        %pre-allocate
        [r c] = size(u_maxtrix_temp);
        %calculation of the value using the formula of Gauss
        for j = 2:1:(r-1)
            for k = 2:1:(c-1)
                u_maxtrix_temp(k,j) = 1/4 *(u_maxtrix_temp(k-1,j)+u_maxtrix(k+1,j)+u_maxtrix_temp(k,j-1)+u_maxtrix(k+1,j+1))-1/4*f_function(k,j);
            end
        end
        %calculate the error
        errors = sqrt(sum(sum((u_maxtrix_temp-u_maxtrix).^2)));
        %saving the matrix
        u_maxtrix = u_maxtrix_temp;
        %increase the interative count
        interative = interative +1;
    end
    %record the matrix to the overrall matrix holder
    matrix_holder_after{count_temp} = u_maxtrix;
    %increase the count of the overrall matrix holder counter
    count_temp = count_temp + 1;
    fprintf('We went thorugh this many round for Gauss-Seidel Method %f at this many nodes: %f\n', interative,nodes_n)
end
%mesh the last plot
mesh(u_maxtrix)
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
f_functtion = @(x,y) sin(pi*(x-ay)/(bx-ax))*cos(pi/2*(2*(y-ay)/by-ay)+1);
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
matrix_holder = cell(1,node_end-4);
matrix_holder_after = cell(1,node_end-4);
length_omega = length(1:0.1:2);

counter_omega = 1;
omega_tracker = zeros(node_start-node_end+1,length_omega);
omega_matrix_holder = cell(1,node_end-4);
omega_and_matrix_holder = cell(1,node_end-4);
omega_and_matrix_counter = 0;

%initlize the first knows
for nodes_n = node_start:node_increment:node_end
    counter_omega = 1;
    for omega_w = 1:0.1:2
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
                    u_maxtrix_temp(k,j) = (1-omega_w)*u_maxtrix(k,j)+omega_w/4 *(u_maxtrix_temp(k-1,j)+u_maxtrix(k+1,j)+u_maxtrix_temp(k,j-1)+u_maxtrix(k+1,j+1))-omega_w/4*f_functtion(k,j);
                end
            end

            errors = sqrt(sum(sum((u_maxtrix_temp-u_maxtrix).^2)));
            u_maxtrix = u_maxtrix_temp;
            interative = interative +1;
           
        end
        %omega tracker
        omega_tracker(counter_omega,nodes_n-4) = interative;
  
        %matrix tracket
        omega_matrix_holder{counter_omega,nodes_n-4} = u_maxtrix;
        counter_omega = counter_omega + 1;
       
        fprintf('We went thorugh this many round for Successive Over-Relaxation Method %f with omega being %f\n', interative,omega_w)
    end
    matrix_holder_after{count_temp} = {omega_matrix_holder(:,count_temp) omega_tracker(:,count_temp)};
    count_temp = count_temp + 1;

end
%mesh the last plot
%get the min omeage value for the iteration
iteration_min = matrix_holder_after{1,length(matrix_holder_after)}{2};
%find the location of the min
[x,y]=find(iteration_min==min(iteration_min));
%get the matrix associated with that omega
matrix_overrall_omega = matrix_holder_after{1,length(matrix_holder_after)}{1};
%assign the associated to u matrix
u_maxtrix = matrix_overrall_omega{x(1)};
mesh(u_maxtrix);





