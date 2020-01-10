% ICA with EASI method
% Title: An Equivariant Source Separation
% To calculate ICA for two sources based on EASI algorithm
%
%      date            programmer             version
%     =======         ===========            =========
%   06/15/2019       Mehrdad Kashefi       original code
%  ..........................................................
% ...........................................................
%% Clearing 
clc;
clear;
close all;
%% Control Parameters
% Perform whittening on data 
sw_white = 0;  % 0  --> No Whittening
               % 1  --> Perform Whittening
normalize = 0;   
serial = 1;
max_iter = 500;  % Maximum number of iterations
rate = 0.0001;
% Non_linear function 
g = @(x) tanh(0.01*x);
% Mixing sources
A = [1,0.3; 0.4,1 ];  % Non Singular mixing matrix
%A = [1,0.7;0.99,1];    % Close to singulat mixing matrix
%% Create two sample sources
T = 10*(1/50);
Fs = 1000;
dt = 1/Fs;
t = 0:dt:T-dt;
s1 = sawtooth(2*pi*50*t);   % A sawtooth source
s2 = sin(2*pi*22*t);        % A sinosiod source
% Plot the original sources
figure(1)
subplot(3,1,1)
plot(t,s1)
grid on
title("Source I")
subplot(3,1,2)
plot(t,s2)
title("Source II")
subplot(3,1,3)
scatter(s1,s2)
grid on
title("Sources Scatter Plot")

x = A*[s1;s2];

figure(2)
subplot(3,1,1)
plot(t,x(1,:))
grid on
title('Mixed x1')
subplot(3,1,2)
plot(t,x(2,:))
grid on
title('Mixed x2')
subplot(3,1,3)
scatter(x(1,:),x(2,:))
grid on
title("Mixed Sources Scatter Plot")
%% ICA Preprocessing
% Zero mean the signals
x_mean = mean(x,2);
x = x - x_mean;
% Whitening the signal
if sw_white==1
    cov_x = x*x';
    [eig_vec,eig_val] = eig(cov_x);
    x_white = (eig_val^-0.5)*eig_vec'*x;
    disp('Data Covariance before whittening');
    disp(x*x');
    disp('Data Covariance after whittening');
    disp(x_white*x_white');
    
    % Plot whitened data
    figure(3)
    subplot(3,1,1)
    plot(t,x_white(1,:))
    grid on
    title('Whittened x1')
    subplot(3,1,2)
    plot(t,x_white(2,:))
    grid on
    title('Whittened x2')
    subplot(3,1,3)
    scatter(x_white(1,:),x_white(2,:))
    grid on
    title("Whittened Scatter Plot")
end

B = eye(2,2);

for iter = 1:max_iter
    y = B*x;
    if normalize
        H = (y*y' - eye(2,2))./(1+rate*(y*y')) + (g(y)*y' - y*g(y)')./(1+rate*det(y'*g(y)));
    else
        H = (y*y' - eye(2,2)) + (g(y)*y' - y*g(y)');
    end

    B = B - rate*H*B';

     figure(4)
    subplot(3,1,1)
    plot(t,y(1,:))
    grid on
    title('Reconstruction y1')
    subplot(3,1,2)
    plot(t,y(2,:))
    grid on
    title('Reconstruction y2')
    subplot(3,1,3)
    scatter(y(1,:),y(2,:))
    grid on
    title("Reconstruction Scatter Plot")

end


