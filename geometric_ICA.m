% ICA HomeWork I
%
% Purpes 
% To Test geometric ICA
%
% Record of revision : 
%       date            programmer            description
%       ====            ==========            ===========
%    02/24/2019      Mehrdad Kashefi         original code 
% 
% Define variables : 
% num_sample         Number of samples
%...........................................................
% ****** Parameter that can be changed *****
% ------> Num Sample 
% ------> Signal type --->  S  ----> Rand-Rand  /  Rand-Sin  /  Sin-Cos
% ------> Transformation -> A  ---->  Three options
%start
% clearing the page and variables
clear;
clc;
close all;
%% Cearting Sources
num_sample = 1000; % Number of samples

% Selecting the source type
% Two randoms
S = [random('Uniform',-10,10,1,num_sample);random('Uniform',-10,10,1,num_sample)];
% Uniform and Sinusoid
%S = [10*sin(linspace(0,2*pi,num_sample));random('Uniform',-10,10,1,num_sample)];
% Two Sinusoid
%S = [10*sin(linspace(0,2*pi,num_sample));10*cos(linspace(0,2*pi,num_sample))];

% Rotation Matrix Selection
%A = [1,-0.05;0.05,1];
%A = [1,0.99;0.99,1];
A = [1,0.28;0.73,1];

X = A*S;  % Mixed matrix


figure(1)
scatter(S(1,:),S(2,:));
title('Original Signals')
xlabel('Signal 1')
ylabel('Signal 2')
figure(2)
scatter(X(1,:),X(2,:));
title('Mixed Signals')
xlabel('Signal 1')
ylabel('Signal 2')
%% Find critical samples (Edge samples)
[~,index] = sort(X(1,:)+X(2,:));
LowSample = X(:,index(1));
HighSample = X(:,index(end));
% Plot the highest and lowest points in red
figure(2)
hold on
scatter(LowSample(1),LowSample(2),'r');
scatter(HighSample(1),HighSample(2),'r');

% Construct the main line
LineGradiant = (HighSample(2)-LowSample(2))/(HighSample(1)-LowSample(1));
Line = @(x)LineGradiant*(x-LowSample(1))+LowSample(2);
LinePoint = @(x,y) abs(-y + LineGradiant*x + (-LineGradiant*LowSample(1) +LowSample(2)))/sqrt(LineGradiant^2 + 1);

% Calculate the distance and position of each point compared to the line
for Sample =1:length(X)
    Point(Sample).position = X(:,Sample);
    if Line(X(1,Sample))>X(2,Sample)
        Point(Sample).state = 'L';
    else 
        Point(Sample).state = 'U';
    end
    Point(Sample).distance = LinePoint(X(1,Sample),X(2,Sample));
end
% Sorting the distance
[~,index] = sort([Point.distance]);
Point = Point(index);
% Sroting up and down
[~,index] = sort([Point.state]);
Point = Point(index);
% Pick the farthest Upper point
MidSample = Point(end).position;

figure(2)
plot([LowSample(1) MidSample(1)], [LowSample(2) MidSample(2)],'r','LineWidth',1.2)
plot([HighSample(1) MidSample(1)], [HighSample(2) MidSample(2)],'b','LineWidth',1.2)
%% Reconstruction
AReconstruct = 1/((MidSample(2)-LowSample(2))/(MidSample(1)-LowSample(1)));
BReconstruct = (MidSample(2)-HighSample(2))/(MidSample(1)-HighSample(1));
XReconstruct = [1,AReconstruct;BReconstruct,1]\X;
%% Plotting the result
figure(1)
hold on
scatter(XReconstruct(1,:),XReconstruct(2,:));
legend('Original','Reconstructed')

figure(3)
plot(S(1,:))
hold on 
plot(XReconstruct(1,:))
title('Signal 1')
legend('Original','Reconstructed')
figure(4)
plot(S(2,:))
hold on 
plot(XReconstruct(2,:))
title('Signal 2')
legend('Original','Reconstructed')

Error = (X-XReconstruct).^2;
MSE = sum(Error(:))/numel(X);
NormMSE = MSE/num_sample;
disp(['Mean Squered Error is ', num2str(MSE)])
disp(['Normalized Mean Squered Error is ', num2str(NormMSE)])
