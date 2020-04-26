%% Author: Lu, Chia-Feng 2013.11.28
clear, close all
clc

%% initialize parameters
NumofVariable=20;
N = 10000; % data points 

%% generate rand data
X=[];
for i=1:NumofVariable
    X(i,:) = rand(1,1000)-0.5;    % uniform rand noise 
%     X(i,:) = sin(2*pi*i*[1:1000]/500);   % sine wave
end

%% mix variables and plot histogram
mixX = mean(X,1);
mixX_kurt=kurtosis(mixX)-3;  % help kurtosis for details

[count,bin]=hist(mixX,100);
figure,plot(bin, count/max(count))
title(['kurtosis (non-Gaussianity) = ' num2str(mixX_kurt)])



