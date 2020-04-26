%% Author: Lu, Chia-Feng 2013.11.28
clear, close all
clc

%% initialize parameters
samplerate=500; % in Hz
N=1000; % data length

freq1=5; % in Hz
freq2=7; % in Hz
taxis=[1:N]/samplerate;

PCnum=2;  % the number of PC used to reconstruct signals

%% generate test signals 
C1 = 0.75*sin(2*pi*freq1*taxis);        % 1st component: a sine wave
C2 = sawtooth(2*pi*freq2*taxis,0.5);    % 2nd component: a triangular wave

% Combine data in different proportions
X(1,:) = 0.5*C1 + 0.5*C2 + 0.1*rand(1,N);
X(2,:) = 0.7*C1 + 0.2*C2 + 0.1*rand(1,N);
X(3,:) = 0.2*C1 + 0.7*C2 + 0.1*rand(1,N);
X(4,:) = -0.3*C1 - 0.6*C2 + 0.3*rand(1,N);
X(5,:) = 0.6*rand(1,N);    % Noise only

% Center data by subtracting mean
X = X - mean(X,2)*ones(1,N);

figure, 
for i=1:size(X,1)
    subplot(size(X,1),1,i)
    plot(taxis,X(i,:)),xlim([taxis(1) taxis(end)])
end

%% Principal Components Analysis using Singular Value Decomposition (SVD)
[U,S,V]=svd(X,0);   % perform SVD on Signals

eigenVal=diag(S).^2;  % calculate eigenvalues
for i=1:size(X,1)
    PC(i,:)=V(:,i)*sqrt(eigenVal(i));   % scale principal components
end

eigenVal=eigenVal/(N-1); % now eigenvalues equal to variances
for i=1:size(X,1)
    eigen_perc(i)=sum(eigenVal(1:i))/sum(eigenVal)*100;  % calculate accumulated percentage of eigenvalues
end

figure, 
for i=1:size(PC,1)
    subplot(size(PC,1),1,i)
    plot(taxis,PC(i,:)),xlim([taxis(1) taxis(end)])
end

figure,plot(eigen_perc,'-o')
xlabel('dimension'),ylabel('percentage of accumulate eigenvalues')

%% Check the covariance of principal components (PC)
cov(PC')  % make sure if the PCs are uncorrelated !
% the off-diagonal terms should be 0

%% Reconstruct the Signal only use the first PCnum PCs
newX = U(:,1:PCnum) * S(1:PCnum,1:PCnum) *V(:,1:PCnum)';

figure, 
for i=1:size(newX,1)
    subplot(size(newX,1),1,i)
    plot(taxis,X(i,:)),hold on
    plot(taxis,newX(i,:),'r'),xlim([taxis(1) taxis(end)])
end

