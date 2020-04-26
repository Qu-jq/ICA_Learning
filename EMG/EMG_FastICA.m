%% Author: Lu, Chia-Feng 2013.11.29
clear, close all
clc

%% initialize parameters
path(path,'.\fastica_25');  % add FastICA path

samplerate = 2000; % in Hz

ICNo=3;  % the specified number of independent components
PCNo=4;  % the preserved number of PCs 

%% load data
load('EMGdata.mat')

figure,
for i=1:size(EMGdata,2)
    subplot(size(EMGdata,2),1,i),plot(EMGdata(:,i))
end

%% perform FastICA
[icasig, A, W] = fastica(EMGdata','numOfIC',ICNo,'displayMode','off','firstEig',1,'lastEig',PCNo); % fast ICA

% [icasig, A, W] = fastica(EMGdata','numOfIC',ICNo,'displayMode','off','interactivePCA','on'); % fast ICA
% fasticag(EMGdata')

figure,
for i=1:size(icasig,1)
    subplot(size(icasig,1),1,i),plot(icasig(i,:))
end

return

%% remove ICs dominated by cardiac signal and noise
rejectICA=[1];

A2=A;
icasig2=icasig;
%%% reconstruct the signal
A2(:,rejectICA)=[];
icasig2(rejectICA,:)=[];

newEMGdata=(A2*icasig2)';

figure,
for i=1:size(newEMGdata,2)
    subplot(size(newEMGdata,2),1,i),plot(newEMGdata(:,i))
end


