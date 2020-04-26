function sigfilter=filter_2sIIR(sig,f,fs,n,type)
% A two-stage IIR filter processing using Butterworth filter design
%
% The  input arguments are identical to those used in BUTTER
%  sig: is the signal to be filtered arranged, each row is a measurement
%    f: a vector of cutoff frequencies in Hz. f is a number for lowpass and
%       highpass filters, and is a two-element vector for bandpass or bandstop
%       filter.
%   fs: the sampling frequency.
%    n: the order of IIR filter.
% type: 'low' for lowpass filter; 'high' for highpass filter; 
%       'bandpass' for bandpass filter; 'stop' for bandstop filter.
%
%  Details of input arguments and examples, please help BUTTER.
%
%  [Author]: Lu, Chia-Feng 2013.11.08

if nargin<4
    error('Not enough input arguments! pleas check') 
elseif nargin<5
    type='low';
end

if f>=fs/2
   error('The sampling frequency is not adequent for the given cutoff frequency, please input a lower f.') 
end

if strcmpi(type,'bandpass')
    %%% highpass
    [b,a]=butter(n,f(1)/(fs/2),'high');
    sigfilter=[];
    for i=1:size(sig,1)
        sigfilter(i,:) = filtfilt(b,a,sig(i,:));
    end
    
    %%% lowpass
    [b,a]=butter(n,f(end)/(fs/2),'low');
    for i=1:size(sigfilter,1)
        sigfilter(i,:) = filtfilt(b,a,sigfilter(i,:));
    end
else
    [b,a]=butter(n,f/(fs/2),type);

    sigfilter=[];
    for i=1:size(sig,1)
        sigfilter(i,:) = filtfilt(b,a,sig(i,:));
    end
end
