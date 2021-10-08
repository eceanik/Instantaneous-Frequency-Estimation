function[signal,t,fVec,snrDb] = inputNonStationary(Fs,sigFreqStart,sigFreqEnd,L,var_noise,time,type)
%% Create non-stationary input signal
% Author: Anik Kumar Samanta (ece.anik@gmail.com)
% Inputs: Signal frequencies, Signal Amplitudes, Length of signal, noise
% variance, and sampling frequency

% Output: Signal and time information
% tic();
Ts = 1/Fs; % Sample time
t = (0:L-1)*Ts; % Time vector

% noise = var_noise*randn(size(t));
noise = normrnd(0,var_noise,1,L);

signal = [];

switch type
    case 'freqLin'
        signal = chirp(t,sigFreqStart,time,sigFreqEnd,'linear');
        beta = (sigFreqEnd - sigFreqStart)/time;
        fVec = sigFreqStart + beta*t;
        
    case 'freqQuad'
        signal = chirp(t,sigFreqStart,time,sigFreqEnd,'quadratic');
        beta = (sigFreqEnd - sigFreqStart)/time^2;
        fVec = sigFreqStart + beta*t.^2;

    case 'amp'
        n = randi(10);              % random no. of signl segments
        l = randi(ceil(L/n),n,1);
        l = [l;L-sum(l)];           % Length of each segment
        amp = rand(1,length(l));
        phase = unifrnd(-pi,pi,1);%phase = 0;
        for i = 1:length(l)
            tNew = (0:l(i)-1)*Ts;
            temp = amp(i)*sin(2*pi*sigFreqStart/Fs*tNew + phase);
            signal = [signal,temp];
        end
        
    case 'phase'
        n = randi(10);              % random no. of signl segments
        l = randi(ceil(L/n),n,1);
        l = [l;L-sum(l)];           % Length of each segment
        phase = unifrnd(-pi,pi,[1,length(l)]);
        for i = 1:length(l)
            tNew = (0:l(i)-1)*Ts;
            temp = 1*sin(2*pi*sigFreqStart/Fs*tNew + phase(i));
            signal = [signal,temp];
        end
end
        
        
        snrDb = snr(signal,noise);
        signal = signal + noise;
        signal = signal';
        
end
