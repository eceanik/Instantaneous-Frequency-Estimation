%%% 
% Instantaneous Frequency Estimator Main File
% Input: Time data (inputNonStationary.m)
% Author: Anik Kumar Samanta


clear all

L = 2000;
Fs = 200;
noiseVar = 0.001/sqrt(2);
epsilon = noiseVar;

s = [10*2*pi/Fs,15*2*pi/Fs,20*2*pi/Fs]';

pPrimary = length(s);
nSin = pPrimary;

sigmaState = 1.2;

sigmaMeas = 0.1;

pVar = 1;


[x,t,freqTh1] = inputNonStationary(Fs,10,20,L,noiseVar,10,'freqLin');
[y,~,freqTh2] = inputNonStationary(Fs,20,30,L,noiseVar,10,'freqLin'); y = [x+y];
[z,~,freqTh3] = inputNonStationary(Fs,15,25,L,noiseVar,10,'freqLin'); y = [y+z];

y = (hilbert(y));
F = eye(pPrimary);
G = eye(pPrimary);

Q = sigmaState*eye(pPrimary);
R = sigmaMeas;

P = pVar*eye(pPrimary);


for i = pPrimary+1:L
    
    dataSt = i - pPrimary;
    dataEn = i - 1;
    
    X = y(dataEn:-1:dataSt);
    
    [hVec, Hvec] = jacobMat(s);
    h = hVec*X;
    H = Hvec*X; 
    
    H = transpose(H);
    s = F*s;
    
    P = F*P*F' + Q;
    
    K = P*H'/(H*P*H' + R);
    
    s = s + K*(y(i) - h);

    s = wrapTo2Pi(abs(s));
    
    P = eye(pPrimary) - K*H*P;
    
    omega(:,i) = s;
   
end
f = (omega*Fs/(2*pi));
 
plot(t,freqTh1,'--',t,freqTh2,'--',t,freqTh3,'--',t,f(1,:),t,f(2,:),t,f(3,:),'linewidth',1.5);

xlabel('Time (s)');
ylabel('Frequency (Hz)');





