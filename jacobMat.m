function[hVec, Hvec] = jacobMat(w)

% f=[10 20 30];
% Fs=200;
% w = 2*pi*f/Fs;


p = length(w);
roots = exp(1i.*w);
A = diag(roots);
Hvec = zeros(p);

coeffs = -1*poly(roots);
hVec = coeffs;
hVec(1) = [];
for i=1:p
    tempRoots = roots;
    tempRoots(i) = [];
    coeff = poly(tempRoots);
    Hvec(i,:) = coeff;
end

Hvec = 1i*A*Hvec;
end

