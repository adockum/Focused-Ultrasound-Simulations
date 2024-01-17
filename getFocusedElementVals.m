function uamp = getFocusedElementVals(kwavenum,xyzvecs,focalPoints,focalPvals)
% Computes amplitude and phases of elements targeting the desired focus or
% multi-focus pattern, and returns a vector of complex numbers representing
% the array encoding. Code based on: Ebbini and Cain, IEEE Trans. Ultrason.
% Ferroelectr. Freq. Control, vol. 36, no. 5, pp. 540-8, Jan. 1989. 

% Inputs: 
%   kwavenum = real or complex-valued wave number, 2pi/lambda. Should be in
%   physical units consistent with xyzvecs. 
%   xyzvecs = Nx3 set of element positions that comprise the array. 
%   focalPoints = Mx3 set of desired foci. 
%   focalPvals = length M of the relative pressure values at each focus 
%   (arb units). 
% Outputs:
%   uamp = complex array of N complex amplitudes.

p = 1000;
c = 1500; 

M = length(focalPvals);
N = length(xyzvecs); 

H = 1i*ones(M,N);

for k = 1:M
    Rmnvec = sqrt(sum((xyzvecs-focalPoints(k,:)).^2,2)); 
    H(k,:) = (1j*p*c*kwavenum/(2*pi))*exp(-1j*kwavenum*Rmnvec)./Rmnvec; 
end

Hadj = conj(H');
HHa_inv = pinv(dot(H,Hadj));
uamp = Hadj.*HHa_inv'.*focalPvals;

end

