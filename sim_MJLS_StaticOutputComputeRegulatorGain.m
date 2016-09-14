
dimX = 2;
dimU = 1;
dimY = 3;
dimW = 2;
numModes = 2;

A = zeros(dimX,dimX,numModes);
B = zeros(dimX,dimU,numModes);
C = zeros(dimY,dimX,numModes);
H = zeros(dimX,dimW,numModes);
J = zeros(dimY,dimV,numModes);
Q = zeros(dimX,dimX,numModes);
R = zeros(dimU,dimU,numModes);

A(:,:,1) = [1.2 1.2;0 1];
A(:,:,2) = [1 .8;0 1];

B(:,:,1) = [0;1];
B(:,:,2) = [0;.2];

C(:,:,1) = [1 .2; 1 .6; .2 1];
C(:,:,2) = [1 1; 1 .6; .2 1];

H(:,:,1) = [.2 0;0 .1 ];
H(:,:,2) = H(:,:,1);

Q(:,:,1) = eye(size(A,1));
Q(:,:,2) = Q(:,:,1);

R(:,:,1) = eye(size(B,2));
R(:,:,2) = R(:,:,1);

Eww = 1^2*eye(size(A,1));

% T = [.7 .3; .6 .4];
T = [.9 .1; .1 .9];
% T  = [.2 .8; .3 .7];

tic
L = MJLS_StaticOutputComputeRegulatorGain(A,B,C,H,Q,R,T);
toc

for i = 1:numModes
  At(:,:,i) = A(:,:,i)+B(:,:,i)*L*C(:,:,i);
end

fprintf('Spectral radius iterative algorithm: %f\n',ComputeSpectralRadius(At,T));