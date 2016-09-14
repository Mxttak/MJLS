% simmain: sumulation script for a reference implementation of the
% control law presented in M. Dolgov, U.D. Hanebeck, "Dynamic Compensation 
% of Markov Jump Linear Systems without Mode Observation", 2015

simTime = 200;

dimX = 2;
dimU = 1;
dimY = 1;
dimW = 2;
dimV = 1;
numModes = 2;

A = zeros(dimX,dimX,numModes);
B = zeros(dimX,dimU,numModes);
C = zeros(dimY,dimX,numModes);
H = zeros(dimX,dimW,numModes);
J = zeros(dimY,dimV,numModes);

Q = zeros(dimX,dimX,numModes);
R = zeros(dimU,dimU,numModes);

A(:,:,1) = [1 1;0 1];
A(:,:,2) = [1.2 1.2;0 1];

B(:,:,1) = [0;1];
B(:,:,2) = [.01;.8];

C(:,:,1) = [1 0];
C(:,:,2) = [.8 0];

H(:,:,1) = [1 0;0 1];
H(:,:,2) = H(:,:,1);

J(:,:,1) = 1;
J(:,:,2) = J(:,:,1);

Q(:,:,1) = eye(size(A,1));
Q(:,:,2) = Q(:,:,1);

R(:,:,1) = eye(size(B,2));
R(:,:,2) = R(:,:,1);

T = [.7 .3; .2 .8];
% T = [.9 .1; .1 .9];

% x0 = [0; 0];
x0 = [2; 0];

Exx0 = .2^2*eye(dimX);
% Exx0 = .01^2*eye(dimX);

Eww = .1^2*eye(size(A,1)); 
% Eww = .01^2*eye(size(A,1)); 

Evv = .05^2*eye(size(C,1));
% Evv = .1^2*eye(size(C,1));


mode0 = ones(numModes,1)/numModes;

controller = MJLS_ClassMeasurementFeedbackConstGains(A,B,C,H,J,T,Q,R,Eww,Evv);

mode = hmmgenerate(simTime+1,T,eye(size(T,1)));
w = mvnrnd(zeros(size(A,1),1),Eww,simTime)';
v = mvnrnd(zeros(size(C,1),1),Evv,simTime)';

controller.Reset(x0);
x = mvnrnd(x0,Exx0)';
traj = zeros(size(A,1),simTime+1);
traj(:,1) = x;
trajEst = zeros(size(A,1),simTime+1);
trajEst(:,1) = x0;
appliedU = zeros(size(B,2),simTime);

for k = 1:simTime
    y = C(:,:,mode(k))*x + J(:,:,mode(k))*v(:,k);

    [u,xe] = controller.DataInDataOut(y);

    x = A(:,:,mode(k))*x + B(:,:,mode(k))*u + H(:,:,mode(k))*w(:,k);

    traj(:,k+1) = x;
    trajEst(:,k+1) = xe;
    appliedU(:,k) = u;
end
