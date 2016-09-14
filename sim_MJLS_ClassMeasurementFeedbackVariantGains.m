simTime = 50;
numRuns = 1e5;

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

A(:,:,1) = [1.2 1.2;0 1];
A(:,:,2) = [1 1.3;0 1];

B(:,:,1) = [0;1];
B(:,:,2) = [0;0.8];

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

T = [.7 .3; .6 .4];     % ti = 1

x0Mean = [5; 0];        % xi = 2

Exx0 = .1^2*eye(size(A,1));

Eww = .2^2*eye(size(A,1)); 
Evv = .1^2*eye(size(C,1));

mode0 = zeros(numModes,1);
mode0(1) = 1;
Exx = zeros(dimX,dimX,numModes);
Ex = zeros(dimX,1,numModes);
for i = 1:numModes
  Exx(:,:,i) = mode0(i)*Exx0;
  Ex(:,:,i) = mode0(i)*x0Mean;
end

controller{1} = MJLS_ClassCostaClosedLoop(A,B,H,C,J,Q,R,Eww,Evv,T,Ex,Exx,x0Mean,mode0,simTime);
controller{2} = MJLS_ClassMeasurementFeedbackConstGains(A,B,C,H,J,T,Q,R,Eww,Evv);
controller{3} = MJLS_ClassMeasurementFeedbackVariantGains(A,B,C,H,J,T,Q,R,Exx0,x0Mean,Eww,Evv,simTime,mode0,controller{2});

costs = zeros(numel(controller),1);

mode = hmmgenerate(simTime+1,T,eye(size(T,1)));
mode = [1, mode];
w = mvnrnd(zeros(size(A,1),1),Eww,simTime)';
v = mvnrnd(zeros(size(C,1),1),Evv,simTime)';
x0 = x0Mean + mvnrnd(zeros(size(A,1),1),Exx0,1)';

for i = 1:numel(controller)
  controller{i}.Reset(x0);
  x{i} = x0;
  traj{i} = zeros(size(A,1),simTime+1);
  traj{i}(:,1) = x0;
  trajEst{i} = zeros(size(A,1),simTime+1);
  trajEst{i}(:,1) = x0Mean;
  appliedU{i} = zeros(size(B,2),simTime);
end

for k = 1:simTime
  for i = 1:numel(controller)
    y{i} = C(:,:,mode(k))*x{i} + J(:,:,mode(k))*v(:,k);

    [u{i},xe{i}] = controller{i}.DataInDataOut(y{i},mode(k));

    x{i} = A(:,:,mode(k))*x{i} + B(:,:,mode(k))*u{i} + H(:,:,mode(k))*w(:,k);

    traj{i}(:,k+1) = x{i};
    trajEst{i}(:,k+1) = xe{i};
    appliedU{i}(:,k) = u{i};
  end
end

for i = 1:numel(controller)
  costs(i) = ComputeCosts(Q,R,mode,traj{i},appliedU{i});
end
costs


colors = {'r','g','b--'};
h = zeros(1,numel(controller));
figure
subplot(4,1,1)
plot(zeros(1,simTime+1),'k')
hold on
for i = 1:numel(controller)
  h(i) = plot(traj{i}(1,:),colors{i},'LineWidth',1.2);
end
axis([0 simTime -5 12])
ylabel('$$x_1$$','Interpreter','Latex','FontSize',12)
l = legend(h,'Optimal [7]','Time-invariant [19]','Proposed');
set(l,'Interpreter','Latex','FontSize',12)

subplot(4,1,2)
plot(zeros(1,simTime+1),'k')
hold on
for i = 1:numel(controller)
  h(i) = plot(traj{i}(2,:),colors{i},'LineWidth',1.2);
end
axis([0 simTime -4 2])
ylabel('$$x_2$$','Interpreter','Latex','FontSize',12)

subplot(4,1,3)
plot(zeros(1,simTime+1),'k')
hold on
for i = 1:numel(controller)
  h(i) = plot(appliedU{i}(1,:),colors{i},'LineWidth',1.2);
end
axis([0 simTime -4 2])
ylabel('$$u$$','Interpreter','Latex','FontSize',12)

subplot(4,1,4)
plot(mode,'bx','LineWidth',1.2)
axis([0 simTime .5 2.5])
ylabel('$$\theta_k$$','Interpreter','Latex','FontSize',12)
xlabel('Time','Interpreter','Latex','FontSize',12)