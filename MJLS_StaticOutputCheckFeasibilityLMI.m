function out = MJLS_StaticOutputCheckFeasibilityLMI(A,B,C,H,Q,R,T)
%==========================================================================
% CheckFeasibilityLMI: evaluates the LMI feasibility conditions for static
% output-feedback control of Markov jump linear systems without mode
% observation.
%
% Input parameters:
% - A: system matrix; (matrix of dimension dimX x dimX x numModes)
% - B: input matrix; (matrix of dimension dimX x dimU x numModes)
% - C: output matrix; (matrix of dimension dimY x dimX x numModes)
% - H: noise matrix; (matrix of dimension dimX x dimW x numModes)
% - Q: state cost matrix; (matrix of dimension dimX x dimX x numModes)
% - R: input cost matrix; (matrix of dimension dimU x dimU x numModes)
% - T: transition matrix; (matrix of dimension numModes x numModes)
% 
% Output parameters:
% - out: out = 1: LMI feasible, i.e., MJLS is stabilizable
%        out = 0: LMI infeasible, i.e., MJLS is not stabilizable
% 
% Dependencies:
% requires Yalmip (http://users.isy.liu.se/johanl/yalmip/)
%
% Author: Maxim Dolgov
% last edited: 15 Sept. 2015, Maxim Dolgov
%==========================================================================

  addpath(genpath('yalmip'))  % change to your Yalmip path

  [dimX,dimU,numModes] = size(B);
  dimY = size(C,1);
  
  mode = zeros(numModes,1);
  mode(1) = 1;
  
  for k = 1:1e6
    mode_ = mode;
    mode = T'*mode_;
    if(sum(abs(mode-mode_)) < 1e6)
      break
    end    
  end
  if(sum(abs(mode-mode_)) > 1e6)
    error('Mode did not converge')
  end
  
  W = sdpvar(dimX+dimU,dimX+dimU,numModes);
  Y = sdpvar(dimX,dimX,numModes);
  deltaY = sdpvar(dimX,dimX,numModes);
  F = sdpvar(dimU,dimY*dimX^2);
  G = sdpvar(dimX,dimX);  
  
  D = zeros(dimX+dimU,dimX,numModes);
  E = zeros(dimX+dimU,dimU,numModes);
  
  dC = zeros(dimY*dimX^2,dimX);
  
  for i = 1:numModes
    D(1:dimX,:,i) = chol(Q(:,:,i));
    E(dimX+1:dimX+dimU,:,i) = chol(R(:,:,i));
    deltaY(:,:,i) = 0*deltaY(:,:,i);
    for j = 1:numModes
      deltaY(:,:,i) = deltaY(:,:,i) + T(j,i)*Y(:,:,j);
    end
    dC(:,:,i) = kron(eye(dimX),vec(C(:,:,i)));
  end
   
  obj = 0;
  constr = set([]);
  for i = 1:numModes
    obj = obj + trace(W(:,:,i));
    constr = constr + set([W(:,:,i), D(:,:,i)*G + E(:,:,i)*F*dC(:,:,i); (D(:,:,i)*G + E(:,:,i)*F*dC(:,:,i))' , G + G' - deltaY(:,:,i)] >= 0) ...
           + set([Y(:,:,i) - mode(i)*H(:,:,i)*H(:,:,i)', A(:,:,i)*G + B(:,:,i)*F*dC(:,:,i); (A(:,:,i)*G + B(:,:,i)*F*dC(:,:,i))' , G + G' - deltaY(:,:,i)] >= 0);
  end
  
  d = solvesdp(constr,obj,sdpsettings('verbose',0,'solver','penalb')); % you can use every other solver
  
  if(d.problem == 0)
    out = 1;
  else 
    out = 0;
  end

end % function MJLS_StaticOutputCheckFeasibilityLMI