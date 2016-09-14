function [L,Exx] = MJLS_StaticOutputComputeRegulatorGain(A,B,C,H,Q,R,T)
%==========================================================================
% MJLS_StaticOutputComputeRegulatorGain: computes the regulator gain for static output-
% feedback control of Markov jump linear systems without mode observation.
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
% - L: regulator gain; (matrix of dimension dimU x dimY)
% 
% Author: Maxim Dolgov
% last edited: 09 Sept. 2015, Maxim Dolgov
%==========================================================================

  [dimX,dimU,numModes] = size(B);
  dimY = size(C,1);
  
  mode = zeros(numModes,1);
  mode(1) = 1;
  
  for k = 1:1e6
    mode_ = mode;
    mode = T'*mode_;
    if(sum(abs(mode-mode_)) < 1e-6)
      break
    end    
  end
  if(sum(abs(mode-mode_)) > 1e-6)
    error('Mode did not converge')
  end
  
  % initialization
  Exx = zeros(dimX,dimX,numModes);
  La = zeros(dimX,dimX,numModes);
  for i = 1:numModes
    Exx(:,:,i) = 1e-4*eye(dimX);
    La(:,:,i) = 1e-4*eye(dimX);
  end  
  oldL = Inf*ones(dimU,dimY);
  % END: initialize
  
  % iteration until convergence
    for k = 1:1e8
      P = OpEpsilon(La,T);
      
      % compute L
      XRBPB = 0;
      BPAX = 0;
      for i = 1:numModes
        XRBPB = XRBPB + kron(C(:,:,i)*Exx(:,:,i)*C(:,:,i)',R(:,:,i)+B(:,:,i)'*P(:,:,i)*B(:,:,i));
        BPAX = BPAX + B(:,:,i)'*P(:,:,i)*A(:,:,i)*Exx(:,:,i)*C(:,:,i)';
      end

      L = reshape(-XRBPB\BPAX(:),[dimU,dimY]);
      % END: compute L

      Exx_ = RecursionExx(Exx,L,A,B,C,H,T,mode);
      La_ = RecursionLambda(La,L,A,B,C,Q,R,T,mode);
      
      % check convergence      
        if(k>1)
          if(sum(sum(abs(oldL-L))) < 1e-8)
            break
          end
        end        
      % END: check convergence
      
      Exx = Exx_;
      La = La_;
      oldL = L;
    end   
  % END: % iteration until convergence
end % function MJLS_StaticOutputComputeRegulatorGain

%% RecursionExx
function outExx = RecursionExx(Exx,L,A,B,C,H,T,mode)
  outExx = zeros(size(Exx));
  numModes = numel(mode);
  
  for j = 1:numModes
    for i = 1:numModes
      outExx(:,:,j) = outExx(:,:,j) + T(i,j)*((A(:,:,i)+ B(:,:,i)*L*C(:,:,i))*Exx(:,:,i)*(A(:,:,i) + B(:,:,i)*L*C(:,:,i))' + mode(i)*H(:,:,i)*H(:,:,i)');
    end
    outExx(:,:,j) = (outExx(:,:,j) + outExx(:,:,j)')/2;
  end
end % function RecursionExx

%% RecursionLambda
function outP = RecursionLambda(P,L,A,B,C,Q,R,T,mode)
  outP = zeros(size(P));
  numModes = numel(mode);
  
  eP = OpEpsilon(P,T);
  
  for j = 1:numModes
    outP(:,:,j) = (A(:,:,j) + B(:,:,j)*L*C(:,:,j))'*eP(:,:,j)*(A(:,:,j)+ B(:,:,j)*L*C(:,:,j)) + Q(:,:,j) + C(:,:,j)'*L'*R(:,:,j)*L*C(:,:,j);
    outP(:,:,j) = (outP(:,:,j) + outP(:,:,j)')/2;
  end
end % function RecursionLambda

%% OpEpsilon
function out = OpEpsilon(M,T)
  [numR,numC,nModes] = size(M);

  out = zeros(numR,numC,size(T,1));
  for i = 1:nModes
    for j = 1:nModes
      out(:,:,i) = out(:,:,i) + T(i,j)*M(:,:,j);
    end
  end
end % function OpEpsilon