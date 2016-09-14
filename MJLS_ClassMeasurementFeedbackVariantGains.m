classdef MJLS_ClassMeasurementFeedbackVariantGains < handle
% MJLS_ClassMeasurementFeedbackVariantGains: a reference implementation of
% the control law presented in M. Dolgov, G. Kurz, U.D. Hanebeck,
% "Finite-horizon Dynamic Compensation of Markov Jump Linear Systems
% without Mode Observation", 2016.
%
% This implementation is provided as is. For academic use only.
%
% AUTHOR: Maxim Dolgov, 09 Mar. 2016
%         maxim.dolgov@kit.edu

properties (SetAccess = protected)
  A = [];
  B = [];
  H = [];
  C = [];
  J = [];
  Q = [];
  R = [];
  Exx0 = [];
  Ex0 = [];
  Eww = [];
  Evv = [];
  T = [];
  simTime = -1;
  mode0 = [];
  
  dimX = -1;
  dimU = -1;
  dimW = -1;
  dimY = -1;
  dimV = -1;
  numModes = -1;
  
  mode = [];
  
  L = [];
  F = [];
  K = [];
  
  x = [];
  t = -1;
end % properties

properties(Constant)
  costTol = 1e-6;
  maxIter = 1e8;
  maxDiff = 1e-8;
end % properties constant


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

methods
  
  %% MJLS_ClassMeasurementFeedbackVariantGains
  function s = MJLS_ClassMeasurementFeedbackVariantGains(A,B,C,H,J,T,Q,R,Exx0,Ex0,Eww,Evv,simTime,mode0)
    
    s.A = A;
    s.B = B;
    s.H = H;
    s.C = C;
    s.J = J;
    s.Q = Q;
    s.R = R;
    s.Exx0 = Exx0;
    s.Ex0 = Ex0;
    s.Eww = Eww;
    s.Evv = Evv;
    s.T = T;
    s.simTime = simTime;
    s.mode0 = mode0;
    
    [s.dimX,s.dimU,s.numModes] = size(B);
    s.dimW = size(Eww,1);
    [s.dimY,s.dimV,~] = size(J);
    
    s.mode = cell(1,simTime+1);
    s.mode{1} = mode0;
    for k = 1:s.simTime
      s.mode{k+1} = s.T'*s.mode{k};
      s.mode{k} = reshape(s.mode{k},[1 1 s.numModes]);
      if(any(isnan(s.mode{k+1})) || any(isinf(s.mode{k+1})) || abs(sum(s.mode{k+1})-1) > s.maxDiff)
        error('Computation of the mode fucked up')
      end
    end
    s.mode{s.simTime+1} = reshape(s.mode{s.simTime+1},[1 1 s.numModes]);  % i-th mode at time step k is now accessible by s.mode{k}(:,:,i)
    
    [s.F,s.K,s.L]  = s.ComputeController();
    
    s.t = 1;
  end % function MJLS_ClassMeasurementFeedbackVariantGains
  
  %% DataInDataOut
  function [u,xe] = DataInDataOut(s,y)
    u = s.L{s.t}*s.x;
    s.x = s.F{s.t}*s.x + s.K{s.t}*y;
    xe = s.x;
    s.t = s.t +1;
  end % function DataInDataOut
  
  %% Reset
  function Reset(s,x0)
    s.x = x0;
    s.t = 1;
  end % function Reset
  
end % methods public

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

methods (Access = protected)
  %% ComputeController
  function [F,K,L] = ComputeController(s)
    % initialization      
      F = cell(1,s.simTime);
      K = cell(1,s.simTime);
      L = cell(1,s.simTime);
      for k = 1:s.simTime
        F{k} = rand(s.dimX)-.5;
        K{k} = rand(s.dimX,s.dimY)-.5;
        L{k} = rand(s.dimU,s.dimX)-.5;
      end
      
    % END: initialization

    
    % iterate until convergence   
    oldCosts = Inf;
    for iter = 1:s.maxIter
      % predict Xu and Xl with current F,K,L
        [EXu,EXl] = s.PredictXuXl(F,K,L);
      % END: predict Xu and Xl with current F,K,L
      
      % initialize Lu and Ll
        Lu = s.Q;
        Ll = zeros(s.dimX,s.dimX,s.numModes);
        
        oCosts = zeros(1,1,s.numModes);
        
      % END: initialize Lu and Ll
      
      for k = s.simTime:-1:1
        Pu = s.OpEpsilon(Lu,s.T);
        Pl = s.OpEpsilon(Ll,s.T);
        Xu = EXu{k};
        Xl = EXl{k};

        % update F, K, L
          XRBPB = 0;        BPuAX = 0;        BPlAX = 0;
          XCBP = 0;         XBP = 0;          XP = 0;
          XCP = 0;          PAX = 0;          JCXC = 0;
          PAXuC = 0;        PAXlC = 0;
          for i = 1:s.numModes
            XRBPB = XRBPB + kron(Xl(:,:,i),s.R(:,:,i) + s.B(:,:,i)'*Pu(:,:,i)*s.B(:,:,i) + s.B(:,:,i)'*Pl(:,:,i)*s.B(:,:,i));
            BPuAX = BPuAX + s.B(:,:,i)'*Pu(:,:,i)*s.A(:,:,i)*Xl(:,:,i);
            BPlAX = BPlAX + s.B(:,:,i)'*Pl(:,:,i)*s.A(:,:,i)*Xl(:,:,i);
            XCBP = XCBP + kron(Xl(:,:,i)*s.C(:,:,i)',s.B(:,:,i)'*Pl(:,:,i));
            XBP = XBP + kron(Xl(:,:,i),s.B(:,:,i)'*Pl(:,:,i));
            XP = XP + kron(Xl(:,:,i),Pl(:,:,i));
            XCP = XCP + kron(Xl(:,:,i)*s.C(:,:,i)',Pl(:,:,i));
            PAX = PAX + Pl(:,:,i)*s.A(:,:,i)*Xl(:,:,i);
            JCXC = JCXC + kron(s.mode{k}(:,:,i)*s.J(:,:,i)*s.J(:,:,i)' + s.C(:,:,i)*Xu(:,:,i)*s.C(:,:,i)' + s.C(:,:,i)*Xl(:,:,i)*s.C(:,:,i)',Pl(:,:,i));
            PAXuC = PAXuC + Pl(:,:,i)*s.A(:,:,i)*Xu(:,:,i)*s.C(:,:,i)';
            PAXlC = PAXlC + Pl(:,:,i)*s.A(:,:,i)*Xl(:,:,i)*s.C(:,:,i)';
          end
          try
            res = pinv([XRBPB - XBP*pinv(XP)*XBP', -XCBP + XBP*pinv(XP)*XCP; -XCBP' + XCP'*pinv(XP)*XBP', JCXC - XCP'*pinv(XP)*XCP])*...
                  [XBP*pinv(XP)*PAX(:) - BPuAX(:) - BPlAX(:); - XCP'*pinv(XP)*PAX(:) + PAXuC(:) + PAXlC(:)];
          catch
            error('Could not solve for K and L')
          end
          
          L_ = res(1:s.dimX*s.dimU);
          K_ = res(s.dimX*s.dimU+1:s.dimX*(s.dimU+s.dimY));

          F_ = pinv(XP)*(PAX(:) + XBP'*L_ - XCP*K_);
         
          tL = reshape(L_,[s.dimU,s.dimX]);
          tK = reshape(K_,[s.dimX,s.dimY]);
          tF = reshape(F_,[s.dimX,s.dimX]);
          
        % END: update F, K, L
        
        % update Xu, Xl, Lu, Ll
          [Lu,Ll,oCosts] = s.RecursionLuLl(Pu,Pl,oCosts,tF,tK,tL);
          
        % END: update Xu, Xl, Lu, Ll      
      end
      
      % compute new costs
      costs = (s.InnerProduct(Lu+Ll,EXu{1}) + s.InnerProduct(Lu,EXl{1}) + s.InnerProduct(s.mode{1},oCosts))/s.simTime;

      % check convergence
        if(abs(costs-oldCosts)< s.costTol)
          fprintf('Converged after %d iterations\n',iter)
          break
        else
          oldCosts = costs;
        end      
      % END: check convergence      
      
      % DEBUG
        if(~mod(iter,20))
          disp('.')
        end
      % END: DEBUG
      
    end
    % END: iterate until convergence
    
  end % function ComputeController
  
  %% PredictXuXl
  function [Xu,Xl] = PredictXuXl(s,F,K,L) % predicts Xu and Xl for the entire optimization horizon
    Xu = cell(1,s.simTime+1);
    Xl = cell(1,s.simTime+1);      
    Xu{1} = zeros(s.dimX,s.dimX,s.numModes);
    Xl{1} = zeros(s.dimX,s.dimX,s.numModes);

    for i = 1:s.numModes
      Xu{1}(:,:,i) = s.mode{1}(i)*s.Exx0;
      Xl{1}(:,:,i) = s.mode{1}(i)*s.Ex0*s.Ex0';
    end
    
    for k = 1:s.simTime
      [Xu{k+1},Xl{k+1}] = RecursionXuXl(s,Xu{k},Xl{k},F{k},K{k},L{k},s.mode{k});
    end
  end % function PredictXuXl
   
  %% RecursionLuLl
  function [Lu_,Ll_,nCosts] = RecursionLuLl(s,Pu,Pl,oCosts,F,K,L)
    Lu_ = zeros(s.dimX,s.dimX,s.numModes);
    Ll_ = zeros(s.dimX,s.dimX,s.numModes);
    nCosts = s.OpEpsilon(oCosts,s.T);
    for i = 1:s.numModes
      Ll_(:,:,i) = L'*s.B(:,:,i)'*Pu(:,:,i)*s.B(:,:,i)*L + (F-s.B(:,:,i)*L)'*Pl(:,:,i)*(F-s.B(:,:,i)*L) + L'*s.R(:,:,i)*L;
      Lu_(:,:,i) = s.Q(:,:,i) + L'*s.R(:,:,i)*L + (s.A(:,:,i)+s.B(:,:,i)*L)'*Pu(:,:,i)*(s.A(:,:,i)+s.B(:,:,i)*L) + (s.A(:,:,i)-F-K*s.C(:,:,i)+s.B(:,:,i)*L)'*Pl(:,:,i)*(s.A(:,:,i)-F-K*s.C(:,:,i)+s.B(:,:,i)*L);
      
      nCosts(:,:,i) = nCosts(:,:,i) + trace(Pu(:,:,i)*s.H(:,:,i)*s.Eww*s.H(:,:,i)) + trace(Pl(:,:,i)*K*s.J(:,:,i)*s.Evv*s.J(:,:,i)'*K') + trace(Pl(:,:,i)*s.H(:,:,i)*s.Eww*s.H(:,:,i));

      Ll_(:,:,i) = (Ll_(:,:,i)+Ll_(:,:,i)')/2;
      Lu_(:,:,i) = (Lu_(:,:,i)+Lu_(:,:,i)')/2;
    end
  end % function RecursionLuLl

end % methods private 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

methods (Static)
  
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
  
  %% InnerProduct
  function out = InnerProduct(in1,in2)
  % computes inner product <in1,in2> = sum_i trace(in1(:,:,i)',in2(:,:,i));
  
    out = 0;
    for i = 1:size(in1,3)
      out = out + trace(in1(:,:,i)'*in2(:,:,i));
    end
  end % function InnerProduct
  
end % methods static

end % classdef