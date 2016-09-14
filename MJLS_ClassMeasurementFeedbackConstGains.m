classdef MJLS_ClassMeasurementFeedbackConstGains < handle
% MJLS_ClassMeasurementFeedbackConstGains : a reference implementation of 
% the  control law presented in M. Dolgov, U.D. Hanebeck, "Dynamic 
% Compensation of Markov Jump Linear Systems without Mode Observation", 
% 2016. 
%
% This implementation is porvided as is. For academic use only.
%
% AUTHOR: Maxim Dolgov, 20. Oct. 2015
%         maxim.dolgov@kit.edu

properties (SetAccess = protected)
  A = [];
  B = [];
  H = [];
  C = [];
  J = [];
  Q = [];
  R = [];
  Eww = [];
  Evv = [];
  T = [];
  
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
end % properties

properties(Constant)
  costTol = 1e-6;
  maxIter = 1e8;
  maxDiff = 1e-8;
end % properties constant


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

methods
  
  %% MJLS_ClassMeasurementFeedbackConstGains
  function s = MJLS_ClassMeasurementFeedbackConstGains(A,B,C,H,J,T,Q,R,Eww,Evv)
    
    s.A = A;
    s.B = B;
    s.H = H;
    s.C = C;
    s.J = J;
    s.Q = Q;
    s.R = R;
    s.Eww = Eww;
    s.Evv = Evv;
    s.T = T;
    
    [s.dimX,s.dimU,s.numModes] = size(B);
    s.dimW = size(Eww,1);
    [s.dimY,s.dimV,~] = size(J);
    
    mode_ = zeros(s.numModes,1);
    mode_(1) = 1;
    for k = 1:s.maxIter
      s.mode = mode_;
      mode_ = s.T'*s.mode;
      if(sum(abs(mode_-s.mode)) < s.maxDiff)
        s.mode = mode_;
        break
      end      
    end
    if(sum(abs(mode_-s.mode)) >= s.maxDiff)
      error('** Mode did not converge **')
    end
    
    [s.F,s.K,s.L]  = s.ComputeController();
  end % function MJLS_ClassMeasurementFeedbackConstGains
  
  %% DataInDataOut
  function [u,xe] = DataInDataOut(s,y)
    u = s.L*s.x;
    s.x = s.F*s.x + s.K*y;
    xe = s.x;
  end % function DataInDataOut
  
  %% Reset
  function Reset(s,x0)
    s.x = x0;
  end % function Reset
    
end % methods public

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

methods (Access = protected)
  
  %% ComputeController
  function [F,K,L] = ComputeController(s)
    % initialization
      Xu = zeros(s.dimX,s.dimX,s.numModes);
      Xl = zeros(s.dimX,s.dimX,s.numModes);
      Lu = zeros(s.dimX,s.dimX,s.numModes);
      Ll = zeros(s.dimX,s.dimX,s.numModes);
      for i = 1:s.numModes
        Xu(:,:,i) = 1e-4*eye(s.dimX);
        Xl(:,:,i) = 1e-4*eye(s.dimX);
        Lu(:,:,i) = 1e-4*eye(s.dimX);
        Ll(:,:,i) = 1e-4*eye(s.dimX);
      end
    % END: initialization

    
    % iterate until convergence
    
    for k = 1:s.maxIter
      Pu = s.OpEpsilon(Lu,s.T);
      Pl = s.OpEpsilon(Ll,s.T);
      
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
          JCXC = JCXC + kron(s.mode(i)*s.J(:,:,i)*s.J(:,:,i)' + s.C(:,:,i)*Xu(:,:,i)*s.C(:,:,i)' + s.C(:,:,i)*Xl(:,:,i)*s.C(:,:,i)',Pl(:,:,i));
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
        
        L = reshape(L_,[s.dimU,s.dimX]);
        K = reshape(K_,[s.dimX,s.dimY]);
        F = reshape(F_,[s.dimX,s.dimX]);
      % END: update F, K, L
      
      % update Xu, Xl, Lu, Ll
        [Xu_,Xl_] = s.RecursionXuXl(Xu,Xl,F,K,L);
        [Lu_,Ll_] = s.RecursionLuLl(Pu,Pl,F,K,L);
        
      % END: update Xu, Xl, Lu, Ll
      
      % check convergence
        if(sum(sum(sum(abs(Xu_-Xu) + abs(Xl_-Xl) + abs(Lu-Lu_) + abs(Ll_-Ll)))) < s.maxDiff)
          break
        end
        
        Xu = Xu_;
        Xl = Xl_;
        Lu = Lu_;
        Ll = Ll_;
      % END: check convergence      
      
    end
    % END: iterate until convergence
    
  end % function ComputeController
  
  %% RecursionXuXl
  function [Xu_,Xl_] = RecursionXuXl(s,Xu,Xl,F,K,L)
    Xu_ = zeros(s.dimX,s.dimX,s.numModes);
    Xl_ = zeros(s.dimX,s.dimX,s.numModes);

    for i = 1:s.numModes
      
      for j = 1:s.numModes
        Xl_(:,:,i) = Xl_(:,:,i) + s.T(j,i)*(s.mode(j)*K*s.J(:,:,j)*s.J(:,:,j)'*K' + K*s.C(:,:,j)*Xu(:,:,j)*s.C(:,:,j)'*K' + (F+K*s.C(:,:,j))*Xl(:,:,j)*(F+K*s.C(:,:,j))');
        Xu_(:,:,i) = Xu_(:,:,i) + s.T(j,i)*(s.mode(j)*s.H(:,:,j)*s.H(:,:,j)' + s.mode(j)*K*s.J(:,:,j)*s.J(:,:,j)'*K' + (s.A(:,:,j)-K*s.C(:,:,j))*Xu(:,:,j)*(s.A(:,:,j)-K*s.C(:,:,j))' + (s.A(:,:,j)-F-K*s.C(:,:,j)+s.B(:,:,j)*L)*Xl(:,:,j)*(s.A(:,:,j)-F-K*s.C(:,:,j)+s.B(:,:,j)*L)');
      end
      
      Xu_(:,:,i) = (Xu_(:,:,i)+Xu_(:,:,i)')/2;
      Xl_(:,:,i) = (Xl_(:,:,i)+Xl_(:,:,i)')/2;
    end

  end % function RecursionXuXl
  
  %% RecursionLuLl
  function [Lu_,Ll_] = RecursionLuLl(s,Pu,Pl,F,K,L)
    Lu_ = zeros(s.dimX,s.dimX,s.numModes);
    Ll_ = zeros(s.dimX,s.dimX,s.numModes);
    for i = 1:s.numModes
      Ll_(:,:,i) = L'*s.B(:,:,i)'*Pu(:,:,i)*s.B(:,:,i)*L + (F-s.B(:,:,i)*L)'*Pl(:,:,i)*(F-s.B(:,:,i)*L) + L'*s.R(:,:,i)*L;
      Lu_(:,:,i) = s.Q(:,:,i) + L'*s.R(:,:,i)*L + (s.A(:,:,i)+s.B(:,:,i)*L)'*Pu(:,:,i)*(s.A(:,:,i)+s.B(:,:,i)*L) + (s.A(:,:,i)-F-K*s.C(:,:,i)+s.B(:,:,i)*L)'*Pl(:,:,i)*(s.A(:,:,i)-F-K*s.C(:,:,i)+s.B(:,:,i)*L);

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
  
end % methods static

end % classdef