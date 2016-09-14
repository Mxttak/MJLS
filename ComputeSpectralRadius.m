function out = ComputeSpectralRadius(A,T)
  [dimX,~,numModes] = size(A);
  M = zeros(numModes*dimX^2);
  
  for i = 1:numModes
    M((i-1)*dimX^2+1:i*dimX^2,(i-1)*dimX^2+1:i*dimX^2) = kron(A(:,:,i),A(:,:,i));
  end
  
  M = kron(T',eye(dimX^2))*M;
  
  out = max(abs(eig(M)));
  
end % function ComputeSpectralRadius