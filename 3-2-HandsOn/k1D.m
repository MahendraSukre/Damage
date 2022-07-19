function K = k1D(K,T,X)
%***************************************************
% k1d: 
%   Creates and assembles "stiffness" matrix
%   of 1D elements.
% Syntax:
%   K = kbar(K,T,X)
% Input:
%   K    :  existing global "stiffness" matrix.
%   T    :  topology matrix for 1D elements.
%   X    :  initial node coordinate matrix. 
% Output:
%   K    :  new global "stiffness" matrix.
% Date:
%   ARF 27-DEC-2003
%***************************************************

for j = 1:rows(T)

  % define element coordinates
  Xe = X(T(j,1:2),:);              

  % evaluate element stiffness
  Ke  = ke1D(Xe);

  % assemble element stiffness into global stiffness
  K   = assmk(K,Ke,T(j,:),cols(X)); 
end

           
