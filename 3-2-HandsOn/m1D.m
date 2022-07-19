function M = m1D(M,T,X)
%***************************************************
% m1D: 
%   Creates and assembles "mass" matrix
%   of 1D elements.
% Syntax:
%   M = m1D(M,T,X)
% Input:
%   M    :  existing global "mass" matrix.
%   T    :  topology matrix for 1D elements.
%   X    :  initial node coordinate matrix. 
%
% Output:
%   M    :  new global "mass" matrix.
%% Date:
%   ARF 27-DEC-2003
%***************************************************

for j = 1:rows(T)

  % define element coordinates
  Xe = X(T(j,1:2),:);              
             
  % evaluate element mass
  Me  = me1D(Xe);

  % assemble element mass into global mass
  M   = assmk(M,Me,T(j,:),cols(X)); 
end

           
