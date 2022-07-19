function K = ksec1D(K,T,X,G,D)
% ARF 7-JAN-2004 Based on ksecbar (IMR)
%***************************************************
% ksec1D: 
%   Creates and assembles secant stiffness matrix
%   for a group of 1D elements with damage model.
% Syntax:
%   K = ksec1D(K,T,X,G,D)
% Input:
%   K    : initial global stiffness matrix.
%   T    : element topology matrix.
%   X    : node coordinate matrix. 
%   G    : material property matrix. 
%   D    : current damage matrix.
% Output:
%   K    : new global stiffness matrix.
%***************************************************

for j = 1:rows(T)               % Loop in elements

  % extract element information from global arrays
  Xe = X(T(j,1:2),:);           % Nodal coordinates 
  Ge = G(T(j,3),:);             % Material properties   

  % select row j and reshape into element format
  De = reshape(D(j,:),1,2)';    % Damage

  % evaluate element stiffness
  Ke = kesec1D(Xe,Ge,De);

  % assemble element stiffness in global stiffness
  K  = assmk(K,Ke,T(j,:),1);
  
end

           
