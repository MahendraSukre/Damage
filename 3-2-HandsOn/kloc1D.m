function K = kloc1D(K,T,X,G,E,D,dDdE)
% ARF 5-MAR-2004
%***************************************************
% kloc1D: 
%   Creates and assembles "local" stiffness matrix
%   for a group of 1D elements with damage model.
% Syntax:
%   K = kloc1D(K,T,X,G,E,D,dDdE)
% Input:
%   K    : initial global stiffness matrix.
%   T    : element topology matrix.
%   X    : node coordinate matrix. 
%   G    : material property matrix.
%   E    : current strain matrix.
%   D    : current damage matrix.
%   dDdE : current damage derivative matrix.
% Output:
%   K    : new global stiffness matrix.
%***************************************************

for j = 1:rows(T)               % Loop in elements

  % extract element information from global arrays
  Xe = X(T(j,1:2),:);           % Nodal coordinates 
  Ge = G(T(j,3),:);             % Material properties   

  % select row j and reshape into element format
  Ee    = reshape(E(j,:),1,2)' ;   % Strain
  De    = reshape(D(j,:),1,2)'  ;  % Damage
  dDdEe = reshape(dDdE(j,:),1,2)'; % Damage derivative
 
  % evaluate element stiffness
  Ke = keloc1D(Xe,Ge,Ee,De,dDdEe) ;

  % assemble element stiffness in global stiffness
  K  = assmk(K,Ke,T(j,:),1);
  
end

           
