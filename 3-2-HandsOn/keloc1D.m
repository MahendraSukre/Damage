function Ke = keloc1D(Xe,Ge,Ee,De,dDdEe)
% ARF 5-MAR-2004 
%***************************************************
% keloc1D:
%   Creates the element "local" tangent stiffness matrix 
%   of a linear 2-node 1D element with damage model
% Syntax: 
%   Ke = keloc1D(Xe,Ge,De,dDdEe)
% Input:
%   Xe   : Coordinates Xe = [x1; x2].
%   Ge   : Element material data Ge = [E , eps0, epsf]
%   Ee   : Strain matrix Ee = [Eegp1 ; Eegp2]
%   De   : Damage matrix De = [Dgp1 ; Dgp2]
%   dDdEe: Damage derivative matrix dDdEe = [dDdEgp1 ; dDdEgp2]

% Output:
%   Ke   : Element stiffness matrix.
%***************************************************

% Gauss abscissae and weights (2-point quadrature)
r = [-1 1]/sqrt(3);
w = [ 1 1]; 

% Young modulus
youn = Ge(1);  

% determine number of nodes per element
nnodes = rows(Xe);

% Initialize stiffness matrix.
Ke = zeros(nnodes);

% Gauss integration of stiffness matrix.

for gp = 1:2  % 2-point quadrature

  % Derivatives of shape functions
  dN = [-1 1]/2; 

  % transform to global coordinates
  Jt = dN*Xe;
  B  = Jt\dN;
 
  % Gauss-point contribution to stiffness matrix
  Ke = Ke - w(gp)*(B'*youn*Ee(gp)*dDdEe(gp)*B)*det(Jt);
  
  %disp(sprintf('Value in local matrix %f',youn*Ee(gp)*dDdEe(gp) ))

end

