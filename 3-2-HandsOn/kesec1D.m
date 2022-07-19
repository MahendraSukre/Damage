function Ke = kesec1D(Xe,Ge,De)
% ARF 7-JAN-2004 Based on kesecbar (IMR) 
%***************************************************
% kesec1D:
%   Creates the element secant stiffness matrix 
%   of a linear 2-node 1D element with damage model
% Syntax: 
%   Ke = kesec1D(Xe,Ge,De)
% Input:
%   Xe   : Coordinates Xe = [x1; x2].
%   Ge   : Element material data Ge = [E , eps0, epsf]
%   De   : Damage matrix De = [Dgp1 ; Dgp2]
%          Maximum damage allowed: 0.99999
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

  % Maximum damage value allowed = 0.99999
  damgp = De(gp,1) ;
  %disp(sprintf('   Damage %f',damgp))

  damgp = min(damgp,0.99999) ;
  % AGHHH! ARF / 24-MAR-2004
  % youn  = (1.0 - damgp)*youn ;

  % Gauss-point contribution to stiffness matrix
  Ke = Ke + w(gp)*( B'*(1.0 - damgp)*youn*B )*det(Jt);
  
  %disp(sprintf('Value in secant matrix %d %f',gp,(1.0 - damgp)*youn))
end

