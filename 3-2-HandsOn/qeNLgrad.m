function [qe,Se,Ee,De,dDdEe] = qeNLgrad(Xe,Ge,Ue,DUe,UNLe,DUNLe,Se,Ee,De)
% ARF 7-JAN-2004 Based on qedbnldi (IMR)
% ARF 5-MAR-2004 Computation of damage derivative
%***************************************************
% qeNLgrad:
% Damage model for 2-node element
% (nonlocal displacements, gradient version, two fields) 
%   1. Evaluates the stress and strain for the current displacements
%   2. Creates the element internal force vector
%   3. Computes dD/dE needed for tangent stiffness matrix
% Syntax: 
%   [qe,Se,Ee,De,dDdEe] = qeNLgrad(Xe,Ge,Ue,UNLe,Se,Ee,De,dDdEe)
% Input:
%   Xe   : Coordinates Xe = [x1; x2]
%   Ge   : Element material data: [E , eps0, epsf]
%   Ue   : Element nodal local displacements at beginning of step [u1; u2] 
%   DUe  : Increment of Ue 
%   UNLe : Element nodal nonlocal displacements at beginning of step [uNL1; uNL2]
%   DUNLe: Increment of UNLe 
%   Se   : Element stress vector [Sgp1;Sgp2]
%   Ee   : Element strain vector [Egp1;Egp2]               
%   De   : Element damage vector [Dgp1;Dgp2]
% Output:  
%   qe   : Element internal force vector.
%   Se   : Updated element stresses.
%   Ee   : Updated element strains.
%   De   : Updated element damage.
%   dDdEe: Updated element damage derivative.
%***************************************************

% Tolerance  
 toleps = 1.e-13 ;

% Gauss abscissae and weights (2-point quadrature)
r = [-1 1]/sqrt(3);
w = [ 1 1]; 

% Material properties
youn = Ge(1) ;  % Young modulus   
eps0 = Ge(2) ;  % Damage threshold 
epsf = Ge(3) ;  % Final strain 
aux1 = epsf/(epsf-eps0) ; 
aux2 = eps0*aux1 ;
aux3 = eps0*epsf ;

% determine number of nodes per element 
nnodes = rows(Xe);

% Initialize internal force vector
qe = zeros(nnodes,1);

% Gauss integration of internal force vector.

for gp = 1:2  % 2-point quadrature

  % Derivatives of shape functions
  dN = [-1 1]/2; 

  % transform to global coordinates
  Jt = dN*Xe;
  B  = Jt\dN; 
    
  % Evaluate local strain
  Ee(gp,:) = (B*(Ue+DUe))' ;   

  % Evaluate the nonlocal strain
  EeNL  = (B*(UNLe+DUNLe))' ;   % Updated value
  DEeNL = (B*DUNLe)'        ;   % Increment
  
  % Evaluate the "current" damage parameter
  if EeNL <= eps0               % Below damage threshold
    dam = 0.0 ;
  elseif EeNL >= epsf           % Above final strain
    dam = 1.0 ;                 % No upper bound for q calculation
  else                          % Softening branch 
    dam = aux1*(1.0-eps0/EeNL) ;
  end  
    
  % Damage cannot decrease
  dam = max([De(gp,1),dam]) ;
  
  % Check: if damage is not in [0,1], there is a problem
  if dam < 0.0 | dam > 1.0 
    disp('Problems with the damage parameter') 
    break
  end  

  % Derivative of damage w.r.t. nonlocal strain
  %disp(sprintf('   Nonlocal strain increment %f',DEeNL))
  
  dDdEgp = 0.0 ;
  if DEeNL > 0.0 ;
    Emax = aux3/(eps0*dam + epsf*(1.0-dam)) ;  
    if (abs(EeNL-Emax) <= toleps)
      dDdEgp = aux2/(EeNL)^2 ;
    end
  end
    
  % Update variables
  De(gp,1) = dam ;                           % Damage
  dDdEe(gp,1) = dDdEgp ;                     % Derivative of damage 
  Se(gp,:)  = (1-De(gp,1))*youn*Ee(gp,:) ;   % Stresses
  qe = qe + w(gp)*( B'*Se(gp,:)')*det(Jt) ;  % Internal forces
  
end


