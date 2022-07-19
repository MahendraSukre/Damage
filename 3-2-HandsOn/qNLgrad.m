function [q,Sn,En,Dn,dDdEn] = qNLgrad(q,T,X,G,u,Du,uNL,DuNL,S,E,D)
% ARF 7-JAN-2004 Based on qdbnldi (IMR)
%***************************************************
% qNLgrad:
%   Creates and assembles current internal force 
%   vector, stress matrix, strain matrix and internal
%   variable matrix for a group of 1D elements with 
%   damage model (nonlocal displacements, gradient version,
%   two fields)
% Syntax:
%   [q,Sn,En,Dn,dEdEn] = qNLgrad(q,T,X,G,u,Du,uNL,DuNL,S,E,D,dDdE)
% Input:
%   q    : current global internal force vector.
%   T    : topology matrix for elements.
%   X    : initial node coordinate matrix. 
%   G    : element property matrix.
%   u    : global local displacement vector at beginning of step.
%   Du   : increment of u.
%   uNL  : global nonlocal displacement vector at beginning of step.
%   DuNL : increment of uNL.
%   S    : current stress matrix.
%   E    : current strain matrix.
%   D    : current damage matrix.
% Output:
%   q    : updated current global internal force vector.
%   Sn   : updated stress matrix.
%   En   : updated strain matrix.
%   Dn   : updated damage matrix.
%   dDdEn: updated damage derivative matrix.
%***************************************************

% reshape global local and nonlocal displacement vectors 
U    = reshape(u,cols(X),rows(X))';
UNL  = reshape(uNL,cols(X),rows(X))';
DU   = reshape(Du,cols(X),rows(X))';
DUNL = reshape(DuNL,cols(X),rows(X))';

for j = 1:rows(T)     % Loop in elements

  % define element arrays 
  Xe    = X(T(j,1:2),:);       % Nodal coordinates
  Ue    = U(T(j,1:2),:);       % Local displacements    
  Ue    = reshape(Ue',2,1);
  DUe   = DU(T(j,1:2),:);      % Increment of local displacements    
  DUe   = reshape(DUe',2,1);
  UNLe  = UNL(T(j,1:2),:);     % Nonlocal displacements        
  UNLe  = reshape(UNLe',2,1);
  DUNLe = DUNL(T(j,1:2),:);    % Increment of nonlocal displacements    
  DUNLe = reshape(DUNLe',2,1);

  Ge   = G(T(j,3),:);           % Material properties

  % select row j and reshape into element format
  Se    = reshape(S(j,:),1,2)';     % Stresses    
  Ee    = reshape(E(j,:),1,2)';     % Strains  
  De    = reshape(D(j,:),1,2)';     % Damage

  % evaluate internal forces for element
  [qe,Sen,Een,Den,dDdEen] = qeNLgrad(Xe,Ge,Ue,DUe,UNLe,DUNLe,Se,Ee,De);  

  % assemble into global arrays
  q = assmq(q,qe,T(j,:),cols(X));     % Internal forces
  Sn(j,:)    = reshape(Sen',1,2);     % Stresses
  En(j,:)    = reshape(Een',1,2);     % Strains
  Dn(j,:)    = reshape(Den',1,2);     % Damage
  dDdEn(j,:) = reshape(dDdEen',1,2);  % Damage derivative
end

           
