function Ke = ke1D(Xe)
%***************************************************
% ke1D: 
%   Creates the element "stiffness" matrix of
%   a 1D element.
% Syntax:
%   Ke = ke1D(Xe)
% Input:
%   Xe  : nodal coordinates = [x1    
%                              x2]    
% Output:
%   Ke   : element "stiffness" matrix.
% Date:
%   ARF 27-DEC-2003
%***************************************************

he = abs(Xe(2)-Xe(1)) ;          % Length of element
Ke = (1/he)*[ 1 -1 
             -1  1] ;            % Element stiffness matrix
         
