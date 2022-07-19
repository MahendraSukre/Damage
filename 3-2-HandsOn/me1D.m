function Me = me1D(Xe)
%***************************************************
% me1D: 
%   Creates the element "mass" matrix of
%   a 1D element.
% Syntax:
%   Me = me1D(Xe)
% Input:
%   Xe  : nodal coordinates = [x1    
%                              x2]    
% Output:
%   Me   : element "mass" matrix.
% Date:
%   ARF 27-DEC-2003
%****************************************************

he = abs(Xe(2)-Xe(1)) ;          % Length of element
Me = he*[1/3 1/6
         1/6 1/3] ;              %  Element mass matrix
     



                     
