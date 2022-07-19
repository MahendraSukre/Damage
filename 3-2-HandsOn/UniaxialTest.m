%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UNIAXIAL TENSION TEST
% Reference: 
% Rodriguez-Ferran A, Morata I, Huerta A (2005), "A new damage model based
% on non-local displacements", Int. J. Numer. Anal. Meth. Geomech., 29:473-493
% DOI: 10.1002/nag.422
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 clear all ; clf ;
 format short e 

% GEOMETRY
 long  = 100 ;                              % Length of bar
 nelem = 100 ;                              % Number of elements
 weak = 0.10 ;                              % Fraction of weakened part  
 nnodes = nelem + 1 ;                       % Number of nodes (linear 1D elements)
 h = long/nelem ;                           % Element size
 
% PARAMETER c = (l_c)^2 OF GRADIENT MODEL                                    
 c = 1 ;
 
% GLOBAL NODE COORDINATE MATRIX AND TOPOLOGY MATRIX
 X = [0:h:long]' ;                          % Nodal coordinates (uniform mesh)
 X = [0:h:long]' ;                          % Nodal coordinates (uniform mesh)
 T = [1:nnodes-1;2:nnodes;ones(1,nelem)]';  % Topology matrix 
                                            % (left node; right node; material number)
 T([nelem/2,nelem/2+1],3) = 2 ;             % Two central elements weaker (material 2)                                    
 
% MASS AND DIFFUSIVITY MATRICES FOR GRADIENT REGULARIZATION
 dof = 1 ;                                  % Degrees of freedom per node 
 M = init(nnodes,dof) ;                     % Initialization of mass matrix
 K = init(nnodes,dof) ;                     % Idem of diffusivity matrix

 M = m1D(M,T,X) ;                           % Global mass matrix 
 K = k1D(K,T,X) ;                           % Global diffusivity matrix   
 B = M +c*K ;

% MATERIAL PARAMETERS (linear softening damage model)
G = [20000 1.e-4  1.25e-2                   % Material parameters: E, eps0, epsf
     18000 1.e-4  1.25e-2]  ;               % for original and weaker materials

% NUMERICAL PARAMETERS
 nsteps  = 2000 ;                           % Number of load steps
 itermax = 50   ;                           % Maximum number of iterations per step 
 tolf    = 0.5e-7 ;                         % Tolerance in forces
 tolx    = 0.5e-7 ;                         % Tolerance in displacements
 Du_pres = 0.0001 ;                         % Increment of prescribed displacement 
 
% COMPUTATION 
 NewtonRaphson ;
 
% POSTPROCESS
 Postprocess
