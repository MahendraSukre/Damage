%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gradient version of the damage model based on nonlocal displacements
% Incremental-iterative solution of the nonlinear mechanical problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scalars
%         SYMBOL        MEANING 
%         dof (=1)      (Primal) Degrees of freedom per node
%         ndof          Total number of degrees of freedom
%         nnodes        Number of nodes
%         nelem         Number of elements
%         nsteps        Number of load steps
%         itermax       Maximum number of iterations per step 
%         errorf        Error in forces (norm of residual forces)
%         errorx        Error in displacements (norm of iterative correction)
%         fref          Reference force (for convergence check)
%         normu         Norm of current displacements (for convergence check) 
%         tolf          Tolerance in forces
%         tolx          Tolerance in displacements
%
% Matrices
%         SYMBOL        DIMENSIONS       MEANING
%         CL            (2,3) (*)        Constraint matrix for du 
%         CNL           (2,3) (*)        Constraint matrix for duNL  
%         C             (4,3) (*)        Constraint matrix 
%         Ksec          (ndof,ndof)      Secant stiffness matrix
%         Kloc          (ndof,ndof)      Local tangent stiffness matrix
%         M             (ndof,ndof)      Mass matrix
%         B             (ndof,ndof)      M+cK matrix (gradient regularization)
%         Kite          (2*ndof,2*ndof)  iteration stiffness matrix

% Nodal fields: displacements
%         SYMBOL        DIMENSIONS       MEANING 
%         u             (ndof,1)         Local displacement field
%         uNL           (ndof,1)         Nonlocal displacement field
%         Du            (ndof,1)         Increment of u
%         DuNL          (ndof,1)         Increment of uNL 
%         du            (ndof,1)         Iterative correction of u
%         duNL          (ndof,1)         Iterative correction of uNL
%         dutot         (2*ndof,1)       Iterative correction of u and uNL
%
%         solu          (ndof,nsteps)    History of local displacements
%         soluNL        (ndof,nsteps)    History of nonlocal displacements   
%
% Nodal fields: forces
%         SYMBOL        DIMENSIONS       MEANING 
%         p             (ndof,1)         Load vector (Femlab notation)
%         Dp            (ndof,1)         Increment of p
%         q             (ndof,1)         Internal force vector (Femlab notation)
%         r             (ndof,1)         Residual force vector 
%         reac          (2,1) *          Reaction force vector

%
% Gauss point fields (two Gauss-points per 1D linear element)
%         SYMBOL        DIMENSIONS       MEANING 
%         S             (nelem,2)        Stresses at beginning of step
%         Sn            (nelem,2)        Stresses at end of step
%         E             (nelem,2)        Strains at beginning of step
%         En            (nelem,2)        Strains at end of step    
%         D             (nelem,2)        Damage parameter at beginning of step
%         Dn            (nelem,2)        Damage parameter at end of step
%         dDdE          (nelem,2)        Derivative of damage
%
%         solS          (nelem,2,nsteps) History of stresses
%         solE          (nelem,2,nsteps) History of strains
%         solD          (nelem,2,nsteps) History of damage
%         soldDdE       (nelem,2,nsteps) History of derivative of damage
%
% (*) PARTICULAR of uniaxial tension test loaded in displacements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION
 dof = 1 ;                                  % One degree of freedom per node
 ndof = nnodes*dof ;                        % Total number of dof
 u    = zeros(ndof,1) ;                     % Local displacements    
 uNL  = zeros(ndof,1) ;                     % Nonlocal displacements
 p    = zeros(ndof,1) ;                     % External loads 
 q    = zeros(ndof,1) ;                     % Internal forces 
 r    = zeros(ndof,1) ;                     % Residual forces  
 S    = zeros(nelem,2) ;                    % Stresses 
 E    = zeros(nelem,2) ;                    % Strains
 En   = E              ;
 D    = zeros(nelem,2) ;                    % Damage
 Dn   = D              ; 
 dDdE = zeros(nelem,2) ;                    % Derivative of damage
 fplot = [0] ;                              % Load (for plotting)
 uplot = [0] ;                              % Displacement (for plotting)
 
 
% BEGINNING OF LOOP IN LOAD INCREMENTS
 for n = 1:nsteps
   disp(sprintf('Load step %d',n))
 
   % Initialization of variables related to convergence checks
   % These initialization values are not shown in convergence plots
    errorf = 1 ; errorx = 1 ; fref = 0   ; normu = 0 ;
   
   % Initialization of incremental displacements
    Du   = zeros(ndof,1) ;
    DuNL = zeros(ndof,1) ;
    
   % Initialization of iteration counter
    i = 0 ;
    
   % Variables for convergence plots
    errorxlist = [] ;
    errorflist = [] ;
    iterations = [] ;
     
   % BEGINNING OF LOOP IN ITERATIONS
    while ((errorf > tolf*fref) | (errorx > tolx*normu)) & (i <= itermax)  
      
      % 2-field iteration matrix
       Ksec = zeros(ndof) ;                  % Initialization of matrices
       Kloc = zeros(ndof) ;                 
       Ksec = ksec1D(Ksec,T,X,G,Dn) ;        % Secant stiffness matrix
       Kloc = kloc1D(Kloc,T,X,G,En,Dn,dDdE) ;% Local tangent stiffness matrix
        Kite = [Ksec   Kloc                   % 2-field iteration matrix:
               -M     B   ] ;                % matrices M and B are constant
              
      % 2-field RHS vector
       RHS = [-r;zeros(ndof,1)] ;  
    
      % 2-field constraint matrix
      % PARTICULAR for uniaxial tension test loaded in displacements
       if i==0                                    % Prediction
          % Prescribed local displacements 
          CL  = [1        1 0               
                 nnodes   1 Du_pres] ;       
          % Prescribed nonlocal displacements 
          CNL = [nnodes+1 1 0               
                 2*nnodes 1 Du_pres] ;
          blx = [0 ; Du_pres ; 0 ; Du_pres] ;   
       elseif i>0                                 % Iterative corrections
          % Restrained local displacements 
           CL  = [1       1 0               
                  nnodes  1 0] ;            
          % Restrained nonlocal displacements 
           CNL = [nnodes+1 1 0               
                  2*nnodes 1 0] ; 
           blx = zeros(4,1) ;    
       end
       C = [CL;CNL] ;                       % 2-field constraint matrix        
    
      % Boundary conditions prescribed via Lagrange multipliers
       Alx = zeros(4,2*nnodes) ;
       Alx(1,1)        = 1 ;
       Alx(2,nnodes)   = 1 ;
       Alx(3,nnodes+1) = 1 ;
       Alx(4,2*nnodes) = 1 ;
       Kite = [Kite   Alx' 
               Alx zeros(4,4)] ;
     
       RHS = [RHS;blx] ;     
       
      % Linear system of equations
       dutot = Kite\RHS ;                   % "Long" vector (2*ndof) 
       du    = dutot(1:ndof) ;              % Local displacements
       duNL  = dutot(ndof+1:2*ndof) ;       % Nonlocal displacements
      
      % Update of incremental displacements 
       Du    = Du + du ;                    % Local displacements
       DuNL  = DuNL + duNL ;                % Nonlocal displacements   
       
      % Internal forces
      q = zeros(ndof,1) ;                   % Initialization of internal forces 
      [q,Sn,En,Dn,dDdE] = qNLgrad(q,T,X,G,u,Du,uNL,DuNL,S,E,D) ;
      
      % Reactions
      reac = reaction(q,CL,dof,1) ;         % Local constraint matrix used
      reac = [reac(1,2);zeros(ndof-2,1);reac(2,2)] ; % "Long" vector of reactions
      
      % Residual forces   
       r = q - p - reac ;
  
      % Variables for convergence check
       errorx = norm(du) ;                  % Norm of iterative local displacements
       errorf = norm(r) ;                   % Norm of residual forces
       normu  = norm(u+Du) ;                % Norm of current local displacements
       fref   = norm(reac) ;                % Norm of reaction forces (*)

       disp(sprintf('   Iteration %d   ex %e    ef %e',i,(errorx/normu),(errorf/fref)))
       
      % Variables for convergence plots
       errorxlist = [errorxlist (errorx/normu)] ;
       errorflist = [errorflist (errorf/fref)] ;
       iterations = [iterations i] ;

       if i==0 ; errorx = 0.0 ; end         % Check only forces at prediction 
       
      % Update of iteration counter 
       i = i + 1 ;  
          
   % END OF LOOP IN ITERATIONS
    end 
     
   % Checking the outcome of the iteration loop
    if (i > itermax)                        % Non-convergence 
        disp(sprintf('   Convergence not attained in %d iterations\n',itermax))
        break
    else                                    % Convergence 
      % Update of state                        
       u    = u+Du ;                        % Local displacements
       uNL  = uNL+DuNL ;                    % Nonlocal displacements
       S    = Sn ;                          % Stresses
       E    = En ;                          % Strains
       D    = Dn ;                          % Damage
       %dDdE = dDdEn ;                      % Derivative of damage
       
      % Storage of state 
       solu(:,n)   = u ;                    % Local displacements
       soluNL(:,n) = uNL ;                  % Nonlocal displacements 
       solS(:,:,n) = S  ;                   % Stresses 
       solE(:,:,n) = E  ;                   % Strains
       solD(:,:,n) = D  ;                   % Damage
       soldDdE(:,:,n) = dDdE ;              % Derivative of damage 

       % Load-displacement curve 
       fplot = [fplot q(ndof)] ;
       uplot = [uplot u(ndof)] ;
       
       disp(sprintf('Convergence at step %d in %d iterations\n',n,i))
    end 
          
% END OF LOOP IN LOAD INCREMENTS
 end 

% SAVING RESULTS
 save uniaxial.mat ;
