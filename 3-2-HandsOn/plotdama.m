function plotdama(T,X,S,comp,nfig,step,axplot,txt,txtcolor)
     % IMR 2003 PlotElBa : Original source : poltt3
%***************************************************
% PlotElementaryVariableBar :
%   Plots a 1D line plot for any elementary variable
%   of a set of bar elements, held in the 
%   topology matrix T, using the coordinate matrix X.
% Syntax:
%   plott3(T,X,S)
% Input:
%   T         :  element topology matrix.
%   X         :  node coordinate matrix. (1D)
%   S         :  global variable held in Gauss Points.
%   comp      :  number of variable to plot
%   nfig      :  number of figure to be plotted
%   axplot    : [x1 x2 y1 y2] plot limits
%   step      :  number of step in the incremental analysis
%***************************************************

% number of elements in T
nel = rows(T) ;
% number of nodes in X
nnodes = rows(X) ;
zerr1 = zeros(nnodes);

% number of available stress/strain components
ncomp = cols(S);

if comp > (ncomp/2.)
  fprintf('\n\nRequested output component No. %d does not exist\n',comp);
  fprintf('Press any key to break plotting routine.\n\n');
  pause;
  return;
end;

% extract component S and reshape into a vector s 
 s1 = S(:,comp) ;
 s2 = S(:,comp+(ncomp/2)) ;

% vectors to bew plotted
 xplot = [] ;
 yplot = [] ;
%
% initial graphics commands
 figure(nfig) ;
 axis([0. 100. 0. 1.]) ;
 title (txt,'FontSize',14) ;
% plot contours for all elements
for i=1:nel
  xe(1:2,1) = X(T(i,1:2),1);
  xm = (xe(1) + xe(2)) / 2. ;
%
  xplot = [xplot;xe(1);xm;xm;xe(2)] ;
  yplot = [yplot;s1(i);s1(i);s2(i);s2(i)] ;
end ;
%
 if step ~= 0
 paso = sprintf('%d',step);
 text(max(X)+1e-1,s2(nel),paso) ;
 end ;
%
 plot(xplot,yplot,txtcolor,'LineWidth',2.5) ;
%
 hold on ;
