% Force-displacement curve
clf 
figure (1) 
plot(uplot ,fplot,'r-','LineWidth',2.5)
xlabel('Displacement','FontSize',24);
ylabel('Force','FontSize',24,'VerticalAlignment','bottom');
axis square ; axis([0. 0.1 0. 2.]) ;
set(gca,'Fontsize',24)

% Final profile of damage
n=2 ;
figure (n)
plotdama(T,X,D,1,n,0,[X(1,1) X(nnodes,1) 0. 1.],' ','b-') ;
xlabel('x','FontSize',24);
ylabel('Damage','FontSize',24,'VerticalAlignment','bottom');
set(gca,'Fontsize',24)
axis square ; axis([0. 100 0. 1.2]) ;
