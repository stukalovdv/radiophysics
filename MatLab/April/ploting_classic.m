subplot(2,6,[1 2])
plot(r(1:length(Er)),Er.*w/125,'Color',[0 0 1],'Linewidth',1)
%ylim([-5e-8 5e-8])
xlabel('r')
ylabel('Er')
title('Er')
grid minor

subplot(2,6,[3 4])
plot(r,(w/125).*Ep(1:length(r)),'Color',[0 0 1],'Linewidth',1)
%ylim([-4e-7 4e-7])
xlabel('r')
ylabel('Ep')
title('Ep')
grid minor                

subplot(2,6,[5 6])
plot(r(1:length(Hz)),(w/125).*Hz,'Color',[0 0 1],'Linewidth',1)
%ylim([-3e-7 3e-7])
xlabel('r')
ylabel('Hz')
title('Hz')
grid minor

subplot(2,6,[8 9])
plot(r2,J,'Color',[0 0 1],'Linewidth',1) 
xlabel('r')
ylabel('J')
title('J')
grid minor

subplot(2,6,[10 11])
plot(r,sigma_pml(1:length(r)),'Color',[0 0 1],'Linewidth',1) 
xlabel('r')
title('pml')
grid minor

% sgtitle({['\nu = ',num2str(nu),';  \omega_{p0}^{2}=',num2str(wp0^2)],['R_{1} = ',num2str(NR1*dr),';   R_{2} = ',num2str(NR2*dr),';   R_{max} = ',num2str(Rmax),';'],['dr=',num2str(dr),';  dt=',num2str(dt)]})