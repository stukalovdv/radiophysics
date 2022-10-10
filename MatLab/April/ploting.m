subplot(2,6,[1 2])
plot(r(1:length(Er)),Er,'Color',[0 0 1],'Linewidth',1)
ylim([-1e-4 1e-4])
xlim([0 Rmax])
xlabel('r')
ylabel('Er')
title('Er')
grid minor

subplot(2,6,[3 4])
plot(r,Ep(1:length(r)),'Color',[0 0 1],'Linewidth',1)
ylim([-1e-5 1e-5])
xlim([0 Rmax])
xlabel('r')
ylabel('Ep')
title('Ep')
grid minor                

subplot(2,6,[5 6])
plot(r(1:length(Hz)),Hz,'Color',[0 0 1],'Linewidth',1)
ylim([-2e-7 1e-7])
xlim([0 Rmax])
xlabel('r')
ylabel('Hz')
title('Hz')
grid minor

subplot(2,6,[8 9])
plot(r,Jr,'Color',[0 0 1],'Linewidth',1) 
%ylim([-10 10])
xlim([0 Rmax])
xlabel('r')
ylabel('Jr')
title('Jr')
grid minor

subplot(2,6,[10 11])
plot(r,Jp,'Color',[0 0 1],'Linewidth',1) 
%ylim([-10 10])
xlim([0 Rmax])
xlabel('r')
ylabel('Jp')
title('Jp')
grid minor

sgtitle({['\nu = ',num2str(nu),';  \omega_{p0}^{2}=',num2str(wp0^2)],['R_{1} = ',num2str(NR1*dr),';   R_{2} = ',num2str(NR2*dr),';   R_{max} = ',num2str(Rmax),';'],['dr=',num2str(dr),';  dt=',num2str(dt)]})