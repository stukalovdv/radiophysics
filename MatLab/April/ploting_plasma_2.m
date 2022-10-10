subplot(4,7,[22 23 24])
plot(T.*wp0,wp0/(4*pi).*Hzt(1:length(T)),'Color',[0 0 0],'Linewidth',1.7) 
xlabel('t','interpreter', 'latex','FontSize',12)
title('$H_{z}(t)$','interpreter', 'latex','FontSize',12)
grid minor

subplot(4,7,[25 26 27])
plot(f_tilda(1:NFn),Y_norm(1:NFn),'Linewidth',1.5,'DisplayName',strcat('\delta=',tit),'Color','black')
xlim([0 10])
title('Спектр','interpreter', 'latex','FontSize',12)
xlabel('\omega / \omega_{p0}','Interpreter', 'tex','FontSize',12);
ylabel('$S(\omega) / S_{0}$','Interpreter', 'latex','FontSize',12);
grid minor