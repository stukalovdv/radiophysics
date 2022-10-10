%Er
subplot(4,7,[1 2])
plot(r_tilda(1:length(Er)),Er.*wp0/(4*pi),'Color',[0 0 0],'Linewidth',1.7)
xlim([0 R2_tilda])
lim=ylim;
hold on
line([R1_tilda R1_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0 0])
line([R2_tilda R2_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0.3 0])
hold off
%xlabel('r', 'interpreter', 'latex','FontSize',12)
title('$E_{r}(r)$','interpreter', 'latex','FontSize',12)
grid minor
%legend('signal','R_{1}','R_{2}','interpreter', 'tex')

%Ephi
subplot(4,7,[3 4])
plot(r_tilda,wp0/(4*pi).*Ep(1:length(r_tilda)),'Color',[0 0 0],'Linewidth',1.7)
xlim([0 R2_tilda])
lim=ylim;
hold on
line([R1_tilda R1_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0 0])
line([R2_tilda R2_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0.3 0])
line([(R_max_tilda-N_pml_tilda) R_max_tilda-N_pml_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0.5 0])
hold off
%xlabel('r','interpreter', 'latex','FontSize',12)
title('$E_{\varphi}(r)$','interpreter', 'latex','FontSize',12)
grid minor
%legend('signal','R_{1}','R_{2}','PML','interpreter', 'tex')

%Ez
subplot(4,7,[5 6])
plot(r_tilda,wp0/(4*pi).*Ez(1:length(r_tilda)),'Color',[0 0 0],'Linewidth',1.7)
xlim([0 R2_tilda])
lim=ylim;
hold on
line([R1_tilda R1_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0 0])
line([R2_tilda R2_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0.3 0])
line([(R_max_tilda-N_pml_tilda) R_max_tilda-N_pml_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0.5 0])
hold off
%xlabel('r','interpreter', 'latex','FontSize',12)
title('$E_{z}(r)$','interpreter', 'latex','FontSize',12)
grid minor

%Hr
subplot(4,7,[8 9])
plot(r_tilda,wp0/(4*pi).*Hr(1:length(r_tilda)),'Color',[0 0 0],'Linewidth',1.7)
xlim([0 R2_tilda])
lim=ylim;
hold on
line([R1_tilda R1_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0 0])
line([R2_tilda R2_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0.3 0])
line([(R_max_tilda-N_pml_tilda) R_max_tilda-N_pml_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0.5 0])
hold off
%xlabel('r','interpreter', 'latex','FontSize',12)
title('$H_{r}(r)$','interpreter', 'latex','FontSize',12)
grid minor

%Hphi
subplot(4,7,[10 11])
plot(r_tilda,wp0/(4*pi).*Hp(1:length(r_tilda)),'Color',[0 0 0],'Linewidth',1.7)
xlim([0 R2_tilda])
lim=ylim;
hold on
line([R1_tilda R1_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0 0])
line([R2_tilda R2_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0.3 0])
line([(R_max_tilda-N_pml_tilda) R_max_tilda-N_pml_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0.5 0])
hold off
%xlabel('r','interpreter', 'latex','FontSize',12)
title('$H_{\varphi}(r)$','interpreter', 'latex','FontSize',12)
grid minor

%Hz
subplot(4,7,[12 13])
plot(r_tilda(1:length(Hz)),wp0/(4*pi).*Hz,'Color',[0 0 0],'Linewidth',1.7)
xlim([0 R2_tilda])
lim=ylim;
hold on
line([R1_tilda R1_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0 0])
line([R2_tilda R2_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0.3 0])
line([(R_max_tilda-N_pml_tilda) R_max_tilda-N_pml_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0.5 0])
hold off
%xlabel('r','interpreter', 'latex','FontSize',12)
title('$H_{z}(r)$','interpreter', 'latex','FontSize',12)
grid minor
%legend('signal','R_{1}','R_{2}','PML','interpreter', 'tex')

%Jz
subplot(4,7,[15 16])
plot(r_tilda,Jr,'Color',[0 0 0],'Linewidth',1.7) 
xlim([0 R2_tilda])
lim=ylim;
ylim([lim(1) lim(2)])
hold on
line([R1_tilda R1_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0 0])
line([R2_tilda R2_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0.3 0])
hold off
%xlabel('r','interpreter', 'latex','FontSize',12)
title('$J_{r}(r)$','interpreter', 'latex','FontSize',12)
grid minor

%Jp
subplot(4,7,[17 18])
plot(r_tilda,Jp,'Color',[0 0 0],'Linewidth',1.7) 
xlim([0 R2_tilda])
lim=ylim;
ylim([lim(1) lim(2)])
hold on
line([R1_tilda R1_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0 0])
line([R2_tilda R2_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0.3 0])
hold off
%xlabel('r','interpreter', 'latex','FontSize',12)
title('$J_{\varphi}(r)$','interpreter', 'latex','FontSize',12)
grid minor

%Jz
subplot(4,7,[19 20])
plot(r_tilda,Jz,'Color',[0 0 0],'Linewidth',1.7) 
xlim([0 R2_tilda])
lim=ylim;
ylim([lim(1) lim(2)])
hold on
line([R1_tilda R1_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0 0])
line([R2_tilda R2_tilda],[lim(1) lim(2)],'LineStyle','--','Color',[1 0.3 0])
hold off
%xlabel('r','interpreter', 'latex','FontSize',12)
title('$J_{z}(r)$','interpreter', 'latex','FontSize',12)
grid minor

%Текст с параметрами
subplot(4,7,[14 21])
set(gca,'Visible', 'off')
text(0.5,0.8,strcat('\nu/\omega_{p0}=',num2str(nu_tilda)),'Interpreter', 'tex','FontSize',18,'HorizontalAlignment','center')
text(0.5,0.6,strcat('\delta=',num2str(Delta)),'Interpreter', 'tex','FontSize',18,'HorizontalAlignment','center')
text(0.5,0.4,strcat('R_{2}=',num2str(R2_tilda)),'Interpreter', 'tex','FontSize',18,'HorizontalAlignment','center')
text(0.5,0.2,strcat('\theta=\pi/',num2str(pi/theta)),'Interpreter', 'tex','FontSize',18,'HorizontalAlignment','center')