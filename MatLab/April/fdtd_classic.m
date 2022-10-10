%clear all;
%close all;
[b, fs]=audioread ('mp3.mp3');
%% Const
% tic
c=3e+10;
w=3e+8;
j0=10;
dr=0.5;
dt=dr/(c*2);
timesteps=25000;
%timesteps=100;
I=0.01*2*pi*c/w;
R1=100;
R2=200;
Rmax=1000*I;
r=(dr:dr:Rmax+dr)';
r2=(dr/2:dr/2:Rmax + dr);
NR11=fix(100/(dr/2));
NR21=fix(200/(dr/2));
Nr1=fix(Rmax/(dr/2));

NR1=fix(100/dr);
NR2=fix(200/dr);
Nr=fix(Rmax/dr);

N_delta=Nr-15*NR1;
Fr=([ones(1,NR11-1),cos(pi*((NR11:NR21)-NR11)/(2*(NR21-NR11))).^2,zeros(1,(Nr1-NR21)+2)])';
sigma_pml=([ones(1,N_delta),cos(pi*((N_delta:Nr)-N_delta)/(30*NR1)),0])';
r_alt=1./r;
a=1;
if a==0 && exist('Er')==0
    a=1;
end

%% PML
pml=1;
if pml==0
    sigma=ones(length(r),1);    
else
    sigma=sigma_pml;
end
%%
field_check_point=15*NR2;
%fig=figure('WindowState','maximized');
%% НУ
T=(dt:dt:(timesteps*dt))';
if a==1
    Er=zeros(length(r),1);
    Ep=zeros(length(r)+1,1);
    Hz=zeros(length(r),1);
    len=length(r);
    n=1;
    Ert=zeros(length(T),1);
    Ept=Ert;
    Hzt=Ert;
end
%J=zeros(length(r),1);
%% Запись файла (гифки)
%fig=figure
% frame = getframe(fig);
% [im,map] = rgb2ind(frame.cdata,4);
% imwrite(im,map,strcat('evolution','.gif'),'DelayTime',0,'Loopcount',0);
%%
for n=n:(n+timesteps)
    %J
    %J=j0.*sin(w*n*dt).*exp((-r2.^2)/(I^2));
    %J_alt (Fr)
    J=j0.*sin(w*n*dt).*Fr;
    
    %J=zeros(length(Fr),1);
    %Er
    for nr=1:length(r)
        Er(nr)=sigma(nr)*(Er(nr)+(c*dt*r_alt(nr))*Hz(nr)-(4*pi*dt)*J(2*nr-1));
    end
    %Ephi
    Ep(1)=sigma(1)*(Ep(1)+(4*pi*dt)*J(2));
    for nr=2:length(r)
        Ep(nr)=sigma(nr)*(Ep(nr)-(c*dt/dr)*(Hz(nr)-Hz(nr-1))+(4*pi*dt)*J(2*nr));
    end
    %Hz
    for nr=1:(length(r)-1)
        Hz(nr)=sigma(nr)*(Hz(nr)-(c*dt*r_alt(nr))*Er(nr)-(c*dt*r_alt(nr)/dr)*(r(nr+1)*Ep(nr+1)-r(nr)*Ep(nr)));
    end
    Hz(length(r))=sigma(length(r))*(Hz(length(r))-(c*dt*r_alt(length(r)))*Er(length(r)));
    %alter_ploting
    %ploting
    %pause(0.00000000000001)
    Ert(n)=Er(ceil((field_check_point)));
    Ept(n)=Ep(ceil((field_check_point)));
    Hzt(n)=Hz(ceil((field_check_point)));
    
%     frame = getframe(fig);     
%     fprintf(1,'Шаг № %3.0f\n',n);
%     [im,map] = rgb2ind(frame.cdata,4);
%     imwrite(im,map,strcat('evolution','.gif'),'DelayTime',0,'WriteMode','Append');   
end
%% Last Plot ans Save
ploting_classic
%alter_ploting
%saveas(gcf,strcat('n=',num2str(n-1),'.png'))
%toc
%sound (b,fs);
%fprintf(1,'Шаг № %3.0f\n',n-1);