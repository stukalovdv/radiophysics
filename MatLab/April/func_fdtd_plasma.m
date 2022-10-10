function [dt,T,Hzt,kpd] = func_fdtd_plasma(nu_tilda,R1_tilda,R2_tilda)
    c=3e+10;
    wp0=3e+9;

    % nu_tilda=0.1;
    % R2_tilda=0.1;
    % R1_tilda=0.09;

    nu=nu_tilda*wp0;
    R1=R1_tilda*c/wp0;
    R2=R2_tilda*c/wp0;

    % [b, fs]=audioread ('mp3_last.mp3');
    %% Const
    %tic
    j0=1;

    dr=0.1*nu_tilda*(R2-R1);
    dt=dr/(c*2);
    T_tilda=100;
    timesteps=fix(T_tilda/(dt*wp0));
    T=(dt:dt:(timesteps*dt))';

    Rmax=R2*60;
    r=(dr:dr:Rmax+dr)';
    len=length(r);
    NR1=ceil(R1/(dr));
    if NR1==0
        NR1=1;
    end
    NR2=ceil(R2/(dr));
    N_pml=45*NR2;%расстояние от R_max до начала PML


    %f(r)
    cosinus=cos(pi*((NR1:NR2)-NR1)/(2*(NR2-NR1)));
    cosinus2=cosinus.^2;
    Fr=ones(len,1);
    Fr(NR1:NR2)=cosinus2;
    Fr((NR2+1):len)=0;

    r_alt=1./r;
    r_tilda=r.*wp0/c;
    nu_tilda=nu/wp0;
    R2_tilda=R2*wp0/c;
    R1_tilda=R1*wp0/c;
    N_pml_tilda=N_pml*dr*wp0/c;
    field_check_point=NR2;
    field_check_point_tilda=field_check_point*dr*wp0/c;
    %% PML
    sigma=ones(len,1);
    cosinus_pml=cos((pi/2)*((((len-N_pml):len)-(len-N_pml))/(len-(len-N_pml))));
    sigma((len+1-length(cosinus_pml)):len)=cosinus_pml;
    %fig=figure('Name','FDTD','NumberTitle','off','WindowState','maximized');
    %% НУ
    Er=zeros(len,1);
    Ep=zeros(len,1);
    Hz=zeros(len,1);
    Jr=j0.*Fr;
    Jp=-j0.*Fr;    
    len=length(r);
    Ert=zeros(length(T),1);
    Ept=Ert;
    Hzt=Ert;
    n=1;
    %% FDTD
    
    tic
    for n=n:(n+timesteps-1)
        %Jr
        Jr=Jr.*(1-dt*nu)+(dt*wp0^2/(4*pi)).*Fr.*Er;
        %Jp
        Jp=Jp.*(1-dt*nu)+(dt*wp0^2/(4*pi)).*Fr.*Ep;
        %Er
        Er(1)=Er(2);
        for nr=2:len
            Er(nr)=sigma(nr)*(Er(nr)+(c*dt*r_alt(nr))*Hz(nr)-(4*pi*dt)*Jr(nr));
        end
        %Ephi
        %Ep(1)=Ep(1)+(4*pi*dt)*Jp(1);
        Ep(1)=Ep(2);
        for nr=2:len
            Ep(nr)=sigma(nr).*(Ep(nr)-(c*dt/dr)*(Hz(nr)-Hz(nr-1))-(4*pi*dt)*Jp(nr));
            Ep(nr)=sigma(nr).*(Ep(nr)+c*dt*(r_alt(nr)*Hz(nr)-(4*pi/c)*Jz(nr)));
        end
        %Hz
        for nr=1:(len-1)
            Hz(nr)=sigma(nr).*(Hz(nr)-(c*dt*r_alt(nr))*Er(nr)-(c*dt*r_alt(nr)/dr)*(r(nr+1)*Ep(nr+1)-r(nr)*Ep(nr)));
        end
        Hz(len)=sigma(len).*(Hz(len)-(c*dt*r_alt(len))*Er(len));

        Ert(n)=Er(ceil((field_check_point)));
        Ept(n)=Ep(ceil((field_check_point)));
        Hzt(n)=Hz(ceil((field_check_point))); 
    end
    time=toc;     
    time_full_min=time/60;
    if time_full_min<1
        time_full_sec=time_full_min*60; 
        fprintf(1,'затраченное время (сек) = %3.3f\n',time_full_sec); 
    else    
    fprintf(1,'затраченное время (мин) = %3.3f\n',time_full_min); 
    end
    %% Энергия
    Wzap=(R2^2+R1^2)*((pi*j0)/wp0)^2; % Запасенная энергия
    Wizl=c*(dt/4)*field_check_point*dr*integral(Ept,Hzt,timesteps);
    kpd=Wizl/Wzap;
end