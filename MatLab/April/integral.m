function I=integral(Ept,Hzt,timesteps)
I=0;
N=timesteps;
    for i=2:timesteps-1
        I=I+(Ept(i).*(Hzt(i)));
    end
I=I+(Ept(1)*Hzt(1)+Ept(N)*Hzt(N))/2;
end