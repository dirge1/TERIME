%---------------------------------------------------------------------------------------------------------------------------
% RIME
% RIME: A physics-based optimization
% Website and codes of RIME:http://www.aliasgharheidari.com/RIME.html
% Hang Su, Dong Zhao, Ali Asghar Heidari, Lei Liu, Xiaoqin Zhang, Majdi Mafarja, Huiling Chen
%  Last update:  Feb 11 2023
%  e-Mail: as_heidari@ut.ac.ir, aliasghar68@gmail.com, chenhuiling.jlu@gmail.com
%
%---------------------------------------------------------------------------------------------------------------------------
%  Authors: Ali Asghar Heidari(as_heidari@ut.ac.ir, aliasghar68@gmail.com),Huiling Chen(chenhuiling.jlu@gmail.com)
%---------------------------------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------------------------------
function [Best_rime_rate,Best_rime,Convergence_curve]=CHSRIME(N,Max_iter,lb,ub,dim,fobj,index)
% disp('RIME is now tackling your problem')
% initialize position

Best_rime=zeros(1,dim);
Best_rime_rate=inf;%change this to -inf for maximization problems

Rimepop=initialization(N,dim,ub,lb);%Initialize the set of random solutions
Lb=lb.*ones(1,dim);% lower boundary
Ub=ub.*ones(1,dim);% upper boundary
it=1;%Number of iterations
Convergence_curve=zeros(1,Max_iter);
Rime_rates=zeros(1,N);%Initialize the fitness value
newRime_rates=zeros(1,N);
W = 5;%Soft-rime parameters, discussed in subsection 4.3.1 of the paper

%Calculate the fitness value of the initial position
for i=1:N
    %     Rime_rates(1,i)=fobj(Rimepop(i,:));%Calculate the fitness value for each search agent

    x=Rimepop(i,:);
    [Im,Vm]=IVload;
    Iph=x(1);
    I0=x(2);
    Rs=x(3);
    Rsh=x(4);
    n=x(5);
    I02=x(6);
    n2=x(7);
    k = 1.380649e-23;
    T = 273.15+45;
    q = 1.602176634e-19;
    Vth= k*T/q;
    Ns=36;
    a=n*Vth*Ns;
    a2=n2*Vth*Ns;
    I = -Vm./(Rs + Rsh) - lambertw(Rs.*I0.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I0 + Vm)./(a.*(Rs + Rsh)))./(a.*(Rs + Rsh))).*a./Rs...
        - lambertw(Rs.*I02.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I02 + Vm)./(a2.*(Rs + Rsh)))./(a2.*(Rs + Rsh))).*a2./Rs + (Rsh.*(I0 + Iph + I02))./(Rs + Rsh);


    fit_I=sqrt(sum((Im-I).^2)/length(Im));

    Rime_rates(1,i)=fit_I;

    %Make greedy selections
    if Rime_rates(1,i)<Best_rime_rate
        Best_rime_rate=Rime_rates(1,i);
        Best_rime=Rimepop(i,:);
    end
end
% Main loop

while it <= Max_iter
    RimeFactor = (rand-0.5)*2*cos((pi*it/(Max_iter*10)))*(1-round(it*W/Max_iter)/W);%Parameters of Eq.(3),(4),(5)
    E =sqrt(it/Max_iter);%Eq.(6)
    newRimepop = Rimepop;%Recording new populations
    normalized_rime_rates=normr(Rime_rates);%Parameters of Eq.(7)
    a=2*(1-it/Max_iter);
    [~,nn]=sort(Rime_rates);
    Xa=Rimepop(nn(1),:);
    Xb=Rimepop(nn(2),:);
    Xc=Rimepop(nn(3),:);
    Pmin=0.6;
    Pmax=1;
    P=Pmax-(Pmax-Pmin)*it/Max_iter;
    for i=1:N
        if rand>0.5
            for j=1:dim
                r1=rand();
                if r1< E
                    newRimepop(i,j)=Best_rime(1,j)+RimeFactor*((Ub(j)-Lb(j))*rand+Lb(j));%Eq.(3)
                end
            end
        else
            if rand>=P
                for j=1:dim
                    X1=Xa(j)+a*(2*rand-1)*abs(2*rand*Xa(j)-newRimepop(i,j));
                    X2=Xb(j)+a*(2*rand-1)*abs(2*rand*Xb(j)-newRimepop(i,j));
                    X3=Xc(j)+a*(2*rand-1)*abs(2*rand*Xc(j)-newRimepop(i,j));
                    newRimepop(i,j)=(X1+X2+X3)/3;
                end
            else
                ww=randperm(dim,1);
                X1=Xa(ww)+a*(2*rand-1)*abs(2*rand*Xa(ww)-newRimepop(i,ww));
                X2=Xb(ww)+a*(2*rand-1)*abs(2*rand*Xb(ww)-newRimepop(i,ww));
                X3=Xc(ww)+a*(2*rand-1)*abs(2*rand*Xc(ww)-newRimepop(i,ww));
            end
        end
        %Hard-rime puncture mechanism
        r2=rand();

        for j=1:dim
            if r2<normalized_rime_rates(i)
                newRimepop(i,j)=Best_rime(1,j);%Eq.(7)
            end
        end
    end

    for i=1:N
        %Boundary absorption
        Flag4ub=newRimepop(i,:)>ub;
        Flag4lb=newRimepop(i,:)<lb;
        newRimepop(i,:)=(newRimepop(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        x=newRimepop(i,:);
        [Im,Vm]=IVload;
        Iph=x(1);
        I0=x(2);
        Rs=x(3);
        Rsh=x(4);
        n=x(5);
        I02=x(6);
        n2=x(7);
        k = 1.380649e-23;
        T = 273.15+45;
        q = 1.602176634e-19;
        Vth= k*T/q;
        Ns=36;
        a=n*Vth*Ns;
        a2=n2*Vth*Ns;
        I = -Vm./(Rs + Rsh) - lambertw(Rs.*I0.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I0 + Vm)./(a.*(Rs + Rsh)))./(a.*(Rs + Rsh))).*a./Rs...
            - lambertw(Rs.*I02.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I02 + Vm)./(a2.*(Rs + Rsh)))./(a2.*(Rs + Rsh))).*a2./Rs + (Rsh.*(I0 + Iph + I02))./(Rs + Rsh);


        fit_I=sqrt(sum((Im-I).^2)/length(Im));
        newRime_rates(1,i)=fit_I;
        %Positive greedy selection mechanism
        if newRime_rates(1,i)<Rime_rates(1,i)
            Rime_rates(1,i) = newRime_rates(1,i);
            Rimepop(i,:) = newRimepop(i,:);
            if newRime_rates(1,i)< Best_rime_rate
                Best_rime_rate=Rime_rates(1,i);
                Best_rime=Rimepop(i,:);
            end
        end
    end

    DR=0.4*(1-it/Max_iter);

    choose_sample=randperm(N,2);
    Xm=Rimepop(choose_sample(1),:);
    Xn=Rimepop(choose_sample(2),:);
    for i=1:N
        for j=1:dim
            alpha=normrnd(0.5,0.1);
            if rand>DR
                K=1;
            else
                K=0;
            end
            Rimepop(i,j)=Rimepop(i,j)+alpha*(Xm(j)-Xn(j))*K;
        end
    end
    for i=1:N
        %Boundary absorption
        Flag4ub=Rimepop(i,:)>ub;
        Flag4lb=Rimepop(i,:)<lb;
        Rimepop(i,:)=(Rimepop(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
    end
    Convergence_curve(it)=Best_rime_rate;
    it=it+1;
end
