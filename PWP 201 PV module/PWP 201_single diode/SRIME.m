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
function [Best_rime_rate,Best_rime,Convergence_curve]=SRIME(N,Max_iter,lb,ub,dim,fobj,index)
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
Wmax=0.3;
Wmin=0.0625;
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
    k = 1.380649e-23;
    T = 273.15+45;
    q = 1.602176634e-19;
    Vth= k*T/q;
    Ns=36;
    a=n*Vth*Ns;
    %电压电流

    I = -Vm./(Rs + Rsh) - lambertw(Rs.*I0.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I0 + Vm)./(a.*(Rs + Rsh)))./(a.*(Rs + Rsh))).*a./Rs...
        + (Rsh.*(I0 + Iph))./(Rs + Rsh);

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
    %     it
    RimeFactor = (rand-0.5)*2*cos((pi*it/(Max_iter*10)))*(1-round(it*W/Max_iter)/W);%Parameters of Eq.(3),(4),(5)
    E =sqrt(it/Max_iter);%Eq.(6)
    newRimepop = Rimepop;%Recording new populations
    normalized_rime_rates=normr(Rime_rates);%Parameters of Eq.(7)
    for i=1:N
        for j=1:dim
            %Soft-rime search strategy
            r1=rand();
            if r1< E
                newRimepop(i,j)=Best_rime(1,j)+RimeFactor*((Ub(j)-Lb(j))*rand+Lb(j));%Eq.(3)
            end

            %Hard-rime puncture mechanism
            r2=rand();
            if r2<normalized_rime_rates(i)
                choose_sample=randperm(N,2);
                x1=newRimepop(choose_sample(1),:);
                x2=newRimepop(choose_sample(2),:);
                newRimepop(i,j)=Best_rime(1,j)+normalized_rime_rates(i)*(x1(j)-x2(j));%Eq.(7)
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
        k = 1.380649e-23;
        T = 273.15+45;
        q = 1.602176634e-19;
        Vth= k*T/q;
        Ns=36;
        a=n*Vth*Ns;
        %电压电流

        I = -Vm./(Rs + Rsh) - lambertw(Rs.*I0.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I0 + Vm)./(a.*(Rs + Rsh)))./(a.*(Rs + Rsh))).*a./Rs...
            + (Rsh.*(I0 + Iph))./(Rs + Rsh);

        fit_I=sqrt(sum((Im-I).^2)/length(Im));
        newRime_rates(1,i)=fit_I;
        %Positive greedy selection mechanism

        delta_f=abs(newRime_rates(1,i)-Rime_rates(1,i));
        dist=0;
        for j=1:dim
            dist=dist+abs(newRimepop(i,j)-Rimepop(i,j));
        end
        %Positive greedy selection mechanism
        if newRime_rates(1,i)<Rime_rates(1,i)
            Rime_rates(1,i) = newRime_rates(1,i);
            Rimepop(i,:) = newRimepop(i,:);
            if newRime_rates(1,i)< Best_rime_rate
                Best_rime_rate=Rime_rates(1,i);
                Best_rime=Rimepop(i,:);

            end
        elseif rand<=exp(-delta_f/dist)
            Rime_rates(1,i) = newRime_rates(1,i);
            Rimepop(i,:) = newRimepop(i,:);
        end

    end
    Convergence_curve(it)=Best_rime_rate;
    it=it+1;
end
