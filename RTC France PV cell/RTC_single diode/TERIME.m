function [Best_rime_rate,Best_rime,Convergence_curve]=TERIME(N,Max_iter,lb,ub,dim)
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
    T = 306.15;
    q = 1.602176634e-19;
    Vth= k*T/q;
    Ns=1;
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
    if rand>0.5
        for i=1:N
            for j=1:dim
                %Soft-rime search strategy
                r1=rand();
                if r1< E
                    newRimepop(i,j)=Best_rime(1,j)+RimeFactor*((Ub(j)-Lb(j))*rand+Lb(j));%Eq.(3)
                end
            end
        end
    else
        for i=1:N
            newRimepop(i,:)=newRimepop(i,:)+rand*(newRimepop(randperm(N,1),:)-newRimepop(randperm(N,1),:));
        end
    end
    for i=1:N
        for j=1:dim
            %Hard-rime puncture mechanism
            r2=rand();
            if r2<normalized_rime_rates(i)
                r3=rand();
                if r3<0.5
                    C=(cos(pi*(it/Max_iter))+1)*(1-it/Max_iter)/2;
                    newRimepop(i,j)=newRimepop(i,j)+C*(newRimepop(randi(N),j)-newRimepop(randi(N),j));
                else
                    newRimepop(i,j)=normrnd(Best_rime(1,j),Best_rime(1,j)*0.001);
                end
            end
        end
    end
    for i=1:N
        for j=1:dim
            %Boundary absorption
            if newRimepop(i,j)>ub(j) || newRimepop(i,j)<lb(j)
                newRimepop(i,j)=(ub(j)+lb(j))/2+(ub(j)-lb(j))/2*(2*rand-1);
            end
        end
    end
    for i=1:N
        x=newRimepop(i,:);
        [Im,Vm]=IVload;
        Iph=x(1);
        I0=x(2);
        Rs=x(3);
        Rsh=x(4);
        n=x(5);
        k = 1.380649e-23;
        T = 306.15;
        q = 1.602176634e-19;
        Vth= k*T/q;
        Ns=1;
        a=n*Vth*Ns;

        I = -Vm./(Rs + Rsh) - lambertw(Rs.*I0.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I0 + Vm)./(a.*(Rs + Rsh)))./(a.*(Rs + Rsh))).*a./Rs...
            + (Rsh.*(I0 + Iph))./(Rs + Rsh);
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
    Convergence_curve(it)=Best_rime_rate;
    it=it+1;
end
