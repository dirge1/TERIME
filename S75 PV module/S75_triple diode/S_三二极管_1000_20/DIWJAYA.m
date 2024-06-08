%---------------------------------------------------------------------------------------------------------------------------
% DIWJAYA
% DIWJAYA: JAYA driven by individual weights for enhanced photovoltaic model parameter estimation
function [Best_jaya_rate,Best_jaya,Convergence_curve]=DIWJAYA(N,Max_iter,lb,ub,dim,fobj,index)
% initialize position

Best_jaya=zeros(1,dim);
Best_jaya_rate=inf;%change this to -inf for maximization problems

Jayapop=initialization(N,dim,ub,lb);%Initialize the set of random solutions
Lb=lb.*ones(1,dim);% lower boundary
Ub=ub.*ones(1,dim);% upper boundary
it=1;%Number of iterations
Convergence_curve=zeros(1,Max_iter);
Jaya_rates=zeros(1,N);%Initialize the fitness value
newJaya_rates=zeros(1,N);

%Calculate the fitness value of the initial position
for i=1:N
    x=Jayapop(i,:);
    [Im,Vm]=IVload;
    Iph=x(1);
    I0=x(2);
    Rs=x(3);
    Rsh=x(4);
    n=x(5);
    I02=x(6);
    n2=x(7);
    I03=x(8);
    n3=x(9);
    k = 1.380649e-23;
    T = 273.15+20;
    q = 1.602176634e-19;
    Vth= k*T/q;
    Ns=36;
    a=n*Vth*Ns;
    a2=n2*Vth*Ns;
    a3=n3*Vth*Ns;
    I = -Vm./(Rs + Rsh) - lambertw(Rs.*I0.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I0 + Vm)./(a.*(Rs + Rsh)))./(a.*(Rs + Rsh))).*a./Rs...
        - lambertw(Rs.*I02.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I02 + Vm)./(a2.*(Rs + Rsh)))./(a2.*(Rs + Rsh))).*a2./Rs...
        - lambertw(Rs.*I03.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I03 + Vm)./(a3.*(Rs + Rsh)))./(a3.*(Rs + Rsh))).*a3./Rs+ (Rsh.*(I0 + Iph + I02 + I03))./(Rs + Rsh);

    fit_I=sqrt(sum((Im-I).^2)/length(Im));

    Jaya_rates(1,i)=fit_I;

    %Make greedy selections
    if Jaya_rates(1,i)<Best_jaya_rate
        Best_jaya_rate=Jaya_rates(1,i);
        Best_jaya=Jayapop(i,:);
    end
end
% Main loop

while it <= Max_iter

    [Java_worst,worst_position]=max(Jaya_rates);

    for i=1:N
        if Jaya_rates(1,i)~=0
            fai(1,i)=Best_jaya_rate/Jaya_rates(1,i);
        else
            fai(1,i)=1;
        end
    end
    mean_fai=mean(fai);


    for i=1:N

        if Jaya_rates(1,i)==Java_worst
            for j=1:dim
                if rand>0.5
                    Jaya_new(i,j)=Best_jaya(1,j)+(Ub(j)-Lb(j))*normrnd(0,1-it/Max_iter);
                else
                    Jaya_new(i,j)=Best_jaya(1,j)-(Ub(j)-Lb(j))*normrnd(0,1-it/Max_iter);
                end
            end
        elseif fai(1,i)<=mean_fai^2
            for j=1:dim
                c1=rand;
                Jaya_new(i,j)=Jayapop(i,j)+rand*(c1*Best_jaya(1,j)+(1-c1)*mean(Jayapop(:,j))-abs(Jayapop(i,j)))-fai(1,i)^2*rand*(Jayapop(worst_position,j)-abs(Jayapop(i,j)));
            end
        else
            good_position=find(fai>mean_fai^2);
            Jaya_good=Jayapop(good_position,:);
            choose_sample=randperm(size(Jaya_good,1),2);
            if fai(choose_sample(1))>fai(choose_sample(2))
                xm=Jaya_good(choose_sample(1),:);
                xn=Jaya_good(choose_sample(2),:);
            else
                xm=Jaya_good(choose_sample(2),:);
                xn=Jaya_good(choose_sample(1),:);
            end
            Jaya_new(i,:)=Jayapop(i,:)+rand*(xm-xn);
        end
    end
    for i=1:N
        for j=1:dim
            if Jaya_new(i,j)>ub(j)||Jaya_new(i,j)<lb(j)
                Jaya_new(i,j)=(ub(j)+lb(j))/2+(ub(j)-lb(j))/2*(2*rand-1);
            end
        end
    end
    for i=1:N
        x=Jaya_new(i,:);
        [Im,Vm]=IVload;
        Iph=x(1);
        I0=x(2);
        Rs=x(3);
        Rsh=x(4);
        n=x(5);
        I02=x(6);
        n2=x(7);
        I03=x(8);
        n3=x(9);
        k = 1.380649e-23;
        T = 273.15+20;
        q = 1.602176634e-19;
        Vth= k*T/q;
        Ns=36;
        a=n*Vth*Ns;
        a2=n2*Vth*Ns;
        a3=n3*Vth*Ns;
        I = -Vm./(Rs + Rsh) - lambertw(Rs.*I0.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I0 + Vm)./(a.*(Rs + Rsh)))./(a.*(Rs + Rsh))).*a./Rs...
            - lambertw(Rs.*I02.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I02 + Vm)./(a2.*(Rs + Rsh)))./(a2.*(Rs + Rsh))).*a2./Rs...
            - lambertw(Rs.*I03.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I03 + Vm)./(a3.*(Rs + Rsh)))./(a3.*(Rs + Rsh))).*a3./Rs+ (Rsh.*(I0 + Iph + I02 + I03))./(Rs + Rsh);

        fit_I=sqrt(sum((Im-I).^2)/length(Im));

        newJaya_rates(1,i)=fit_I;


        if newJaya_rates(1,i)<Jaya_rates(1,i)
            Jaya_rates(1,i) = newJaya_rates(1,i);
            Jayapop(i,:) = Jaya_new(i,:);
            if newJaya_rates(1,i)< Best_jaya_rate
                Best_jaya_rate=Jaya_rates(1,i);
                Best_jaya=Jayapop(i,:);
            end
        end
    end
    Convergence_curve(it)=Best_jaya_rate;
    it=it+1;
end

