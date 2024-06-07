%---------------------------------------------------------------------------------------------------------------------------
% DIWNgo
% DIWNgo: Ngo driven by individual weights for enhanced photovoltaic model parameter estimation
function [Best_Ngo_rate,Best_Ngo,Convergence_curve]=HSOA(N,Max_iter,lb,ub,dim,fobj,index)
% initialize position

Best_Ngo=zeros(1,dim);
Best_Ngo_rate=inf;%change this to -inf for maximization problems
Local_best=ones(1,N)*inf;
Local_best_position=zeros(N,dim);
Ngopop=initialization(N,dim,ub,lb);%Initialize the set of random solutions
Lb=lb.*ones(1,dim);% lower boundary
Ub=ub.*ones(1,dim);% upper boundary
it=1;%Number of iterations
Convergence_curve=zeros(1,Max_iter);
Ngo_rates=zeros(1,N);%Initialize the fitness value
newNgo_rates=zeros(1,N);

%Calculate the fitness value of the initial position
for i=1:N
    x=Ngopop(i,:);
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
    T = 306.15;
    q = 1.602176634e-19;
    Vth= k*T/q;
    Ns=1;
    a=n*Vth*Ns;
    a2=n2*Vth*Ns;
    a3=n3*Vth*Ns;
    I = -Vm./(Rs + Rsh) - lambertw(Rs.*I0.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I0 + Vm)./(a.*(Rs + Rsh)))./(a.*(Rs + Rsh))).*a./Rs...
        - lambertw(Rs.*I02.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I02 + Vm)./(a2.*(Rs + Rsh)))./(a2.*(Rs + Rsh))).*a2./Rs...
        - lambertw(Rs.*I03.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I03 + Vm)./(a3.*(Rs + Rsh)))./(a3.*(Rs + Rsh))).*a3./Rs+ (Rsh.*(I0 + Iph + I02 + I03))./(Rs + Rsh);

    fit_I=sqrt(sum((Im-I).^2)/length(Im));

    Ngo_rates(1,i)=fit_I;

    %Make greedy selections
    if Ngo_rates(1,i)<Best_Ngo_rate
        Best_Ngo_rate=Ngo_rates(1,i);
        Best_Ngo=Ngopop(i,:);
    end
    if Ngo_rates(1,i)<Local_best(i)
        Local_best(i)=Ngo_rates(1,i);
        Local_best_position(i,:)=Ngopop(i,:);
    end
end
% Main loop

while it <= Max_iter
    for i=1:N
        fcmax=2;
        fcmin=0;
        F=0.5;
        u=1;
        v=1;
        kk=rand*2*pi;
        r=u*exp(kk*v);
        xx=r*cos(kk);
        yy=r*sin(kk);
        zz=r*kk;
        A=fcmax-(fcmax-fcmin)*cos(pi/2*(1-it/Max_iter));

        B=2*A^2*rand;
        Ngo_new(i,:)=Local_best_position(i,:)+abs(A*Ngopop(i,:)+B*(Best_Ngo(1,:)-Ngopop(i,:)))*xx*yy*zz;
        Flag4ub=Ngo_new(i,:)>ub;
        Flag4lb=Ngo_new(i,:)<lb;
        Ngo_new(i,:)=(Ngo_new(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

        x=Ngo_new(i,:);
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
        T = 306.15;
        q = 1.602176634e-19;
        Vth= k*T/q;
        Ns=1;
        a=n*Vth*Ns;
        a2=n2*Vth*Ns;
        a3=n3*Vth*Ns;
        I = -Vm./(Rs + Rsh) - lambertw(Rs.*I0.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I0 + Vm)./(a.*(Rs + Rsh)))./(a.*(Rs + Rsh))).*a./Rs...
            - lambertw(Rs.*I02.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I02 + Vm)./(a2.*(Rs + Rsh)))./(a2.*(Rs + Rsh))).*a2./Rs...
            - lambertw(Rs.*I03.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I03 + Vm)./(a3.*(Rs + Rsh)))./(a3.*(Rs + Rsh))).*a3./Rs+ (Rsh.*(I0 + Iph + I02 + I03))./(Rs + Rsh);


        fit_I=sqrt(sum((Im-I).^2)/length(Im));

        newNgo_rates(1,i)=fit_I;
        if newNgo_rates(1,i)<Ngo_rates(1,i)
            Ngo_rates(1,i) = newNgo_rates(1,i);
            Ngopop(i,:) = Ngo_new(i,:);

            if newNgo_rates(1,i)< Best_Ngo_rate
                Best_Ngo_rate=Ngo_rates(1,i);
                Best_Ngo=Ngopop(i,:);
            end
            if Ngo_rates(1,i)<Local_best(i)
                Local_best(i)=Ngo_rates(1,i);
                Local_best_position(i,:)=Ngopop(i,:);
            end
        end
    end

    RP=((N-(1:N))/N).^2;
    [~,sort_position]=sort(Ngo_rates);
    for i=1:N
        if rand<RP(i)
            for o=1:10000
                WW=randperm(N,3);
                if WW(1)~=i || WW(2)~=i || WW(3)~=i
                    break;
                end
            end
            X1=Ngopop(WW(1),:);
            X2=Ngopop(WW(2),:);
            X3=Ngopop(WW(3),:);
            Ngo_new(sort_position(i),:)=X1(1,:)+F*(X2(1,:)-X3(1,:));
        end
    end
    Flag4ub=Ngo_new(i,:)>ub;
    Flag4lb=Ngo_new(i,:)<lb;
    Ngo_new(i,:)=(Ngo_new(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
    x=Ngo_new(i,:);
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
    T = 306.15;
    q = 1.602176634e-19;
    Vth= k*T/q;
    Ns=1;
    a=n*Vth*Ns;
    a2=n2*Vth*Ns;
    a3=n3*Vth*Ns;
    I = -Vm./(Rs + Rsh) - lambertw(Rs.*I0.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I0 + Vm)./(a.*(Rs + Rsh)))./(a.*(Rs + Rsh))).*a./Rs...
        - lambertw(Rs.*I02.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I02 + Vm)./(a2.*(Rs + Rsh)))./(a2.*(Rs + Rsh))).*a2./Rs...
        - lambertw(Rs.*I03.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I03 + Vm)./(a3.*(Rs + Rsh)))./(a3.*(Rs + Rsh))).*a3./Rs+ (Rsh.*(I0 + Iph + I02 + I03))./(Rs + Rsh);

    fit_I=sqrt(sum((Im-I).^2)/length(Im));

    newNgo_rates(1,i)=fit_I;
    if newNgo_rates(1,i)<Ngo_rates(1,i)
        Ngo_rates(1,i) = newNgo_rates(1,i);
        Ngopop(i,:) = Ngo_new(i,:);

        if newNgo_rates(1,i)< Best_Ngo_rate
            Best_Ngo_rate=Ngo_rates(1,i);
            Best_Ngo=Ngopop(i,:);
        end
        if Ngo_rates(1,i)<Local_best(i)
            Local_best(i)=Ngo_rates(1,i);
            Local_best_position(i,:)=Ngopop(i,:);
        end
    end
    Convergence_curve(it)=Best_Ngo_rate;
    it=it+1;
end

