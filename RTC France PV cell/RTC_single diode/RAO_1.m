function [Best_Ngo_rate,Best_Ngo,Convergence_curve]=RAO_1(N,Max_iter,lb,ub,dim)
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
    k = 1.380649e-23;
    T = 306.15;
    q = 1.602176634e-19;
    Vth= k*T/q;
    Ns=1;
    a=n*Vth*Ns;
    
    I = -Vm./(Rs + Rsh) - lambertw(Rs.*I0.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I0 + Vm)./(a.*(Rs + Rsh)))./(a.*(Rs + Rsh))).*a./Rs...
        + (Rsh.*(I0 + Iph))./(Rs + Rsh);
    
    fit_I=sqrt(sum((Im-I).^2)/length(Im));
    
    Ngo_rates(1,i)=fit_I;
    
    %Make greedy selections
    if Ngo_rates(1,i)<Best_Ngo_rate
        Best_Ngo_rate=Ngo_rates(1,i);
        Best_Ngo=Ngopop(i,:);
    end
end
% Main loop

while it <= Max_iter
    [Ngo_worst,worst_position]=max(Ngo_rates);
    for i=1:N
        for j=1:dim
            Ngo_new(i,j)=Ngopop(i,j)+rand*(Best_Ngo(j)-Ngopop(worst_position,j));
        end
    end
    
    for i=1:N
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
    k = 1.380649e-23;
    T = 306.15;
    q = 1.602176634e-19;
    Vth= k*T/q;
    Ns=1;
    a=n*Vth*Ns;
    
    I = -Vm./(Rs + Rsh) - lambertw(Rs.*I0.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I0 + Vm)./(a.*(Rs + Rsh)))./(a.*(Rs + Rsh))).*a./Rs...
        + (Rsh.*(I0 + Iph))./(Rs + Rsh);
    
    fit_I=sqrt(sum((Im-I).^2)/length(Im));
    
    newNgo_rates(1,i)=fit_I;
    if newNgo_rates(1,i)<Ngo_rates(1,i)
        Ngo_rates(1,i) = newNgo_rates(1,i);
        Ngopop(i,:) = Ngo_new(i,:);
        
        if newNgo_rates(1,i)< Best_Ngo_rate
            Best_Ngo_rate=Ngo_rates(1,i);
            Best_Ngo=Ngopop(i,:);
        end
    end
    end
    
    for i=1:N
        for j=1:dim
            alpha=rand;
            beita=rand;
            index=rand;
            if index<alpha
                Ngo_new(i,j)=Ngopop(i,j)+rand*(Ngopop(i,j)-Best_Ngo(j)).*Levy(1);
            elseif index>beita
                Ngo_new(i,j)=Ngopop(i,j)+exp(rand)*(Best_Ngo(j)-Ngopop(i,j));
            else
                Ngo_new(i,j)=mean(Ngopop(:,j))+rand*(mean(Ngopop(:,j))-Best_Ngo(j));
            end
        end
    end
        for i=1:N
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
    k = 1.380649e-23;
    T = 306.15;
    q = 1.602176634e-19;
    Vth= k*T/q;
    Ns=1;
    a=n*Vth*Ns;
    
    I = -Vm./(Rs + Rsh) - lambertw(Rs.*I0.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I0 + Vm)./(a.*(Rs + Rsh)))./(a.*(Rs + Rsh))).*a./Rs...
        + (Rsh.*(I0 + Iph))./(Rs + Rsh);
    
    fit_I=sqrt(sum((Im-I).^2)/length(Im));
    
    newNgo_rates(1,i)=fit_I;
    if newNgo_rates(1,i)<Ngo_rates(1,i)
        Ngo_rates(1,i) = newNgo_rates(1,i);
        Ngopop(i,:) = Ngo_new(i,:);
        
        if newNgo_rates(1,i)< Best_Ngo_rate
            Best_Ngo_rate=Ngo_rates(1,i);
            Best_Ngo=Ngopop(i,:);
        end
    end
    end
    
    Convergence_curve(it)=Best_Ngo_rate;
    it=it+1;
end
end

% ___________________________________
function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;
end

