%---------------------------------------------------------------------------------------------------------------------------
% DIWNgo
% DIWNgo: Ngo driven by individual weights for enhanced photovoltaic model parameter estimation
function [Best_Ngo_rate,Best_Ngo,Convergence_curve]=DO(N,Max_iter,lb,ub,dim,fobj,index)
% initialize position

Best_Ngo=zeros(1,dim);
Best_Ngo_rate=inf;%change this to -inf for maximization problems
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
    T = 273.15+45;
    q = 1.602176634e-19;
    Vth= k*T/q;
    Ns=36;
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
    beitat=randn(N,dim);
    alpha=rand*(1/Max_iter^2*it^2-2/Max_iter*it+1);
    b=-2*alpha;
    c=1-alpha-b;
    kk=1-rand()*(c+alpha*it^2+b*it);
    if randn<1.5
        for i=1:N



            sita=rand*2*pi-pi;
            r=1/exp(sita);
            vx=r*cos(sita);
            vy=r*sin(sita);
            y=rand;
            if y<0
                YY=0;
            else
                YY=1/(y*sqrt(2*pi))*exp(-1/2*log(y)^2);
            end
            Xs=rand([1,dim]).*(Ub-Lb)+Lb;
            Ngo_new(i,:)=Ngopop(i,:)+alpha*vx*vy*YY*(Xs-Ngopop(i,:));
        end
    else
        for i=1:N
            q = 1 / (Max_iter^2 - 2*Max_iter + 1) * it^2 - 2 / (Max_iter^2 - 2*Max_iter + 1) * it + 1 + 1 / (Max_iter^2 - 2*Max_iter + 1);

            kk = 1 - rand() * q;

            Ngo_new(i,:)=Ngopop(i,:)*kk;
        end
    end
    for i=1:N
        Flag4ub=Ngo_new(i,:)>ub;
        Flag4lb=Ngo_new(i,:)<lb;
        Ngo_new(i,:)=(Ngo_new(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
    end
    for i=1:N
        for j=1:dim
            delta =  2* it /Max_iter;
            Ngo_new(i,j)=Ngo_new(i,j)-alpha*beitat(i,j)*(mean(Ngo_new(:,j))-alpha*beitat(i,j)*Ngo_new(i,j));

        end
        for i=1:N
            Flag4ub=Ngo_new(i,:)>ub;
            Flag4lb=Ngo_new(i,:)<lb;
            Ngo_new(i,:)=(Ngo_new(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        end

        for i=1:N
            for j=1:dim
                Ngo_new(i,j)=Best_Ngo(j)+Levy(1).*alpha.*(Best_Ngo(j)-Ngo_new(i,j)*delta);
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
            T = 273.15+45;
            q = 1.602176634e-19;
            Vth= k*T/q;
            Ns=36;
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
end

% ___________________________________
function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;
end

