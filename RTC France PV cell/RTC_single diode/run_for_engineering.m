clc;clear;close all
lb=[0 0 0 0 1];
ub=[1 1e-6 0.5 100 2];
dim=length(lb);

nPop=20; % ��Ⱥ��
Max_iter=10e4-1; % ����������

fobj=[];
result_pa={}; % ����Optimal results
 
tic
for w=1:100

    [Best_score,Best_x,cg_curve]=TERIME(nPop,Max_iter,lb,ub,dim);
    
    result_pa{1,w}='RIME';
    result_pa{2,w}=cg_curve;
    result_pa{3,w}=Best_score;
    result_pa{4,w}=Best_x;
    result_pa{5,w}=toc;
    
    [Im,Vm]=IVload;
    Iph=Best_x(1);
    I0=Best_x(2);
    Rs=Best_x(3);
    Rsh=Best_x(4);
    n=Best_x(5);
    k = 1.380649e-23;
    T = 306.15;
    q = 1.602176634e-19;
    Vth= k*T/q;
    Ns=1;
    a=n*Vth*Ns;
    I = -Vm./(Rs + Rsh) - lambertw(Rs.*I0.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I0 + Vm)./(a.*(Rs + Rsh)))./(a.*(Rs + Rsh))).*a./Rs...
        + (Rsh.*(I0 + Iph))./(Rs + Rsh);
    
    result_RMSE(w)=sqrt(sum((Im-I).^2)/length(I));

end

obj_mean=mean(result_RMSE)
obj_min=min(result_RMSE)
obj_max=max(result_RMSE)

% save('result_nPop','result_RMSE','result_pa','Time','num','obj_mean','obj_min','obj_max')