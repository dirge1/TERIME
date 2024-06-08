clc;clear;close all
lb=[0 0 0 0 1 0 1 0 1];
ub=[1 1e-6 0.5 100 2 1e-6 2 1e-6 2];
dim=length(lb);

nPop=20; % 种群数
Max_iter=10e4-1; % 最大迭代次数

fobj=[];
result_pa={}; % 保存Optimal results

tic
for w=1:100
    
    [Best_score,Best_x,cg_curve]=RIME_improve_boundary2(nPop,Max_iter,lb,ub,dim,fobj,2);
    
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
    I02=Best_x(6);
    n2=Best_x(7);
    I03=Best_x(8);
    n3=Best_x(9);
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
    
    result_RMSE(w)=sqrt(sum((Im-I).^2)/length(I));
    
end

obj_mean=mean(result_RMSE)
obj_min=min(result_RMSE)
obj_max=max(result_RMSE)

% save('result_nPop','result_RMSE','result_pa','Time','num','obj_mean','obj_min','obj_max')