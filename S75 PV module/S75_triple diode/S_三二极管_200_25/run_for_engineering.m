clc;clear;close all
type = 7;
[lb,ub,dim,fobj] = Engineering_Problems(type);
% rng(10000)
%% 调用算法

nPop_m=20; % 种群数
Max_iter_m=10e4-1; % 最大迭代次数
M_m=1;  %初始的优化次数
N_m=10;   %后续的优化次数
O_m=0;     %循环次数
U_m=10;   %优化初始值的选择区域

r=1;
for A=1:length(nPop_m)
    for B=1:length(Max_iter_m)
        for C=1:length(M_m)
            for D=1:length(N_m)
                for E=1:length(O_m)
                    for F=1:length(U_m)
                        nPop=nPop_m(A); % 种群数
                        Max_iter=Max_iter_m(B); % 最大迭代次数
                        M=M_m(C);  %初始的优化次数
                        N=N_m(D);   %后续的优化次数
                        O=O_m(E);    %循环次数
                        U=U_m(F);   %优化初始值的选择区域
                        (O*N*(Max_iter+1))*nPop+M*(Max_iter+1)*nPop
                        tic
                        for w=1:100
                            w
                            index = 1;
                            fit_re2=[];
                            Optimal_results={}; % 保存Optimal results
                            for i=1:M
%                                 [Best_score,Best_x,cg_curve]=RIME(nPop,Max_iter,lb,ub,dim,fobj,2);
%                                 [Best_score,Best_x,cg_curve]=RIME_cross(nPop,Max_iter,lb,ub,dim,fobj,2);
%                                 [Best_score,Best_x,cg_curve]=RIME_improve(nPop,Max_iter,lb,ub,dim,fobj,2);
                                [Best_score,Best_x,cg_curve]=DIWJAYA(nPop,Max_iter,lb,ub,dim,fobj,2);
                                Optimal_results{1,index}="RIME";
                                Optimal_results{2,index}=cg_curve;
                                Optimal_results{3,index}=Best_score;
                                Optimal_results{4,index}=Best_x;
                                Optimal_results{5,index}=toc;

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
                                T = 273.15+25;
                                q = 1.602176634e-19;
                                Vth= k*T/q;
                                Ns=36;
                                a=n*Vth*Ns;
                                a2=n2*Vth*Ns;
                                a3=n3*Vth*Ns;
                                I = -Vm./(Rs + Rsh) - lambertw(Rs.*I0.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I0 + Vm)./(a.*(Rs + Rsh)))./(a.*(Rs + Rsh))).*a./Rs...
                                    - lambertw(Rs.*I02.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I02 + Vm)./(a2.*(Rs + Rsh)))./(a2.*(Rs + Rsh))).*a2./Rs...
                                    - lambertw(Rs.*I03.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I03 + Vm)./(a3.*(Rs + Rsh)))./(a3.*(Rs + Rsh))).*a3./Rs+ (Rsh.*(I0 + Iph + I02 + I03))./(Rs + Rsh);
                                fit_re2(index)=sqrt(sum((Im-I).^2)/length(I));
                                index = index +1;
                            end
                            min(fit_re2)
                            [result_RMSE(r,w),b]=min(fit_re2);
                            result_pa{r}{1,w}=Optimal_results{1,b};
                            result_pa{r}{2,w}=Optimal_results{2,b};
                            result_pa{r}{3,w}=Optimal_results{3,b};
                            result_pa{r}{4,w}=Optimal_results{4,b};
                            result_pa{r}{5,w}=Optimal_results{5,b};
                            %                             toc
                        end
                        num(r,:)=[nPop Max_iter M N O U];
                        Time(r)=toc;
                        r=r+1;
                        r
                    end
                end
            end
        end
    end
end
obj_mean=mean(result_RMSE)
obj_min=min(result_RMSE)
obj_max=max(result_RMSE)

% save('result_nPop','result_RMSE','result_pa','Time','num','obj_mean','obj_min','obj_max')