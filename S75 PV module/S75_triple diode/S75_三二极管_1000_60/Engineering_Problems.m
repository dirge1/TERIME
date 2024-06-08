function [lb,ub,dim,fobj] = Engineering_Problems(type)
% type：问题类型
% 不同数字 对应 不同问题
% 比如，type = 1 ： 选择优化 Tension/compression spring design problem
% type = 2 ： 选择优化 Pressure vessel design problem
switch type
    case 7
        fobj = @ pv2;
        % France
        Isc0=4.7;
        kk=2e-3;
        G=1000;
        TT=60;
        Isc = Isc0*G/1000+kk*(TT-25);
        lb=[0 0 0 0 1 0 1 0 1];
        ub=[2*Isc 1e-6 2 5000 4 1e-6 4 1e-6 4];
        dim=length(lb);
end

    function fitness=pv2(x)
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
        for i=1:length(Vm)
            I = -Vm./(Rs + Rsh) - lambertw(Rs.*I0.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I0 + Vm)./(a.*(Rs + Rsh)))./(a.*(Rs + Rsh))).*a./Rs...
                + (Rsh.*(I0 + Iph))./(Rs + Rsh);
        end
        V_nn=Vm(1:3);
        I_nn=I(1:3);
        Im_nn=Im(1:3);
        
        V_nn2=Vm(end-4:end);
        I_nn2=I(end-4:end);
        Im_nn2=Im(end-4:end);
        
        coeforig=polyfit(I_nn,V_nn,1);
        coefnew=polyfit(Im_nn,V_nn,1);
        
        coeforig2=polyfit(V_nn2,I_nn2,1);
        coefnew2=polyfit(V_nn2,Im_nn2,1);
        
        slporig=coeforig(1);
        slpnew=coefnew(1);
        slporig2=coeforig2(1);
        slpnew2=coefnew2(1);
        
        fit_re=sqrt(sum((Im-I).^2)/length(Im));
        fit_slp=((slporig-slpnew)^2)*0.001;        
        fit_slp2=(((slporig2-slpnew2)/slpnew2)^2);  
%         fitness=fit_re;
        fitness=fit_re+fit_slp2;
%         fitness=fit_re+fit_slp+fit_slp2;
    end    
end