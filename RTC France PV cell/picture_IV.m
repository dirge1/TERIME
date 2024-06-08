clear;clc;close all
[Im,Vm]=IVload;
Best_x=[0.760787966627552,0.310684576235422e-6,0.0365469455783781,52.8897869226465,1.47726933130755];
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
I1 = -Vm./(Rs + Rsh) - lambertw(Rs.*I0.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I0 + Vm)./(a.*(Rs + Rsh)))./(a.*(Rs + Rsh))).*a./Rs...
    + (Rsh.*(I0 + Iph))./(Rs + Rsh);
%% 
Best_x=[0.761192232661397,0.0198438377877356e-6,0.0653137435647713,56.5273389901431,1.31094015805559,0.999999999999905e-6,1.84427052514872];
Iph=Best_x(1);
I0=Best_x(2);
Rs=Best_x(3);
Rsh=Best_x(4);
n=Best_x(5);
I02=Best_x(6);
n2=Best_x(7);
k = 1.380649e-23;
T = 306.15;
q = 1.602176634e-19;
Vth= k*T/q;
Ns=1;
a=n*Vth*Ns;
a2=n2*Vth*Ns;
I2 = -Vm./(Rs + Rsh) - lambertw(Rs.*I0.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I0 + Vm)./(a.*(Rs + Rsh)))./(a.*(Rs + Rsh))).*a./Rs...
    - lambertw(Rs.*I02.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I02 + Vm)./(a2.*(Rs + Rsh)))./(a2.*(Rs + Rsh))).*a2./Rs + (Rsh.*(I0 + Iph + I02))./(Rs + Rsh);
%% 
Best_x=[0.761368655564523,0.999999999999486e-6,0.0809919804107505,62.0588417010747,1.97566867412311,0.00276181997840802e-6,1.19705065484424,0.999999999999834e-6,1.97566878745587];

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
I3 = -Vm./(Rs + Rsh) - lambertw(Rs.*I0.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I0 + Vm)./(a.*(Rs + Rsh)))./(a.*(Rs + Rsh))).*a./Rs...
    - lambertw(Rs.*I02.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I02 + Vm)./(a2.*(Rs + Rsh)))./(a2.*(Rs + Rsh))).*a2./Rs...
    - lambertw(Rs.*I03.*Rsh.*exp(Rsh.*(Rs.*Iph + Rs.*I03 + Vm)./(a3.*(Rs + Rsh)))./(a3.*(Rs + Rsh))).*a3./Rs+ (Rsh.*(I0 + Iph + I02 + I03))./(Rs + Rsh);
%% 
error1=abs(Im-I1);
error2=abs(Im-I2);
error3=abs(Im-I3);
figure
scatter(Vm,error1,'o')
hold on
scatter(Vm,error2,'s')
scatter(Vm,error3,'*')
grid on

ax = gca;

% 设置网格线样式为虚线
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = '--';  % 使用 '--' 代表虚线

% 添加每列的竖线
unique_voltages = unique(Vm);
for i = 1:length(unique_voltages)
    voltage = unique_voltages(i);
    % 获取当前电压值下的最大误差
    max_error = max([error1(Vm == voltage); error2(Vm == voltage); error3(Vm == voltage)]);
    % 添加竖线
    plot([voltage voltage], [0 max_error], 'Color', [0.7 0.7 0.7], 'LineStyle', '-', 'LineWidth', 0.01);
end

fosize=20;
leg=legend('SDM','DDM','TDM');
set(leg,'Location','northwest','fontsize',fosize,'fontname','Times New Roman');
set(gca, 'fontsize',fosize,'fontname','Times New Roman');
set(gcf,'unit','centimeters','position',[10 5 20 15]);
xlabel('Voltage(V)','fontsize',fosize,'fontname','Times New Roman');
ylabel('IAE(A)','fontsize',fosize,'fontname','Times New Roman');

xlim([-0.25,0.6])


%% 
figure
scatter(Vm,Im)
hold on
plot(Vm,I1)
plot(Vm,I2)
plot(Vm,I3)

grid on

ax = gca;

% 设置网格线样式为虚线
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = '--';  % 使用 '--' 代表虚线

fosize=20;
leg=legend('I-V data','SDM','DDM','TDM');
set(leg,'Location','southwest','fontsize',fosize,'fontname','Times New Roman');
set(gca, 'fontsize',fosize,'fontname','Times New Roman');
set(gcf,'unit','centimeters','position',[10 5 20 15]);
xlabel('Voltage(V)','fontsize',fosize,'fontname','Times New Roman');
ylabel('Current(A)','fontsize',fosize,'fontname','Times New Roman');

% xlim([-0.21,0.38])
ylim([0.756,0.762])
%% 
figure
scatter(Vm,Im)
hold on
plot(Vm,I1)
plot(Vm,I2)
plot(Vm,I3)
xlim([0.38,0.46])
grid on
ax = gca;
% 设置网格线样式为虚线
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = '--';  % 使用 '--' 代表虚线

fosize=20;
leg=legend('I-V data','SDM','DDM','TDM');
set(leg,'Location','southwest','fontsize',fosize,'fontname','Times New Roman');
set(gca, 'fontsize',fosize,'fontname','Times New Roman');
set(gcf,'unit','centimeters','position',[10 5 20 15]);
xlabel('Voltage(V)','fontsize',fosize,'fontname','Times New Roman');
ylabel('Current(A)','fontsize',fosize,'fontname','Times New Roman');
%% 
figure
scatter(Vm,Im)
hold on
plot(Vm,I1)
plot(Vm,I2)
plot(Vm,I3)
xlim([0.46,0.53])
grid on
ax = gca;
% 设置网格线样式为虚线
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = '--';  % 使用 '--' 代表虚线

fosize=20;
leg=legend('I-V data','SDM','DDM','TDM');
set(leg,'Location','southwest','fontsize',fosize,'fontname','Times New Roman');
set(gca, 'fontsize',fosize,'fontname','Times New Roman');
set(gcf,'unit','centimeters','position',[10 5 20 15]);
xlabel('Voltage(V)','fontsize',fosize,'fontname','Times New Roman');
ylabel('Current(A)','fontsize',fosize,'fontname','Times New Roman');
%% 
figure
scatter(Vm,Im)
hold on
plot(Vm,I1)
plot(Vm,I2)
plot(Vm,I3)
xlim([0.53,0.59])
grid on
ax = gca;
% 设置网格线样式为虚线
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = '--';  % 使用 '--' 代表虚线

fosize=20;
leg=legend('I-V data','SDM','DDM','TDM');
set(leg,'Location','southwest','fontsize',fosize,'fontname','Times New Roman');
set(gca, 'fontsize',fosize,'fontname','Times New Roman');
set(gcf,'unit','centimeters','position',[10 5 20 15]);
xlabel('Voltage(V)','fontsize',fosize,'fontname','Times New Roman');
ylabel('Current(A)','fontsize',fosize,'fontname','Times New Roman');
%% 
figure
scatter(Vm,Im)
hold on
plot(Vm,I1)
plot(Vm,I2)
plot(Vm,I3)
grid on
ax = gca;
% 设置网格线样式为虚线
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = '--';  % 使用 '--' 代表虚线

fosize=20;
leg=legend('I-V data','SDM','DDM','TDM');
set(leg,'Location','southwest','fontsize',fosize,'fontname','Times New Roman');
set(gca, 'fontsize',fosize,'fontname','Times New Roman');
set(gcf,'unit','centimeters','position',[10 5 20 15]);
xlabel('Voltage(V)','fontsize',fosize,'fontname','Times New Roman');
ylabel('Current(A)','fontsize',fosize,'fontname','Times New Roman');