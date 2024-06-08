clear;clc;close all
TT=1:2000:99999;
TT=[TT 99999];

load result_RIME.mat
for i=1:100
c_curve(i,:)=result_pa{1}{2,i};
end
A=mean(c_curve)';
semilogy(TT,A(TT),'-*','color',[255 127 39]/255)
hold on
c_curve=[];

load result_RIME_cross.mat
for i=1:100
c_curve(i,:)=result_pa{1}{2,i};
end
A=mean(c_curve)';
semilogy(TT,A(TT),'->','color',[126 47 142]/255)
c_curve=[];

load result_CDRIME.mat
for i=1:100
c_curve(i,:)=result_pa{1}{2,i};
end
A=mean(c_curve)';
semilogy(TT,A(TT),'-<','color',[119 172 48]/255)

c_curve=[];

load result_SRIME.mat
for i=1:100
c_curve(i,:)=result_pa{1}{2,i};
end
A=mean(c_curve)';
semilogy(TT,A(TT),'-p','color',[198 0 0]/255)
c_curve=[];

load result_RIME_differential_learning.mat
for i=1:100
c_curve(i,:)=result_pa{1}{2,i};
end
A=mean(c_curve)';
plot(TT,A(TT),'-s','color',[217 83 25]/255)

load result_RIME_improve_boundary2.mat
for i=1:100
c_curve(i,:)=result_pa{1}{2,i};
end
A=mean(c_curve)';
semilogy(TT,A(TT),'-o','color',[0 114 189]/255)

load result_DIWJAYA.mat
for i=1:100
c_curve(i,:)=result_pa{1}{2,i};
end
A=mean(c_curve)';
semilogy(TT,A(TT),'-+','color',[255 201 14]/255)
hold on
c_curve=[];

load result_DO.mat
for i=1:100
c_curve(i,:)=result_pa{1}{2,i};
end
A=mean(c_curve)';
semilogy(TT,A(TT),'-x','color',[128 128 128]/255)
c_curve=[];

load result_NGO.mat
for i=1:100
c_curve(i,:)=result_pa{1}{2,i};
end
A=mean(c_curve)';
semilogy(TT,A(TT),'-v','color',[204 204 0]/255)

c_curve=[];

load result_RAO.mat
for i=1:100
c_curve(i,:)=result_pa{1}{2,i};
end
A=mean(c_curve)';
semilogy(TT,A(TT),'-d','color',[153 0 51]/255)

legend('RIME','SLCRIME','CDRIME','SRIME','MRIME','TERIME','DIWJAYA','DO','NGO','CLRAO-1', 'NumColumns', 2)
fosize=26;
set(gcf,'unit','centimeters','position',[10 5 20 15]);
xlabel('Iteration number','FontWeight','bold','fontsize',fosize,'fontname','Times New Roman');
ylabel('Mean RMSE','FontWeight','bold','fontsize',fosize,'fontname','Times New Roman');
set(gca, 'fontsize',fosize,'fontname','Times New Roman');