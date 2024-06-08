%% ��ע΢�Ź��ںţ��Ż��㷨��   Swarm-Opti
% https://mbd.pub/o/author-a2mVmGpsYw==
function mx=mutations(x,x_all,xbest,lb,ub,dim,iter,max_iter,index,Best_rime_rate,j,p)
% ���룺
% x        ����ǰ�ĸ��壬һά
% x_all   ���и���
%xbest  ���Ÿ��壬һά
% lb       �������½�
% ub      �������Ͻ�
% dim    ����ά��
% iter     ��ǰ����
% max_iter ����������
% index  ����ѡ��������
% �����
% mx     �����ĸ���

% ��������ֵĲ��� �ɸ��������޸�

lb=lb.*ones(1,dim);
ub=ub.*ones(1,dim);

switch index

    %-----------------------------------------��˹����---------------------------------------------------
    case 1
        p=0.5; % ѡ�����
        if rand>p
            mu=(lb+ub)/2; % ��˹�����ľ�ֵ�������޸�
            sigma=(ub-lb)/6; % ��˹�����ķ�������޸�
            mx=normrnd(mu,sigma);  % % ��˹����������ֵ
        else
            mx=x;
        end
        % �߽��飬��ֹ���� ������Χ
        Flag4ub=mx>ub;
        Flag4lb=mx<lb;
        mx=(mx.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

        %-----------------------------------------��˹��Ӣ����---------------------------------------------
    case 2
%         p=0.4; % ѡ�����
        mu=xbest; % �����Ž� ��Ϊ��˹�����ľ�ֵ
        if rand>p
            sigma=xbest*0.001;
            % 0.005
            r=normrnd(mu,sigma);  % % ��˹����������ֵ
            mx=r; % �����ĸ���
        else
            mx=xbest;
        end
        % �߽��飬��ֹ���� ������Χ
        %         Flag4ub=mx>ub;
        %         Flag4lb=mx<lb;
        %         mx=(mx.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        %-----------------------------------------��������---------------------------------------------------
    case 3
        p=0.5; % ѡ�����
        if rand>p
            mu=0; % �����Ž� �� ��ǰ�� ��Ϊ��˹�����ľ�ֵ
            sigma=1; % �����޸�
            r=(sigma./((x-mu).^2+sigma^2))/pi;% �����ֲ��ĸ����ܶȺ���������ֵ
            mx=xbest+r.*xbest; % �� xbest���Ž���Ϊ�ο� ��ȡ�����ĸ���
            %         mx=x+r.*x; % �������� x ��Ϊ�ο� ��ȡ�����ĸ���
        else
            mx=x;
        end
        % �߽��飬��ֹ���� ������Χ
        Flag4ub=mx>ub;
        Flag4lb=mx<lb;
        mx=(mx.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        %-----------------------------------------�������ۼƷֲ�����-------------------------------------
    case 4
        p=0.5; % ѡ�����
        if rand>p
            a=0; %
            b=1; %
            rp=randn(1,dim);%
            mx=a+b*tan(pi*(rp-1/2)); %�������ۼƷֲ�����������ֵ
        else
            mx=x;
        end
        % �߽��飬��ֹ���� ������Χ
        Flag4ub=mx>ub;
        Flag4lb=mx<lb;
        mx=(mx.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        %-----------------------------------------t�ֲ�����-------------------------------------
    case 5
        p=0.2; % ѡ�����
        if rand>p
            mx=xbest+xbest.*trnd(iter); % �� xbest ���Ž���Ϊ�ο� ��ȡ�����ĸ���
            %             mx=x+x.*trnd(iter); % �� x ��Ϊ�ο� ��ȡ�����ĸ���
        else
            mx=xbest;
        end
        % �߽��飬��ֹ���� ������Χ
        %         Flag4ub=mx>ub;
        %         Flag4lb=mx<lb;
        %         mx=(mx.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        %-----------------------------------------����Ӧt�ֲ�����-------------------------------------
    case 6
        % ʹ���㷨�ĵ������� ��Ϊt�ֲ������ɶȲ���
        p=0.5; % ѡ�����
        tf = exp((iter/max_iter)^2); % ���ɶȲ��� ��������� �ʷ���������
        % ע�⣺����ֻ�򵥾����˷����Թ�ʽ�÷��������Թ�ʽ�кܶ࣬�������ɴ���
        if rand>p
            mx=xbest+xbest.*trnd(tf); % �� xbest ���Ž���Ϊ�ο� ��ȡ�����ĸ���
            %             mx=x+x.*trnd(tf); % �� x ��Ϊ�ο� ��ȡ�����ĸ���
        else
            mx=x;
        end
        % �߽��飬��ֹ���� ������Χ
        Flag4ub=mx>ub;
        Flag4lb=mx<lb;
        mx=(mx.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        %-----------------------------------------��̬�Ʊ���-------------------------------------
        % �ο����ף�����,������,��Сƽ.�Ľ������Ż��㷨������ͬ������������ʶ[J].��������ѧ��,2022,26(10):119-129.
    case 7

        p=0.5; % ѡ�����
        Ex=xbest;               %��̬������
        En=exp(iter/max_iter);   %��̬����
        He=En/10^(-3);      %��̬�Ƴ���
        if rand>p
            E_n = normrnd(En,He); % ������̬�ֲ������
            ra= normrnd(Ex,abs(E_n)); % ������̬�����
            mx=exp(-(ra-Ex).^2/(2*E_n^2)); % ���������Ⱥ���
        else
            mx=x;
        end
        % �߽��飬��ֹ���� ������Χ
        Flag4ub=mx>ub;
        Flag4lb=mx<lb;
        mx=(mx.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

        %-----------------------------------------���ڱ���-------------------------------------
        % �ο����ף�DOI�� 10.1109/TEVC.2012.2196047
    case 8

        A=1; % �������
        AT=5; % ��������
        if mod(iter,AT)==0
            mx=x.*(1+A*(0.5-rand(1,dim)));
        else
            mx=x;
        end
        % �߽��飬��ֹ���� ������Χ
        Flag4ub=mx>ub;
        Flag4lb=mx<lb;
        mx=(mx.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

        %-----------------------------------------DE/best/1-------------------------------------
    case 9

        p=0.5; % ѡ�����
        F=rand; % ��������,rand�����������
        npop = size(x_all,1); % ��Ⱥ��
        if rand>p
            mx=xbest+F*(x_all(randi(npop),:)-x_all(randi(npop),:));
        else
            mx=x;
        end
        % �߽��飬��ֹ���� ������Χ
        Flag4ub=mx>ub;
        Flag4lb=mx<lb;
        mx=(mx.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

        %-----------------------------------------DE/rand-to-best/1-------------------------------------
    case 10

        p=0.5; % ѡ�����
        F=rand; % ��������,rand�����������
        npop = size(x_all,1); % ��Ⱥ��
        if rand>p
            mx=x+F*(xbest-x)+F*(x_all(randi(npop),:)-x_all(randi(npop),:));
        else
            mx=x;
        end
        % �߽��飬��ֹ���� ������Χ
        Flag4ub=mx>ub;
        Flag4lb=mx<lb;
        mx=(mx.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

        %-----------------------------------------DE/rand/2-------------------------------------
    case 11

        p=0.5; % ѡ�����
        F=rand; % ��������,rand�����������
        npop = size(x_all,1); % ��Ⱥ��
        if rand>p
            mx=x_all(randi(npop),:)+F*(x_all(randi(npop),:)-x_all(randi(npop),:))+F*(x_all(randi(npop),:)-x_all(randi(npop),:));
        else
            mx=x;
        end
        % �߽��飬��ֹ���� ������Χ
        Flag4ub=mx>ub;
        Flag4lb=mx<lb;
        mx=(mx.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

        %-----------------------------------------DE/best/2-------------------------------------
    case 12

        p=0.5; % ѡ�����
        F=rand; % ��������,rand�����������
        npop = size(x_all,1); % ��Ⱥ��
        if rand>p
            mx=xbest+F*(x_all(randi(npop),:)-x_all(randi(npop),:))+F*(x_all(randi(npop),:)-x_all(randi(npop),:));
        else
            mx=x;
        end
        % �߽��飬��ֹ���� ������Χ
        Flag4ub=mx>ub;
        Flag4lb=mx<lb;
        mx=(mx.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

        %-----------------------------------------�Ǿ��ȱ���-------------------------------------
        % �ο����ף�http://kns.cnki.net/kcms/detail/23.1538.TP.20211217.1715.008.html
    case 13

        %         p=1-iter/max_iter; % ѡ�����
        p=0.1;
        F=round(rand); % �����ȡ 0 �� 1
        b=3; %  bΪϵͳ������������������Ŷ��Ե������� �������̶ȣ�ȡֵһ��Ϊ 2 �� 5
        if rand>p
            if F==0
                mx=xbest+(ub-xbest)*(1-rand^(1-iter/max_iter)^b);
            else
                mx=xbest-(xbest-lb)*(1-rand^(1-iter/max_iter)^b);
            end
        else
            mx=xbest;
        end
        % �߽��飬��ֹ���� ������Χ
        %         Flag4ub=mx>ub;
        %         Flag4lb=mx<lb;
        %         mx=(mx.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

end
% ax = gca;
% set(ax,'Tag',char([100,105,115,112,40,39,20316,32773,58,...
%     83,119,97,114,109,45,79,112,116,105,39,41]));
% eval(ax.Tag)
end
