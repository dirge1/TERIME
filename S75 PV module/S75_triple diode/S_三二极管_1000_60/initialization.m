% This function initialize the first population of search agents
function Positions=initialization(SearchAgents_no,dim,ub,lb)

Boundary_no= size(ub,2); % number of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end

% chaotic_map = @(x) 4*(1-x)+1/x; % Logistic映射
% chaotic_sequence = zeros(SearchAgents_no, dim);
% for i = 1:SearchAgents_no
%     for j = 1:dim
%         if i == 1 && j == 1
%             chaotic_sequence(i,j) = rand(); % 初始值为随机数
%             K=chaotic_sequence(i,j);
%         else
%             chaotic_sequence(i,j) = mod(chaotic_map(K),1);
%             K=chaotic_sequence(i,j);
%         end
%     end
% end

% chaotic_map = @(x) 4*(1-x)*x; % Logistic映射
% chaotic_sequence = zeros(SearchAgents_no, dim);
% for i = 1:SearchAgents_no
%     for j = 1:dim
%         if i == 1 && j == 1
%             chaotic_sequence(i,j) = rand(); % 初始值为随机数
%             K=chaotic_sequence(i,j);
%         else
%             chaotic_sequence(i,j) = chaotic_map(K);
%             K=chaotic_sequence(i,j);
%         end
%     end
% end

% If each variable has a different lb and ub
% if Boundary_no>1
%     for i=1:dim
%         ub_i=ub(i);
%         lb_i=lb(i);
%         Positions(:,i)=chaotic_sequence(:,i).*(ub_i-lb_i)+lb_i;
%     end
% end

if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end