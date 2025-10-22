% function Fitness = calFitness(PopObj)
% % Calculate the fitness by shift-based density
% N      = size(PopObj,1);
% fmax   = max(PopObj,[],1);
% fmin   = min(PopObj,[],1);
% PopObj = (PopObj-repmat(fmin,N,1))./repmat(fmax-fmin,N,1);
% Dis    = zeros(N);
% M = size(PopObj,2);
% for i = 1 : N
%     SPopObj = min(PopObj,repmat(PopObj(i,:),N,1));
%     for j = [1:i-1,i+1:N]
%         U = PopObj(i,:)-SPopObj(j,:);
%         V = PopObj(i,:)-PopObj(j,:);
%         c = length(find(V>0));
%         C = c/M;
%         D = sum(U);
%         Dis(i,j) = C * D;
%     end
% end
% Fitness = sum(Dis,2);
% end
function Fitness = calFitness(PopObj)
% Calculate the fitness by shift-based density

N      = size(PopObj,1);
fmax   = max(PopObj,[],1);
fmin   = min(PopObj,[],1);
PopObj = (PopObj-repmat(fmin,N,1))./repmat(fmax-fmin,N,1);
Dis    = inf(N);
for i = 1 : N
    SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
    for j = [1:i-1,i+1:N]
        Dis(i,j) = norm(PopObj(i,:)-SPopObj(j,:));
    end
end
Fitness = min(Dis,[],2);
end