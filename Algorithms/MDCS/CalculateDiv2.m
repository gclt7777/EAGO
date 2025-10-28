function [Distance] = CalculateDiv2(PopObj,N)
   %Caculate parallel distance
   [~,M] = size(PopObj);
   PopObj = (PopObj-repmat(min(PopObj),N,1))./(repmat(max(PopObj)-min(PopObj),N,1));

   Distance = zeros(N,N);
   for i = 1:N
       Fi = PopObj(i,:);
       Fdelta = PopObj - repmat(Fi,N,1);
       Distance(i,:) = sqrt(sum(Fdelta.^2,2)-(sum(Fdelta,2)).^2./M);
       Distance(i,i) = Inf;
   end

   % Distance = inf(N);
   % for i = 1 : N
   %     SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
   %     for j = [1:i-1,i+1:N]
   %         Distance(i,j) = norm(PopObj(i,:)-SPopObj(j,:));
   %     end
   % end
end