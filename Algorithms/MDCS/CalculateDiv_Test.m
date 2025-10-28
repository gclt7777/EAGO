function [divValue] = CalculateDiv_Test(PopObj,q)
   %多样性计算函数

   %获得邻居个数
   [N,M] = size(PopObj);
   %q = floor(N / length(UQpi));
   %计算每个个体的多样性值divValue
   Distance = zeros(N,N);
   divValue = zeros(1,N);
   %angle = acos(1-pdist2(PopObj,PopObj,'cosine'));
   
   for i = 1:N
       for j = 1: N
            Distance(i,j) = sqrt(sum((PopObj(i,:)-PopObj(j,:)).^2)-sum(PopObj(i,:)-PopObj(j,:)).^2/M); 
       end
       [sortDis,index] = sort(Distance(i,:),'ascend');
       %sortDis = sort(angle(i,:),'ascend');
       %divValue(i) = sortDis(2)+(1-q)*sortDis(3);
       %计算divValue值
       for j = 1:q
           divValue(i) = divValue(i) + exp(-(j-1))*sortDis(j+1);
       end
   end
end