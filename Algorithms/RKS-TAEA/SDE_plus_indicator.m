function [SDE_fitness,SDEP_fitness, Distance] = SDE_plus_indicator(PopObj, f, p)
    N = size(PopObj,1);
    PopObj = (PopObj-repmat(min(PopObj,[],1),N,1))./repmat(max(PopObj,[],1)-min(PopObj,[],1),N,1);
    Distance = inf(N);
    for i = 1 : N
        SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
        for j = 1 : N
            if(i == j)
                continue;
            end
%             d = 1/norm(PopObj(i,:),p);
            Distance(i,j) = norm(PopObj(i,:)-SPopObj(j,:),p);
        end
    end
    SDEP_fitness = zeros(1, N);
    if f == 1
        [a,Remain] = sortrows(sort(Distance, 2),'descend');
        SDEP_fitness(Remain) = N:-1:1;
    end
    SDE_fitness = min(Distance,[],2); 
end