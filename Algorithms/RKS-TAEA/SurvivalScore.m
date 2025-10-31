function [CrowdDis, p, normalization] = SurvivalScore(front, IdealPoint)

    [m,n] = size(front);
    CrowdDis = zeros(1,m) ;
    
    if m<n
        p = 1;
        normalization = max(front,[],1)';
        return 
    end
    
    % shift the ideal point to the origin
    front = front - IdealPoint;

    % Detect the extreme points and normalize the front
    Extreme = FindCornerSolutions(front);
    [front, normalization] = Normalize(front, Extreme);
    
    % set the distance for the extreme solutions
    CrowdDis(Extreme) = Inf;
    selected = false(1,m);
    selected(Extreme) = true;
    
    % approximate p (norm)
    d = Point2LineDistance(front, zeros(1,n), ones(1,n));
    d(Extreme) = Inf;
    [~, index] = min(d);
    %selected(index) = true;
    %CrowdDis(index) = Inf;
    p = log(n) / log(1 / mean(front(index,:)));
    
    if(isnan(p) || p<=0.1)
        p=1;
    end
    
    nn = vecnorm(front, p, 2);
    %distances = pdist2(front, front, 'minkowski', p);
    distances = zeros(m);%用SPEA的拥挤距离优化器改改？
    for i = 1 : m
        SPopObj = max(front,repmat(front(i,:),m,1));
        for j = 1 : m
            if(i == j)
                continue;
            end
            distances(i,j) = norm(front(i,:)-SPopObj(j,:),p);
        end
    end
    
    distances = distances ./ repmat(nn, 1, m);
 
    neighbors = 2;%可改动根据目标数量来设置邻居数量
    remaining = 1:m;
    remaining = remaining(~selected);
    for i=1:m-sum(selected)-1
%         maxim = mink(distances(remaining, selected),neighbors,2);
%         [d, index] = max(sum(maxim,2));
%         best = remaining(index);
%         remaining(index) = [];
%         selected(best)=true;
%         CrowdDis(1,best) = d;
        maxim = mink(distances(remaining, selected),neighbors,2);
        d = sum(maxim,2);
        CrowdDis1 = ConvergenceScore(front(remaining,:), p);
        [c,index] = max(d.*CrowdDis1');
        best = remaining(index);
        remaining(index) = [];
        selected(best)=true;
        CrowdDis(1,best) = c;
    end
    
end

function [front, normalization] = Normalize(front, Extreme)
[m,n] = size(front);

if(length(Extreme)~=length(unique(Extreme)))
    normalization = max(front,[],1)';
    front = front./repmat(normalization',m,1);
    return
end
    
% Calculate the intercepts of the hyperplane constructed by the extreme
% points and the axes
Hyperplane = front(Extreme,:)\ones(n,1);
if any(isnan(Hyperplane)) || any(isinf(Hyperplane)) || any(Hyperplane<0)
     normalization = max(front,[],1)';
else
    normalization = 1./Hyperplane;
    if any(isnan(normalization)) || any(isinf(normalization))
        normalization = max(front,[],1)';
    end
end
% Normalization
front = front./repmat(normalization',m,1);
end