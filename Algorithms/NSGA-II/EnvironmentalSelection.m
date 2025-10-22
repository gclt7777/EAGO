function [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Population,N)
% The environmental selection of NSGA-II
%
% This implementation matches the behaviour of the original function-based
% PlatEMO code, making it suitable for legacy environments.

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
    Next = FrontNo < MaxFNo;

    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);

    %% Select the solutions in the last front based on their crowding distances
    Last = find(FrontNo == MaxFNo);
    if ~isempty(Last)
        [~,Rank] = sort(CrowdDis(Last),'descend');
        picks = min(length(Last), max(0, N - sum(Next)));
        if picks > 0
            Next(Last(Rank(1:picks))) = true;
        end
    end

    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end

function CrowdDis = CrowdingDistance(PopObj,FrontNo)
% Calculate the crowding distance of each solution front by front

    [N,M]    = size(PopObj);
    CrowdDis = zeros(1,N);
    Fronts   = setdiff(unique(FrontNo),inf);
    for f = 1 : numel(Fronts)
        Front = find(FrontNo == Fronts(f));
        if isempty(Front)
            continue;
        end
        Fmax = max(PopObj(Front,:),[],1);
        Fmin = min(PopObj(Front,:),[],1);
        for i = 1 : M
            [~,Rank] = sortrows(PopObj(Front,i));
            CrowdDis(Front(Rank(1)))   = inf;
            CrowdDis(Front(Rank(end))) = inf;
            span = Fmax(i) - Fmin(i);
            if span <= 0
                continue;
            end
            for j = 2 : length(Front)-1
                CrowdDis(Front(Rank(j))) = CrowdDis(Front(Rank(j))) + ...
                    (PopObj(Front(Rank(j+1)),i) - PopObj(Front(Rank(j-1)),i)) / span;
            end
        end
    end
end
