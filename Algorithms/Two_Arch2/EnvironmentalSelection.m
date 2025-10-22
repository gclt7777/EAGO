function [CA,DA,FrontNo,CrowdDis] = EnvironmentalSelection(Population,N)
% The environmental selection procedure of Two_Archive2
%
% This version retains the structure of the original function-based
% implementation so that historical PlatEMO releases can execute the
% algorithm without relying on class-based helpers.

    if isempty(Population)
        CA = Population;
        DA = Population;
        FrontNo  = [];
        CrowdDis = [];
        return;
    end

    %% Update the convergence archive (CA)
    CA = UpdateCA(Population,N);
    if isempty(CA)
        FrontNo  = [];
        CrowdDis = [];
    else
        [FrontNo,~] = NDSort(CA.objs,CA.cons,N);
        CrowdDis    = CrowdingDistance(CA.objs,FrontNo);
    end

    %% Update the diversity archive (DA)
    DA = UpdateDA(Population,N);
end

function CA = UpdateCA(Population,N)
% Update the convergence archive by non-dominated sorting and crowding

    if isempty(Population)
        CA = Population;
        return;
    end
    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
    Next = FrontNo < MaxFNo;
    CA   = Population(Next);
    if numel(CA) < N
        Last = find(FrontNo == MaxFNo);
        if ~isempty(Last)
            CrowdDis = CrowdingDistance(Population.objs,FrontNo);
            [~,rank] = sort(CrowdDis(Last),'descend');
            picked   = Last(rank(1:min(N-numel(CA),numel(Last))));
            CA       = [CA,Population(picked)];
        end
    elseif numel(CA) > N
        CrowdDis = CrowdingDistance(Population.objs,FrontNo);
        chosen   = find(Next);
        [~,rank] = sort(CrowdDis(chosen),'descend');
        CA       = Population(chosen(rank(1:N)));
    end
end

function DA = UpdateDA(Population,N)
% Update the diversity archive using reference vector niching

    PopObj = Population.objs;
    [NP,M] = size(PopObj);
    [W,~]  = UniformPoint(N,M);

    % Normalise objective values
    fmin   = min(PopObj,[],1);
    fmax   = max(PopObj,[],1);
    range  = max(fmax - fmin, 1e-12);
    NormObj = (PopObj - repmat(fmin,NP,1))./repmat(range,NP,1);
    NormObj(~isfinite(NormObj)) = 0;

    % Associate each solution to a reference vector
    Cosine = 1 - pdist2(NormObj,W,'cosine');
    Cosine(isnan(Cosine)) = 0;
    Cosine = max(min(Cosine,1),-1);
    Angle  = acos(Cosine);
    [~,associate] = min(Angle,[],2);
    Distance = sqrt(sum(NormObj.^2,2));

    Selected = false(1,NP);
    for i = 1 : size(W,1)
        current = find(associate==i);
        if isempty(current)
            continue;
        end
        [~,best] = max(Distance(current));
        Selected(current(best)) = true;
    end

    chosen = find(Selected);
    if numel(chosen) < N
        remain = setdiff(1:NP,chosen);
        if ~isempty(remain)
            [~,rank] = sort(Distance(remain),'descend');
            fill = remain(rank(1:min(N-numel(chosen),numel(remain))));
            chosen = [chosen,fill]; %#ok<AGROW>
        end
    elseif numel(chosen) > N
        [~,rank] = sort(Distance(chosen),'descend');
        chosen = chosen(rank(1:N));
    end
    DA = Population(chosen);
end

function CrowdDis = CrowdingDistance(PopObj,FrontNo)
% Calculate crowding distances per non-dominated front

    [N, M] = size(PopObj);
    CrowdDis = zeros(1,N);
    Fronts   = setdiff(unique(FrontNo),inf);
    for f = 1 : numel(Fronts)
        front = find(FrontNo == Fronts(f));
        if numel(front) < 2
            CrowdDis(front) = inf;
            continue;
        end
        Fmax = max(PopObj(front,:),[],1);
        Fmin = min(PopObj(front,:),[],1);
        span = Fmax - Fmin;
        span(span<=0) = 1;
        for i = 1 : M
            [~,rank] = sortrows(PopObj(front,i));
            CrowdDis(front(rank(1)))   = inf;
            CrowdDis(front(rank(end))) = inf;
            for j = 2 : numel(front)-1
                CrowdDis(front(rank(j))) = CrowdDis(front(rank(j))) + ...
                    (PopObj(front(rank(j+1)),i) - PopObj(front(rank(j-1)),i)) / span(i);
            end
        end
    end
end
