function DVCOEA(Global)
%DVCOEA Multi-population MOEA based on decision variable contribution.
%
%   This implementation follows the idea of the decision-variable
%   contribution driven co-evolutionary algorithm (DVCOEA). Decision
%   variables are first analysed to determine their convergence and
%   diversity contributions to each objective. Convergence-related
%   variables are further assigned to the objectives that benefit most from
%   their perturbations, which yields multiple sub-populations that operate
%   in dedicated subspaces. Each sub-population performs differential
%   evolution restricted to its active variables, while an external archive
%   maintains elite solutions and supports redistribution of individuals.
%
% <algorithm> <L>
%
%------------------------------- Reference --------------------------------
% K. Li, Q. Zhang, G. G. Yen, and S. Kwong, A Novel Multi-Population Based
% Evolutionary Algorithm for Large Scale Multiobjective Optimization,
% IEEE Transactions on Evolutionary Computation, 2016, 20(5): 757-772.
%
% J. Zou, Y. Jin, Q. Zhang, H. Wang, and J. Zhang, A Decision Variable
% Contribution Based Evolutionary Algorithm for Large-Scale Many-Objective
% Optimization, IEEE Transactions on Cybernetics, 2021, 51(1): 27-40.
%--------------------------------------------------------------------------

    %% Parameter setting
    [nSample,nPerturb,F,CR,archiveFactor,reclassGap,subPopSize] = ...
        Global.ParameterSet(3,2,0.5,0.9,1.5,20,0);

    %% Initial population and archive
    Population = Global.Initialization();
    archiveSize = max(Global.N, round(archiveFactor*Global.N));
    Archive     = UpdateArchive(Population,INDIVIDUAL.empty,archiveSize);

    %% Decision variable contribution analysis
    analysedSet = Population;
    [CV,DV,NV,objectiveImpact,dominantObj,impactScore] = ...
        VariableContributionAnalysis(Global,analysedSet,nSample,nPerturb);

    subsets    = BuildObjectiveSubspaces(CV,DV,dominantObj,impactScore,Global.M,Global.D);
    subPopSize = DetermineSubpopulationSize(subPopSize,Global.N,Global.M);
    SubPops    = InitializeSubpopulations(Population,Archive,subPopSize,Global);

    generation = 0;

    %% Optimization loop
    while Global.NotTermination(Population)
        generation = generation + 1;

        Offspring = INDIVIDUAL.empty;
        for m = 1 : Global.M
            [offs,SubPops{m}] = EvolveSubpopulation(SubPops{m},subsets{m},Archive,F,CR,m,subPopSize,Global);
            Offspring = [Offspring,offs];
        end

        Union       = [Population,Offspring];
        Population  = DVCOEASelection([Union,Archive],Global.N);
        Archive     = UpdateArchive([Population,Offspring],Archive,archiveSize);
        SubPops     = RedistributeSubpopulations(Population,Archive,subPopSize,Global);

        if reclassGap > 0 && mod(generation,reclassGap) == 0
            analysedSet = [Population,Archive];
            [CV,DV,NV,objectiveImpact,dominantObj,impactScore] = ...
                VariableContributionAnalysis(Global,analysedSet,nSample,nPerturb);
            subsets = BuildObjectiveSubspaces(CV,DV,dominantObj,impactScore,Global.M,Global.D);
            SubPops = RedistributeSubpopulations(Population,Archive,subPopSize,Global);
        end
    end
end

function [CV,DV,NV,objImpact,dominantObj,impactScore] = ...
        VariableContributionAnalysis(Global,Population,nSample,nPerturb)
% Analyse the contribution of each decision variable using perturbations.

    D = Global.D;
    M = Global.M;
    if isempty(Population)
        CV = [];
        DV = [];
        NV = 1:D;
        objImpact   = zeros(D,M);
        dominantObj = zeros(D,1);
        impactScore = zeros(D,1);
        return;
    end

    nSample = max(1,min(nSample,length(Population)));
    sampleIndex = randperm(length(Population),nSample);
    Samples     = Population(sampleIndex);

    Lower = Global.lower;
    Upper = Global.upper;
    if isempty(Lower)
        Lower = zeros(1,D);
    end
    if isempty(Upper)
        Upper = ones(1,D);
    end
    step = max(1e-3,(Upper-Lower))*0.05;

    objImpact   = zeros(D,M);
    impactScore = zeros(D,1);
    aggDelta    = zeros(D,1);
    aggDeltaSq  = zeros(D,1);
    hitCount    = zeros(D,1);

    for s = 1 : nSample
        baseDec = Samples(s).decs;
        baseObj = Samples(s).objs;
        for i = 1 : D
            for p = 1 : max(1,nPerturb)
                direction = 2*rand()-1;
                perturbedDec       = baseDec;
                perturbedDec(i)    = baseDec(i) + direction*step(i);
                perturbedDec(i)    = min(max(perturbedDec(i),Lower(i)),Upper(i));
                if perturbedDec(i) == baseDec(i)
                    continue;
                end
                perturbed = INDIVIDUAL(perturbedDec);
                delta     = perturbed.objs - baseObj;
                objImpact(i,:)   = objImpact(i,:) + abs(delta);
                change           = sum(delta);
                impactScore(i)   = impactScore(i) + sum(abs(delta));
                aggDelta(i)      = aggDelta(i) + change;
                aggDeltaSq(i)    = aggDeltaSq(i) + change.^2;
                hitCount(i)      = hitCount(i) + 1;
            end
        end
    end

    valid = hitCount > 0;
    avgImpact = zeros(D,1);
    avgImpact(valid) = impactScore(valid)./hitCount(valid);
    impactScore = avgImpact;

    avgAgg = zeros(D,1);
    avgAgg(valid) = aggDelta(valid)./hitCount(valid);
    stdAgg = zeros(D,1);
    stdAgg(valid) = sqrt(max(0,aggDeltaSq(valid)./hitCount(valid) - avgAgg(valid).^2));

    if any(valid)
        baseThreshold = median(avgImpact(valid));
    else
        baseThreshold = 0;
    end
    threshold = max(1e-6,baseThreshold * 0.1);
    NV = find(~valid | avgImpact <= threshold);

    convScore = zeros(D,1);
    convScore(valid) = max(0,-avgAgg(valid));
    divScore  = zeros(D,1);
    divScore(valid)  = stdAgg(valid);

    CV = find(avgImpact > threshold & convScore >= divScore);
    DV = find(avgImpact > threshold & convScore < divScore);

    objImpact(~valid,:) = 0;
    dominantObj = zeros(D,1);
    for i = reshape(CV,1,[])
        [~,idx] = max(objImpact(i,:));
        dominantObj(i) = idx;
    end
end

function subsets = BuildObjectiveSubspaces(CV,DV,dominantObj,impactScore,M,D)
% Build subspaces for each objective based on dominant contributions.

    subsets = cell(1,M);
    for m = 1 : M
        subsets{m} = [];
    end

    for i = reshape(CV,1,[])
        obj = dominantObj(i);
        if obj < 1 || obj > M
            continue;
        end
        subsets{obj}(end+1) = i; %#ok<AGROW>
    end

    % Assign high-impact diversity variables if some objectives lack support
    remainingDV = DV(:)';
    [~,order] = sort(impactScore(remainingDV),'descend');
    remainingDV = remainingDV(order);

    for m = 1 : M
        subsets{m} = unique(subsets{m});
        if isempty(subsets{m}) && ~isempty(remainingDV)
            take = min(max(1,ceil(numel(remainingDV)/(M-m+1))),numel(remainingDV));
            subsets{m} = remainingDV(1:take);
            remainingDV(1:take) = [];
        end
    end

    % Fallback: ensure every subset contains at least one index
    for m = 1 : M
        if isempty(subsets{m})
            subsets{m} = 1:D;
        end
    end
end

function subPopSize = DetermineSubpopulationSize(requested,N,M)
% Determine the working size of each subpopulation.

    if requested <= 0
        subPopSize = max(5,ceil(N/M));
    else
        subPopSize = requested;
    end
end

function SubPops = InitializeSubpopulations(Population,Archive,subPopSize,Global)
% Build the initial subpopulations from the current population and archive.

    M = Global.M;
    pool = [Population,Archive];
    if isempty(pool)
        SubPops = repmat({INDIVIDUAL.empty},1,M);
        return;
    end

    SubPops = cell(1,M);
    for m = 1 : M
        SubPops{m} = SelectByObjective(pool,m,subPopSize);
    end
end

function SubPops = RedistributeSubpopulations(Population,Archive,subPopSize,Global)
% Redistribute individuals to subpopulations after population update.

    M = Global.M;
    pool = [Population,Archive];
    SubPops = cell(1,M);
    for m = 1 : M
        if isempty(pool)
            SubPops{m} = INDIVIDUAL.empty;
        else
            SubPops{m} = SelectByObjective(pool,m,subPopSize);
        end
    end
end

function [Offspring,SubPop] = EvolveSubpopulation(SubPop,subset,Archive,F,CR,objIndex,subPopSize,Global)
% Perform subspace-restricted DE within a subpopulation.

    if isempty(SubPop)
        Offspring = INDIVIDUAL.empty;
        return;
    end
    if isempty(subset)
        subset = 1:Global.D;
    end

    baseDecs = SubPop.decs;
    N = size(baseDecs,1);

    pool = [SubPop,Archive];
    poolDecs = pool.decs;
    if size(poolDecs,1) < 3
        deficit = 3 - size(poolDecs,1);
        replicate = SubPop(randi(length(SubPop),1,deficit));
        pool      = [pool,replicate];
        poolDecs  = pool.decs;
    end

    idx2 = randi(size(poolDecs,1),N,1);
    idx3 = randi(size(poolDecs,1),N,1);
    same = idx2 == idx3;
    if any(same)
        idx3(same) = mod(idx3(same)+1,size(poolDecs,1));
        idx3(idx3 == 0) = size(poolDecs,1);
    end

    donorSub = baseDecs(:,subset) + F*(poolDecs(idx2,subset) - poolDecs(idx3,subset));
    mask     = rand(N,length(subset)) < CR;
    for i = 1 : N
        if ~any(mask(i,:))
            mask(i,randi(length(subset))) = true;
        end
    end
    trial = baseDecs;
    baseSubset = baseDecs(:,subset);
    baseSubset(mask) = donorSub(mask);
    trial(:,subset) = baseSubset;

    Lower = repmat(Global.lower,N,1);
    Upper = repmat(Global.upper,N,1);
    if ~isempty(Lower)
        trial = max(trial,Lower);
    end
    if ~isempty(Upper)
        trial = min(trial,Upper);
    end

    Offspring = INDIVIDUAL(trial);

    Combined = [SubPop,Offspring];
    SubPop   = SelectByObjective(Combined,objIndex,subPopSize);
end

function Selected = SelectByObjective(Pop,objIndex,maxSize)
% Select individuals with respect to a particular objective and feasibility.

    if isempty(Pop)
        Selected = INDIVIDUAL.empty;
        return;
    end

    cons = Pop.cons;
    if isempty(cons)
        violation = zeros(length(Pop),1);
    else
        violation = sum(max(0,cons),2);
    end

    objs = Pop.objs(:,objIndex);
    [~,order] = sortrows([violation,objs]);
    pick = min(maxSize,length(order));
    Selected = Pop(order(1:pick));
end

function Population = DVCOEASelection(Population,N)
% NSGA-II style environmental selection with crowding preservation.

    if isempty(Population)
        return;
    end

    PopObj  = Population.objs;
    PopCons = Population.cons;
    if isempty(PopCons)
        PopCons = [];
    end

    [FrontNo,MaxFront] = NDSort(PopObj,PopCons,N);
    Next = FrontNo < MaxFront;

    CrowdDis = CalculateCrowdingDistance(PopObj,FrontNo);

    Last = find(FrontNo == MaxFront);
    Remain = N - sum(Next);
    if Remain > 0 && ~isempty(Last)
        [~,rank] = sort(CrowdDis(Last),'descend');
        picks    = min(Remain,length(rank));
        Next(Last(rank(1:picks))) = true;
    end

    Population = Population(Next);
end

function CrowdDis = CalculateCrowdingDistance(PopObj,FrontNo)
% Calculate crowding distances front by front.

    [N,M]    = size(PopObj);
    CrowdDis = zeros(1,N);
    fronts   = setdiff(unique(FrontNo),inf);
    for f = 1 : numel(fronts)
        front = find(FrontNo == fronts(f));
        if numel(front) < 2
            CrowdDis(front) = inf;
            continue;
        end
        Fmax = max(PopObj(front,:),[],1);
        Fmin = min(PopObj(front,:),[],1);
        span = Fmax - Fmin;
        for i = 1 : M
            [~,rank] = sort(PopObj(front,i));
            CrowdDis(front(rank(1)))   = inf;
            CrowdDis(front(rank(end))) = inf;
            if span(i) <= 0
                continue;
            end
            for j = 2 : length(front)-1
                CrowdDis(front(rank(j))) = CrowdDis(front(rank(j))) + ...
                    (PopObj(front(rank(j+1)),i) - PopObj(front(rank(j-1)),i)) / span(i);
            end
        end
    end
end

function Archive = UpdateArchive(Population,Archive,maxArchive)
% Update the external archive via non-dominated sorting and crowding.

    Combined = [Archive,Population];
    if isempty(Combined)
        Archive = Combined;
        return;
    end

    Decs = Combined.decs;
    rounded = round(Decs,8);
    [~,uniqueIdx] = unique(rounded,'rows','stable');
    Combined = Combined(uniqueIdx);

    PopObj  = Combined.objs;
    PopCons = Combined.cons;
    if isempty(PopCons)
        PopCons = [];
    end

    [FrontNo,MaxFront] = NDSort(PopObj,PopCons,maxArchive);
    Next    = FrontNo < MaxFront;
    Archive = Combined(Next);

    Last  = find(FrontNo == MaxFront);
    Remain = maxArchive - sum(Next);
    if Remain > 0 && ~isempty(Last)
        CrowdDis = CalculateCrowdingDistance(PopObj,FrontNo);
        [~,rank] = sort(CrowdDis(Last),'descend');
        picks    = min(Remain,length(rank));
        Archive  = [Archive,Combined(Last(rank(1:picks)))];
    end
end
