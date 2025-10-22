function MOEADDRA(Global)
% <algorithm> <M>
% MOEA/D-DRA (function-based wrapper for legacy PlatEMO)
%
% This implementation closely mirrors the original MOEA/D-DRA workflow
% while exposing a function entry point so that legacy PlatEMO releases,
% which expect algorithms to be defined as functions, can execute the
% algorithm without modification.
%
%------------------------------- Reference --------------------------------
% H. Li and Q. Zhang, Multiobjective optimization problems with complicated
% Pareto sets, MOEA/D and NSGA-II, IEEE Transactions on Evolutionary
% Computation, 2009, 13(2): 284-302.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference
% "Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB
% platform for evolutionary multi-objective optimization [educational
% forum], IEEE Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting (kept identical to the class-based implementation)
    [delta,T,nr,updateInterval,type] = Global.ParameterSet(0.9,20,2,30,1);

    %% Generate the weight vectors
    [W,Global.N] = UniformPoint(Global.N,Global.M);
    T = min(max(2,T),Global.N);

    %% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T);

    %% Generate random population
    Population = Global.Initialization();
    Z = min(Population.objs,[],1);

    %% Utility estimation state for dynamic resource allocation
    subprobUtility  = ones(Global.N,1);
    currentFitness  = calcScalarFitness(Population,W,Z,type);
    storedFitness   = currentFitness;
    generation      = 0;

    %% Optimization
    while Global.NotTermination(Population)
        generation = generation + 1;

        % Select the order of subproblems according to their utilities
        order = prioritySelection(subprobUtility);

        % For each selected subproblem
        for t = 1 : length(order)
            i = order(t);

            % Choose the parents
            if rand < delta
                P = B(i,randperm(size(B,2)));
            else
                P = randperm(Global.N);
            end

            % Generate an offspring by differential evolution
            Offspring = DE(Population(i),Population(P(1)),Population(P(2)));

            % Update the ideal point
            Z = min(Z,Offspring.obj);

            % Update the neighbouring solutions (or the whole population)
            if rand < delta
                neighbors = B(i,randperm(size(B,2)));
            else
                neighbors = randperm(Global.N);
            end
            [Population,currentFitness] = updateSubproblem(Population,Offspring,W,Z,type,neighbors,nr,currentFitness);
        end

        % Periodically update the resource allocation utilities
        if mod(generation,updateInterval) == 0
            [subprobUtility,storedFitness,currentFitness] = updateUtility(Population,W,Z,type,subprobUtility,storedFitness,currentFitness);
        end
    end
end

%% ---- Helper functions -------------------------------------------------
function fitness = calcScalarFitness(Population,W,Z,type)
    PopObj = Population.objs;
    N      = size(PopObj,1);
    fitness = zeros(N,1);
    for i = 1 : N
        fitness(i) = decomposition(PopObj(i,:),W(i,:),Z,type,PopObj);
    end
end

function [Population,fitness] = updateSubproblem(Population,Offspring,W,Z,type,indices,nr,fitness)
    PopObj = Population.objs;
    replaced = 0;
    for k = 1 : length(indices)
        idx = indices(k);
        g_old = decomposition(PopObj(idx,:),W(idx,:),Z,type,PopObj);
        g_new = decomposition(Offspring.obj,W(idx,:),Z,type,PopObj,Offspring.obj);
        if g_new <= g_old
            Population(idx) = Offspring;
            PopObj(idx,:)   = Offspring.obj;
            fitness(idx)    = g_new;
            replaced = replaced + 1;
        else
            fitness(idx) = g_old;
        end
        if replaced >= nr
            break;
        end
    end
end

function [utility,storedFitness,currentFitness] = updateUtility(Population,W,Z,type,utility,storedFitness,currentFitness)
    PopObj = Population.objs;
    refreshedFitness = zeros(length(Population),1);
    for i = 1 : length(Population)
        refreshedFitness(i) = decomposition(PopObj(i,:),W(i,:),Z,type,PopObj);
    end
    currentFitness = refreshedFitness;
    improvement = (storedFitness - currentFitness)./max(storedFitness,1e-12);
    improvement(~isfinite(improvement)) = 0;
    improvement(improvement < 0) = 0;
    avgImprove = mean(improvement);
    if avgImprove <= 0
        utility = ones(size(utility));
    else
        learningRate = 0.05;
        utility = (1-learningRate).*utility + learningRate.*(improvement/avgImprove);
    end
    utility = max(utility,1e-3);
    storedFitness = currentFitness;
end

function order = prioritySelection(utility)
    N = length(utility);
    order = zeros(1,N);
    available = true(1,N);
    weights   = utility(:)';
    if ~any(weights>0)
        order = randperm(N);
        return;
    end
    for i = 1 : N
        weights(~available) = 0;
        total = sum(weights);
        if total <= 0
            remaining = find(available);
            order(i:end) = remaining(randperm(length(remaining)));
            return;
        end
        r = rand()*total;
        cum = cumsum(weights);
        idx = find(cum >= r,1);
        if isempty(idx)
            idx = find(available,1);
        end
        order(i) = idx;
        available(idx) = false;
    end
end

function g = decomposition(obj,weight,Z,type,PopObj,candidate)
    if nargin < 6
        candidate = [];
    end
    switch type
        case 1
            % PBI approach
            diff = obj - Z;
            normW = sqrt(sum(weight.^2));
            normW = max(normW,1e-12);
            normDiff = sqrt(sum(diff.^2));
            if normDiff == 0
                g = 0;
                return;
            end
            Cosine = sum(diff.*weight)/(normW*normDiff);
            Cosine = max(min(Cosine,1),-1);
            g = normDiff*Cosine + 5*normDiff*sqrt(max(0,1-Cosine.^2));
        case 2
            % Tchebycheff approach
            g = max(abs(obj - Z).*weight);
        case 3
            % Normalized Tchebycheff approach
            Zmax = max(PopObj,[],1);
            if ~isempty(candidate)
                Zmax = max(Zmax,candidate);
            end
            scale = Zmax - Z;
            scale(scale <= 1e-12) = 1e-12;
            g = max(abs(obj - Z)./scale.*weight);
        case 4
            % Modified Tchebycheff approach
            weight(weight <= 1e-12) = 1e-12;
            g = max(abs(obj - Z)./weight);
        otherwise
            g = max(abs(obj - Z).*weight);
    end
end
