function NSGAII(Global)
% <algorithm> <M>
% NSGA-II (old-API wrapper)
% This function-based implementation mirrors the classic PlatEMO NSGA-II
% workflow so that the algorithm can be executed inside legacy releases
% that expect function handles instead of ALGORITHM subclasses.
%
%------------------------------- Reference --------------------------------
% K. Deb, A. Pratap, S. Agarwal, and T. Meyarivan, A fast and elitist
% multiobjective genetic algorithm: NSGA-II, IEEE Transactions on
% Evolutionary Computation, 2002, 6(2): 182-197.
%--------------------------------------------------------------------------

    %% Generate random population
    Population = Global.Initialization();
    [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Global.N);

    %% Optimization
    while Global.NotTermination(Population)
        % Binary tournament mating based on rank and crowding distance
        MatingPool = TournamentSelection(2,2*Global.N,FrontNo,-CrowdDis);
        Offspring  = GA(Population(MatingPool));

        % Environmental selection for the next generation
        [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N);
    end
end
