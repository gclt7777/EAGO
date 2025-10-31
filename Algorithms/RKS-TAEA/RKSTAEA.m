classdef RKSTAEA < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            CAsize= Algorithm.ParameterSet(Problem.N);

            %% Generate random population
            Population = Problem.Initialization();
            [DA,p] = UpdateDA([],Population,Problem.N);
            CA = UpdateCA([],Population,CAsize,p);


            %% Optimization
            while Algorithm.NotTerminated(DA)
     
                [~,SDEP_fitness,~] = SDE_plus_indicator(CA.objs,1,p);         
                MatingPool = TournamentSelection(2,Problem.N,-SDEP_fitness);
                Offspring  = OperatorGA(Problem,CA(MatingPool)); 
                [DA,p] = UpdateDA(DA,Offspring,Problem.N);
                CA = UpdateCA(CA,Offspring,CAsize,p);

            end
        end
    end
end