 classdef PaperTest2 < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Two-archive algorithm 2
% CAsize --- --- Convergence archive size
% p      --- --- The parameter of fractional distance

%------------------------------- Reference --------------------------------
% H. Wang, L. Jiao, and X. Yao, Two_Arch2: An improved two-archive
% algorithm for many-objective optimization, IEEE Transactions on
% Evolutionary Computation, 2015, 19(4): 524-541.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            [CAsize] = Algorithm.ParameterSet(Problem.N);

            %% Generate random population
            Population = Problem.Initialization();
            [CA,Fit] = UpdateCA([],Population,CAsize);
            DA = UpdateDA([],Population,Problem.N,W);


            %% Optimization
            while Algorithm.NotTerminated(DA)
               [ParentC,ParentM] = MatingSelection(CA,DA,Fit,Problem.FE,Problem.maxFE,Problem.N);
               Offspring         = [OperatorGA(Problem,ParentC,{1,15,0,0}),OperatorGA(Problem,ParentM,{0,0,1,15})];
               [CA,Fit] = UpdateCA(CA,Offspring,CAsize);
               DA = UpdateDA(DA,Offspring,Problem.N,W);
            end
        end
    end
end