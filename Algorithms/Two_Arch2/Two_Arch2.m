function Two_Arch2(Global)
% <algorithm> <M>
% Two_Archive2 (old-API wrapper)
% This function-based variant mirrors the legacy PlatEMO workflow so that
% Two_Archive2 can run inside historical releases expecting function
% handles rather than ALGORITHM subclasses.
%
%------------------------------- Reference --------------------------------
% M. Li, S. Yang, and X. Liu, An improved two-archive algorithm for
% many-objective optimization, European Journal of Operational Research,
% 2014, 243(2): 556-567.
%--------------------------------------------------------------------------

    %% Generate initial population and the two archives
    Population = Global.Initialization();
    [CA,DA,FrontNo,CrowdDis] = EnvironmentalSelection(Population, Global.N);

    %% Optimization
    while Global.NotTermination(CA)
        %% Exploitative mating from the convergence archive (CA)
        numGA = 2 * ceil(Global.N/2);
        if isempty(CA)
            OffspringCA = INDIVIDUAL.empty();
        else
            % Tournament selection prefers lower ranks and larger crowding
            if numel(CA) == 1
                parents = repmat(1,1,numGA);
            else
                parents = TournamentSelection(2, numGA, FrontNo, -CrowdDis);
            end
            OffspringCA = GAhalf(CA(parents));
        end

        %% Explorative mating from the diversity archive (DA)
        numDE = max(Global.N - numel(OffspringCA), 0);
        if numDE > 0 && ~isempty(DA) && ~isempty(CA)
            baseIdx   = randi(numel(DA), 1, numDE);
            guideIdx1 = randi(numel(CA), 1, numDE);
            guideIdx2 = randi(numel(CA), 1, numDE);
            OffspringDA = DE(DA(baseIdx), CA(guideIdx1), CA(guideIdx2));
        else
            OffspringDA = INDIVIDUAL.empty();
        end

        %% Update both archives using the newly generated offspring
        Merged = [CA, DA, OffspringCA, OffspringDA];
        [CA,DA,FrontNo,CrowdDis] = EnvironmentalSelection(Merged, Global.N);
    end
end
