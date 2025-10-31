function CA = UpdateCA(CA,New,MaxSize,p)
% Update CA


    CA = [CA,New];
    N  = length(CA);
    if N <= MaxSize
        return;
    end
    
    %% Calculate the fitness of each solution
       [SDE_fitness, ~, SDE_Distance] = SDE_plus_indicator(CA.objs, 0, p);
    
    %% Delete part of the solutions by their fitnesses
     if(sum(SDE_fitness ~= 0) <= MaxSize)
        %  One-time evaluation
        [~,Remain] = sortrows(sort(SDE_Distance, 2),'descend');
        Remain = Remain(1:MaxSize);
        else
        % Dynamic evaluation
        Choose = SDE_fitness ~= 0;
        while sum(Choose) ~= MaxSize
            Remain = find(Choose);
            mat = min(SDE_Distance(Choose,Choose),[],2);
            [~,x]  = min(mat);
%             mat = min(SDE_Distance(Choose,Choose),2);
%             [~,x]  = min(min(mat,[],2));
            Choose(Remain(x)) = false;
        end
    end
        CA = CA(Remain);


end