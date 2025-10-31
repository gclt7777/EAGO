function [DA,p] = UpdateDA(DA,New,MaxSize)
% Update DA

    DA = [DA,New];
    ND = NDSort(DA.objs,1);
    DA = DA(ND==1);
    N  = length(DA);
    [nInd, ~]    = size(DA.objs);
    CrowdDis = zeros(1,nInd);
    front1 = DA.objs;
    if (size(front1,1)>1)
        IdealPoint = min(front1);
    else
        IdealPoint = front1;
    end
    [CrowdDis, p, ~] = SurvivalScore(front1, IdealPoint);
    if N <= MaxSize
        return;
    end
    
    
    [~,Rank] = sort(CrowdDis,'descend');
    DA = DA(Rank(1:MaxSize));
end