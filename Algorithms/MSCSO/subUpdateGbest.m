function [Gbest,CrowdDis] = subUpdateGbest(Gbest,N)
% Update the global best set
Gbest    = Gbest(NDSort(Gbest.objs,1)==1);
CrowdDis = subCrowdingDistance(Gbest.objs);
[~,rank] = sort(CrowdDis,'descend');
Gbest    = Gbest(rank(1:min(N,length(Gbest))));
CrowdDis = CrowdDis(rank(1:min(N,length(Gbest))));
end