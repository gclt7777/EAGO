function [ParentC,ParentM] = MatingSelection(CA,DA,Fit,FE,maxFE,N)
% 匹配算子选择
    % CAParent1 = randi(length(CA),1,ceil(N/2));
    % CAParent2 = randi(length(CA),1,ceil(N/2));
    % Dominate  = any(CA(CAParent1).objs<CA(CAParent2).objs,2) - any(CA(CAParent1).objs>CA(CAParent2).objs,2);  
    % %从收敛性档案中选择N/2个非支配个体，从DA档案中选择N/2个个体
    % ParentC   = [CA([CAParent1(Dominate==1),CAParent2(Dominate~=1)]),...
    %              DA(randi(length(DA),1,ceil(N/2)))];
    % %随机打乱收敛性档案内的个体
    % ParentM   = CA(randi(length(CA),1,N));
    if FE/maxFE < 0.5 || (FE/maxFE > 0.5 && rand <0.5)
       MatingPool1 = TournamentSelection(2,ceil(N/2),Fit);
       ParentC = [CA(MatingPool1),DA(randi(length(DA),1,ceil(N/2)))];
       MatingPool2 = TournamentSelection(2,N,Fit);
       ParentM   = CA(MatingPool2);
    else
       DAParent1 = randi(length(DA),1,ceil(N/2));
       DAParent2 = randi(length(DA),1,ceil(N/2));
       Dominate  = any(DA(DAParent1).objs<DA(DAParent2).objs,2) - any(DA(DAParent1).objs>DA(DAParent2).objs,2);
       MatingPool1 = TournamentSelection(2,ceil(N/2),Fit);
       ParentC   = [DA([DAParent1(Dominate==1),DAParent2(Dominate~=1)]),...
                  CA(MatingPool1)];
       ParentM   = DA(randi(length(DA),1,N));
    end
end