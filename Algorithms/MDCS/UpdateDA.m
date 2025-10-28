function DA = UpdateDA(DA,New,MaxSize,W)
% Update DA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Find the non-dominated solutions
    DA = [DA,New];
    ND = NDSort(DA.objs,1); 
    DA = DA(ND==1);  %保留非支配解
    N  = length(DA);
    if N <= MaxSize
        return;
    end
    Popobj = DA.objs; %第一个非支配层级保留的解
    
    %% 归一化
    Zmin   = min(Popobj,[],1);
    Zmax   = max(Popobj,[],1);
    Popobj = (Popobj-repmat(Zmin,size(Popobj,1),1))./repmat(Zmax-Zmin,size(Popobj,1),1);
    Choose = false(1,N);

    %% 关联参考向量
    NZ = size(W,1);
    Cosine = 1 - pdist2(Popobj,W,'cosine');
    Distance = repmat(sqrt(sum(Popobj.^2,2)),1,NZ).*sqrt(1-Cosine.^2);
    %每个解对应的参考向量pi及其与参考向量之间的距离d
    [d,pi] = min(Distance',[],1);
    %获取每个参考向量关联的个体
    [~,index] = min(Distance,[],1);
    %设置每个参考向量上的密度为0
    %rho = zeros(1,NZ);
    

    %% 选择极值点
    M = size(Popobj,2);    
    Extreme = zeros(1,M);
    w       = zeros(M)+1e-6+eye(M);
    for i = 1 : M
        objmatrix = DA(~Choose).objs;
        [~,Extreme(i)] = min(max(DA(~Choose).objs./repmat(w(i,:),sum(~Choose),1),[],2)+0.1*objmatrix(:,i)/(1e-6)); % 带惩罚的AFS函数值
        Choose(Extreme(i)) = true; 
        %rho(pi(Extreme(i))) = rho(pi(Extreme(i)))+1;
    end

    %遍历每个参考向量
    for i =1:NZ
        Choose(index(i)) = true;
    end

    if sum(Choose) > MaxSize
        %随机删除一些解
        Choosed = find(Choose);
        k = randperm(sum(Choose),sum(Choose)-MaxSize);
        Choose(Choosed(k)) = false;
    elseif sum(Choose) < MaxSize
        %% 计算种群信息熵
        UQpi = unique(pi);
        ConnectZNum = length(UQpi);
        EntropyQ = CalculateEntropy(UQpi,ConnectZNum,pi,N);
        % a = floor(N / ConnectZNum);
        q = floor(N /ConnectZNum* (1-EntropyQ));
        if q < 1
            q =1;
        end
        k = MaxSize- sum(Choose);
        divValue = CalculateDiv_Test(Popobj,q);
        Unselected = find(Choose == 0);
        UnselectedDiv = divValue(Unselected);
        [~,indexDiv] = sort(UnselectedDiv,'descend');
        needSelect = indexDiv(1:k);
        Choose(Unselected(needSelect)) = true;
    end


    




    % %删除或者增加解
    % if sum(Choose) > MaxSize
    %     %随机删除一些解
    %     Choosed = find(Choose);
    %     k = randperm(sum(Choose),sum(Choose)-MaxSize);
    %     Choose(Choosed(k)) = false;
    % elseif sum(Choose) < MaxSize
    %     %参考向量信息
    %     Zchoose = false(1,NZ);
    %     Zchoose(UQpi) = true;
    %     %计算解的多样性
    %     divValue = CalculateDiv_Test(Popobj,q); 
    %     div = d ./(divValue+1);
    %     %进行环境选择
    %     while sum(Choose) < MaxSize
    %         %选择一个最不拥挤的参考向量
    %         Temp = find(Zchoose);
    %         Jmin = find(rho(Temp)==min(rho(Temp)));
    %         j = Temp(Jmin(randi(length(Jmin))));
    %         I = find(pi==j & Choose==0);
    %         if ~isempty(I)
    %             [~,s] = min(div(I));
    %             Choose(I(s)) = true;
    %             rho(j) = rho(j)+1;
    %         else
    %             Zchoose(j) = false;
    %         end
    %     end
    % end   
    DA = DA(Choose);
end