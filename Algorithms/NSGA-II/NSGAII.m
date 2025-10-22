function NSGAII(Global)
% <2002> <multi> <real/integer> <constrained/none>
% NSGA-II for legacy PlatEMO (GLOBAL style), with robust operator calling.

    %% 初始化
    Population = Global.Initialization();
    [Population,FrontNo,CrowdDis] = NSGAII_EnvironmentalSelection(Population,Global.N);

    %% 主循环
    while Global.NotTermination(Population)
        % 竞赛选择（按前沿号与拥挤距离）
        MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);

        % 产生子代：优先调用平台外部算子；失败则内置 GA
        Parents   = Population(MatingPool);
        Offspring = nsga2_make_offspring(Global,Parents);

        % 环境选择
        [Population,FrontNo,CrowdDis] = NSGAII_EnvironmentalSelection([Population,Offspring],Global.N);
    end
end

% ==================== 环境选择 ====================
function [Population,FrontNo,CrowdDis] = NSGAII_EnvironmentalSelection(Population,N)
% 快速非支配排序 + 拥挤距离（兼容是否含约束）
    try
        cons = Population.cons;
    catch
        cons = [];
    end
    [FrontNo,MaxFNo] = NDSort(Population.objs,cons,N);
    Next = FrontNo < MaxFNo;

    CrowdDis = CrowdingDistance(Population.objs,FrontNo);

    Last = find(FrontNo==MaxFNo);
    if ~isempty(Last)
        [~,Rank] = sort(CrowdDis(Last),'descend');
        need = N - sum(Next);
        if need > 0
            Next(Last(Rank(1:need))) = true;
        end
    end
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end

% ==================== 子代生成（自适应外部/内置） ====================
function Offspring = nsga2_make_offspring(Global,Parents)
% 先试 OperatorGA / Operator 的不同参数顺序；都不行就用内置 SBX+PM
    % ---------- 1) 优先 OperatorGA ----------
    if exist('OperatorGA','file')==2 || exist('OperatorGA','builtin')==5
        % (Global,Parents)
        try, Offspring = OperatorGA(Global,Parents); return; end
        % (Parents,Global)
        try, Offspring = OperatorGA(Parents,Global); return; end
        % (Parents)
        try, Offspring = OperatorGA(Parents);        return; end
    end
    % ---------- 2) 其次 Operator ----------
    if exist('Operator','file')==2 || exist('Operator','builtin')==5
        % (Global,Parents)
        try, Offspring = Operator(Global,Parents);   return; end
        % (Parents,Global)
        try, Offspring = Operator(Parents,Global);   return; end
        % (Parents)
        try, Offspring = Operator(Parents);          return; end
    end

    % ---------- 3) 兜底：内置 SBX + 多项式变异 ----------
    Decs = Parents.decs;
    Np   = size(Decs,1);
    D    = size(Decs,2);

    % 成对，奇数补齐
    if mod(Np,2)==1
        Decs = [Decs; Decs(end,:)];
        Np   = Np + 1;
    end

    lb = Global.lower(:)';  ub = Global.upper(:)';
    % 简单识别整数维度（边界均为整数）
    isInt = (abs(lb - round(lb))<eps) & (abs(ub - round(ub))<eps);

    % 参数
    proC = 1.0;  disC = 20;     % SBX
    proM = 1.0;  disM = 20;     % 多项式变异
    pmVar = proM / D;

    OffDec = zeros(Np,D);
    for i = 1:2:Np
        p1 = Decs(i,:);  p2 = Decs(i+1,:);
        c1 = p1;         c2 = p2;

        % SBX 交叉
        if rand < proC
            u    = rand(1,D);
            beta = zeros(1,D);
            mask = rand(1,D) <= 0.5;
            beta(mask)  = (2.*u(mask)).^(1/(disC+1));
            beta(~mask) = (2.*(1-u(~mask))).^(-1/(disC+1));
            c1 = 0.5*((1+beta).*p1 + (1-beta).*p2);
            c2 = 0.5*((1-beta).*p1 + (1+beta).*p2);
        end

        % 多项式变异
        c1 = poly_mutation(c1,lb,ub,pmVar,disM);
        c2 = poly_mutation(c2,lb,ub,pmVar,disM);

        % 截断 + 整数取整
        c1 = min(max(c1,lb),ub);
        c2 = min(max(c2,lb),ub);
        if any(isInt)
            c1(isInt) = round(c1(isInt));
            c2(isInt) = round(c2(isInt));
        end

        OffDec(i,:)   = c1;
        OffDec(i+1,:) = c2;
    end

    if length(Parents) < Np
        OffDec(end,:) = [];
    end
    Offspring = INDIVIDUAL(OffDec);
end

% ==================== 多项式变异（单个个体一行） ====================
function x = poly_mutation(x,lb,ub,pmVar,disM)
    D  = numel(x);
    do = rand(1,D) < pmVar;
    if ~any(do), return; end

    xl = lb(do); xu = ub(do);
    xi = x(do);

    span = max(xu - xl, eps);
    delta1 = (xi - xl) ./ span;
    delta2 = (xu - xi) ./ span;

    rnd = rand(1,sum(do));
    mut_pow = 1/(disM+1);
    mask = rnd <= 0.5;

    % 左侧
    xy = 1 - delta1(mask);
    val = 2.*rnd(mask) + (1-2.*rnd(mask)).*(xy.^(disM+1));
    deltaq = val.^mut_pow - 1;
    xi(mask) = xi(mask) + deltaq .* span(mask);

    % 右侧
    xy = 1 - delta2(~mask);
    val = 2.*(1 - rnd(~mask)) + 2.*(rnd(~mask)-0.5).*(xy.^(disM+1));
    deltaq = 1 - val.^mut_pow;
    xi(~mask) = xi(~mask) + deltaq .* span(~mask);

    x(do) = min(max(xi, xl), xu);
end
