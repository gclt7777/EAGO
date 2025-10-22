function MSCSO(Global)
% <algorithm> <L>
% MSCSO (old-API wrapper)
% 说明：由 classdef(ALGORITHM) 版本改写为老版 PlatEMO 的 function 版本。
%
% 依赖：UniformPoint, EnvironmentalSelection, LMOCSO_fuzzy, SE,
%       calFitness, Operator（这些保持原函数名不变）

    %% 参数读取（与原版一致）
    % Rate, Acc 对应原 Algorithm.ParameterSet(0.2,2/5)
    [Rate,Acc] = Global.ParameterSet(0.2,2/5);

    %% 参考向量与种群初始化
    [V,Global.N] = UniformPoint(Global.N,Global.M);
    Particles    = Global.Initialization();

    % 进度（老版用 Global.evaluated / Global.evaluation）
    prog = progress_ratio(Global);
    Particles = EnvironmentalSelection(Particles, V, prog^2);

    %% 主循环
    while Global.NotTermination(Particles)
        phi = 0.7;
        prog = progress_ratio(Global);

        if prog < Rate
            % —— 模糊分治阶段 ——
            [Offspring,Winner,ymax,ymin] = LMOCSO_fuzzy(Global,Particles,Rate,Acc);
            Offspring1 = SE(Global,Particles(Winner),ymax,ymin);
            Particles  = EnvironmentalSelection([Particles,Offspring,Offspring1], V, prog^2);
        else
            % —— 竞争-学习阶段 ——
            Fitness = calFitness(Particles.objs);

            if length(Particles) >= 2
                % 保证偶数个用于两两配对
                Rank  = randperm(length(Particles), floor(length(Particles)/2)*2);
            else
                Rank = [1,1];
            end

            Loser  = Rank(1:end/2);
            Winner = Rank(end/2+1:end);

            % 让 Winner 的适应度更优（若相反则交换）
            Change = Fitness(Loser) >= Fitness(Winner);
            Temp   = Winner(Change);
            Winner(Change) = Loser(Change);
            Loser(Change)  = Temp;

            % 产生子代
            Offspring = Operator(Global, Particles, Particles(Loser), Particles(Winner), phi);

            % 环境选择
            Particles = EnvironmentalSelection([Particles,Offspring], V, prog^2);
        end
    end
end

%% ---- 辅助：计算进度比（老版 PlatEMO 字段名适配） ----
function p = progress_ratio(Global)
    % 防越界：初始化后 Global.evaluated >= Global.N
    if isprop(Global,'evaluated')
        fe = Global.evaluated;
    else
        % 极老版本兜底（若字段名不同，可在此改）
        fe = Global.FE;
    end
    if isprop(Global,'evaluation')
        feMax = Global.evaluation;
    else
        feMax = Global.maxFE;
    end
    if feMax <= 0
        p = 0;
    else
        p = max(0,min(1, fe/feMax));
    end
end
