function OffDec = Fuzzy_Search(Problem,Rate,Acc,Particles)
Total = 1;
S = floor(sqrt(2*Rate*Total/Acc));
Step = zeros(1,S+2);  % Step(1)=0，Step(S+2) is the compensation step
for i = 1 : S
    Step(i+1) = (S*i-i*i/2)*Acc;
end
Step(S+2) = Rate*Total;  % compensation step

%% Fuzzy Operation
ParticleDec = Particles.decs;
R    = Problem.upper-Problem.lower;
iter = get_progress_ratio(Problem);  % step=[0,0.6,0.8,0.8]
for i = 1 : S+1
    if iter>Step(i) && iter<=Step(i+1)
        gamma_a = R*10^-i.*floor(10^i*R.^-1.*(ParticleDec-Problem.lower)) + Problem.lower;
        gamma_b = R*10^-i.*ceil(10^i*R.^-1.*(ParticleDec-Problem.lower)) + Problem.lower;
                x1      = gamma_a(:,1:floor(Problem.D/2));
                y1      = gamma_a(:,(floor(Problem.D/2)+1):end);
                x2      = gamma_b(:,1:floor(Problem.D/2));
                y2      = gamma_b(:,(floor(Problem.D/2)+1):end);
                mu1     = (x1+3.*x2)/4;
                sigma1  = (x2-x1)/12;
                mu2     = (y1+3.*y2)/4;
                sigma2  = (y2-y1)/12;%+
                mu3     = (3.*x1+x2)/4;
                sigma3  = (x2-x1)/12;
                mu4     = (3.*y1+y2)/4;
                sigma4  = (y2-y1)/12;%-
                quadrant1_x = normrnd(mu1,sigma1);
                quadrant1_y = normrnd(mu2,sigma2);
                quadrant1   = [quadrant1_x,quadrant1_y];
                quadrant2_x = normrnd(mu3,sigma3);
                quadrant2_y = normrnd(mu2,sigma2);
                quadrant2   = [quadrant2_x,quadrant2_y];
                quadrant3_x = normrnd(mu3,sigma3);
                quadrant3_y = normrnd(mu4,sigma4);
                quadrant3   = [quadrant3_x,quadrant3_y];
                quadrant4_x = normrnd(mu1,sigma1);
                quadrant4_y = normrnd(mu4,sigma4);
                quadrant4   = [quadrant4_x,quadrant4_y];
    end
end
OffDec = [quadrant1;quadrant2;quadrant3;quadrant4];
end

function ratio = get_progress_ratio(Problem)
    [currentFE,maxFE] = get_evaluation_info(Problem);
    if maxFE <= 0
        ratio = 0;
    else
        ratio = max(0,min(1,currentFE./maxFE));
    end
end

function [currentFE,maxFE] = get_evaluation_info(Problem)
    currentFE = get_prop_or_field(Problem,{"FE","evaluated","Evaluated","currentFE","CurrentFE"});
    maxFE     = get_prop_or_field(Problem,{"maxFE","MaxFE","evaluation","Evaluation","maxEvaluations","MaxEvaluations","maxEvaluation","MaxEvaluation"});
    if isempty(currentFE)
        currentFE = 0;
    end
    if isempty(maxFE)
        maxFE = 0;
    end
end

function value = get_prop_or_field(Problem,candidates)
    value = [];
    for i = 1:numel(candidates)
        name = candidates{i};
        if isstruct(Problem) && isfield(Problem,name)
            candidate = Problem.(name);
        elseif isprop_safe(Problem,name)
            candidate = Problem.(name);
        else
            candidate = [];
        end
        if ~isempty(candidate)
            value = candidate;
            return;
        end
    end
end

function tf = isprop_safe(obj,name)
    try
        tf = isprop(obj,name);
    catch
        tf = false;
    end
end
