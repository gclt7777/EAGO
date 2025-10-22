function Offspring = SE(Problem,SEPop,ymax,ymin)
Fitness = calFitness(SEPop.objs);
SEPopDec = SEPop.decs;
tiny = eps;
remainingRatio = max(0,1-get_progress_ratio(Problem));
W = ((Fitness-ymin+tiny)/(ymax-ymin+tiny)).*remainingRatio;
OffDec1 = SEPopDec+rand().*W.*(Problem.upper+Problem.lower-2*SEPopDec);
OffVel1 = OffDec1-SEPopDec;
OffDec2 = SEPopDec-rand().*W.*SEPopDec;
OffVel2 = OffDec2-SEPopDec;
OffDec = [OffDec1;OffDec1];
OffVel = [OffVel1;OffVel2];
%% Polynomial mutation
[N,D] = size(OffDec);
Lower  = repmat(Problem.lower,N,1);
Upper  = repmat(Problem.upper,N,1);
disM   = 20;
Site   = rand(N,D) < 1/D;
mu     = rand(N,D);
temp   = Site & mu<=0.5;
OffDec       = max(min(OffDec,Upper),Lower);
OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
    (1-(OffDec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
temp  = Site & mu>0.5;
OffDec(temp) = OffDec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
    (1-(Upper(temp)-OffDec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
Offspring = call_evaluation(Problem,OffDec,OffVel);
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

function Offspring = call_evaluation(Problem,OffDec,OffVel)
    if nargin < 3
        OffVel = [];
    end
    if has_evaluation_method(Problem)
        Offspring = Problem.Evaluation(OffDec,OffVel);
    else
        Offspring = INDIVIDUAL(OffDec,OffVel);
    end
end

function tf = has_evaluation_method(obj)
    tf = false;
    if isempty(obj)
        return;
    end
    try
        tf = ismethod(obj,'Evaluation');
        if tf
            return;
        end
    catch
        tf = false;
    end
    try
        m = methods(obj);
        tf = any(strcmp(m,'Evaluation'));
    catch
        tf = false;
    end
end
