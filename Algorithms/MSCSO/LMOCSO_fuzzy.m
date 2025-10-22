function [Offspring,Winner,ymax,ymin] = LMOCSO_fuzzy(Problem,Particles,Rate,Acc)
Fitness = calFitness(Particles.objs);
ymax = max(Fitness);
ymin = min(Fitness);
if length(Particles) >= 2
    Rank = randperm(length(Particles),floor(length(Particles)/2)*2);
else
    Rank = [1,1];
end
Loser  = Rank(1:end/2);
Winner = Rank(end/2+1:end);
Change = Fitness(Loser) >= Fitness(Winner);
Temp   = Winner(Change);
Winner(Change) = Loser(Change);
Loser(Change)  = Temp;
WinnerDec  = Particles(Winner).decs;
LoserDec   = Particles(Loser).decs;
Winner_fuzzyDec = repmat(WinnerDec,4,1);
Loser_fuzzyDec  = Fuzzy_Search(Problem,Rate,Acc,Particles(Loser));
[N1,D1]      = size(LoserDec);
[N2,D2]      = size(Loser_fuzzyDec);
WinnerVel  = Particles(Winner).adds(zeros(N1,D1));
LoserVel  = Particles(Loser).adds(zeros(N1,D1));
Loser_fuzzyVel = repmat(LoserVel,4,1);
r1     = repmat(rand(N1,1),1,D1);
r2     = repmat(rand(N1,1),1,D1);
v     = r1.*LoserVel + r2.*(WinnerDec-LoserDec);
x     = LoserDec + v + r1.*(v-LoserVel);
OffDec = [WinnerDec;x];
OffVel = [WinnerVel;v];
r1_fuzzy     = repmat(rand(N2,1),1,D2);
r2_fuzzy     = repmat(rand(N2,1),1,D2);
v_fuzzy     = r1_fuzzy.*Loser_fuzzyVel + r2_fuzzy.*(Winner_fuzzyDec-Loser_fuzzyDec);
x_fuzzy     = Loser_fuzzyDec + v_fuzzy + r1_fuzzy.*(v_fuzzy-Loser_fuzzyVel);
OffDec = [OffDec;x_fuzzy];
OffVel = [OffVel;v_fuzzy];
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
