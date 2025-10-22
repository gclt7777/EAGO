function Offspring = SE(Problem,SEPop,ymax,ymin)
Fitness = calFitness(SEPop.objs);
SEPopDec = SEPop.decs;
W = ((Fitness-ymin+1/inf)/(ymax-ymin+1/inf))*(Problem.maxFE-Problem.FE)/Problem.maxFE;
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
Offspring = Problem.Evaluation(OffDec,OffVel);
end