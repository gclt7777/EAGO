function [DDC] = Calculate_DDC(PopObj,N)
    zmax = max(PopObj,[],1);
    zmin = min(PopObj,[],1);
    DI = zeros(1,N);
    DN = zeros(1,N);
    DDC = zeros(1,N);
    for i = 1:N
        DI(i) = sqrt(sum((PopObj(i,:)-zmin).^2));
        DN(i) = sqrt(sum((zmax-PopObj(i,:)).^2));
    end
    maxC = min(DI);
    maxD = max(DN);
    for i = 1:N
        DDC(i) = max(DI(i)-maxC,maxD-DN(i));
    end
    DDI = DI-maxC;
    DDN = maxD-DN;
    [~,aindx] = sort(DDI,'ascend');
    [~,bindex] = sort(DDN,'descend');
    [~,cindex] = sort(DDC,'ascend');
end