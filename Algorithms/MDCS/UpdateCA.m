function [CA,Fit] = UpdateCA(CA,New,MaxSize)
% Update CA
    CA = [CA,New];
    %ND = NDSort(CA.objs,1);
    %CA = CA(ND == 1);
    N  = length(CA);
    if N < MaxSize
        return;
    end
    
    % 只选择非支配解，并且在非支配解里面选择排名前N的个体
    CAObj = CA.objs;
    Choose = false(1,N);
    ddC = Calculate_DDC(CAObj,N);
    [~,index] = sort(ddC,'ascend');
    Choose(index(1:MaxSize)) = true;
    CA = CA(Choose);
    Fit = ddC(Choose);
end