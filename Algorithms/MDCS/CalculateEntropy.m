function [EntropyQ] = CalculateEntropy(UQpi,ConnectZNum,pi,N)
    
    EachZnum = zeros(1,ConnectZNum);
    %计算种群现在的信息熵
    for i =1:ConnectZNum
        a = find(pi == UQpi(i));
        EachZnum(i) = length(a);
    end
    EntropyNow = 0;
    for i = 1:length(EachZnum)
        EntropyNow = EntropyNow - (EachZnum(i)/N) *log2((EachZnum(i)/N));
    end

    %计算种群最佳信息熵
    bestZnum = floor(N/ConnectZNum);
    remainZnum = N - bestZnum*ConnectZNum;
    EntropyBest = -(ConnectZNum-remainZnum)*(bestZnum/N)*log2(bestZnum/N)-remainZnum*((bestZnum+1)/N)*log2(((bestZnum+1)/N));
    % bestZnumG = floor(N/NZ);
    % remainZnumG = N -bestZnumG*NZ;
    % EntropyBestG = -(NZ-remainZnumG)*(bestZnumG/N)*log2(bestZnumG/N)-remainZnumG*((bestZnumG+1)/N)*log2((bestZnumG+1)/N);
    % EntropyQG = EntropyNow /EntropyBestG;
    %利用信息熵比值计算领域数量
    EntropyQ = EntropyNow/EntropyBest;
end