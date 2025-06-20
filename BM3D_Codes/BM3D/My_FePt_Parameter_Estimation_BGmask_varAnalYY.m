function [Alpha, Sigma, curr_mean] = My_FePt_Parameter_Estimation_BGmask_varAnalYY(Data,AlpBinSize,AlpThFact,SigBinSize,Mask)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input
    %  Data: N*M*L array, N is 1st dim of image, M is 2nd dim of image, L
    %        is # of copies
    %  AlpBinSize: Bin size for determining alpha
    %  AlpThFact: a factor to determine what threshold will be used for
    %             determining alpha. For example, if AlpThFact is 1, all
    %             data will be used. If AlpThFact is 2, only high 50% of
    %             data will be used. If AlpThFact is 4, only high 25% of
    %             data will be used.
    %  SigBinSize: Bin size for determining sigma
    
    %% here is changed
    %  BGbound: boundary for background. if BGbound = [20 80], then the
    %           signal part will be assumed to be box bounded by four
    %           points (20,20), (20,80), (80,20), (80,80), and every point
    %           outside this box will be used for sigma determination.
    %  SigPrec: what precision will be used for sigma determination. if
    %           SigPrec is 1, the scan will stop if sigma is determined in
    %           accuracy of 1.
    %  SigNumHist: number of bins for intensity histogram used for sigma
    %              determination. 100 will be a good number.
    %  SigResampleFact: Poisson distribution will be resampled to make it
    %                  quasi-continuous, not integer descrete. The Poisson
    %                  distribution will be resampled with
    %                  2^(SigResampleFact) times finer than integer step.

    % Author: Yongsoo Yang
    %         Coherent Imaging Group, Dept. of Physics and Astronomy, UCLA
    
    % Date : 2015.10.9.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% alpha determination

    BinSize = AlpBinSize;

    Dim1 = size(Data,1);
    Dim2 = size(Data,2);
    NumCopy = size(Data,3);

    % binned data
    currSampleData = zeros(Dim1/BinSize,Dim2/BinSize,BinSize^2*NumCopy);

    Dim1Ind = 1:BinSize:Dim1;
    Dim2Ind = 1:BinSize:Dim2;

    % bin data
    currInd = 0;
    for j=1:BinSize	
        for k=1:BinSize
            tmpData = circshift(Data,-1*[j-1 k-1]);
            currSampleData(:,:,(currInd*NumCopy+1):((currInd+1)*NumCopy)) = tmpData(Dim1Ind,Dim2Ind,:);
            currInd = currInd + 1;
        end
    end


    currMeanData = squeeze(mean(currSampleData,3));

    AT = sort(currMeanData(:),'descend');

    Thresh = AT(round(length(AT)/AlpThFact));


    numPxl = sum(currMeanData(:) > Thresh);
    statsData = zeros(size(currSampleData,3),numPxl);
    curr_ind = 0;
    for i=1:Dim1/BinSize
        for j=1:Dim2/BinSize
            if currMeanData(i,j) > Thresh
                curr_ind = curr_ind + 1;
                statsData(:,curr_ind) = reshape(currSampleData(i,j,:),[size(currSampleData,3), 1]);
            end
        end
    end

    init_guess = 1;

    opts = optimset('Display','off');
    fun = @(p,xdata) sum((var(xdata / p(1),0,1)).*(mean(xdata / p(1),1)))/sum(mean(xdata / p(1),1).^2);
    %fun = @(p,xdata) sum(abs(var(xdata / p(1),0,1) - mean(xdata / p(1),1)));
    [p,~] = lsqcurvefit(fun,init_guess,statsData,1,[],[],opts);
    Alpha = p;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% sigma determination
    BinSize = SigBinSize;
    
    currSampleData = zeros(Dim1/BinSize,Dim2/BinSize,BinSize^2*NumCopy);

    Dim1Ind = 1:BinSize:Dim1;
    Dim2Ind = 1:BinSize:Dim2;

    currInd = 0;
    for j=1:BinSize	
        for k=1:BinSize
            tmpData = circshift(Data,-1*[j-1 k-1]);
            currSampleData(:,:,(currInd*NumCopy+1):((currInd+1)*NumCopy)) = tmpData(Dim1Ind,Dim2Ind,:);
            currInd = currInd + 1;
        end
    end

    
    binnedBigMask = zeros(size(Mask,1),size(Mask,2));
    
    for j=1:BinSize
        for k=1:BinSize
            binnedBigMask = binnedBigMask + circshift(Mask,-1*[j-1 k-1]);            
        end
    end
    binnedMask = binnedBigMask(1:BinSize:size(Mask,1),1:BinSize:size(Mask,2)) / BinSize^2;
    
    binnedTotMask = repmat(binnedMask,[1 1 size(currSampleData,3)]);
    
    %%%%%%%%%%%%%%This is where mask is applied
    bgMasked = currSampleData.*binnedTotMask;
    
    bgData = bgMasked(bgMasked~=0); %into linear array data
    
    curr_mean = mean(bgData(:));
    sample_var = var(bgData(:));
    meanalph = Alpha;
    if sample_var <= curr_mean*meanalph
      Sigma = 0;
    else
      Sigma = sqrt(sample_var - curr_mean*meanalph);
    end

end

