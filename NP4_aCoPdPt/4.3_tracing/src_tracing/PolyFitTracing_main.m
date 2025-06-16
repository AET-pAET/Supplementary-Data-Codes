%% load data

% load reconstructed 3D volume
FinalVol = final_Rec;

% maximum number of iteration (the fit will be assumed to be failed if not
% converged within this number of iteration)
MaxIter = 15;

% critical iteration (this number of consecutive iterations should be
% successful, then it will be recorded as good peak)
CritIter = 5;

% threshold for minimum intensity to become a peak
Th = 0.00000001;

% radius for cropping the peak region for fitting
SearchRad = 3;

% pixel resolution
Res = 0.469;

ourputstring = 'Peak_Trace_NiPt_Nat_MS_test_interp2_Rad3_iter5';



%% initilize fit coefficients
fitCoeff = [];
for i=0:4
    for j=0:4
        for k=0:4
            if i+j+k <= 4                
                if max([i j k]) == 4
                    fitCoeff(end+1,:) = [i j k -1];                
                else                
                    fitCoeff(end+1,:) = [i j k 0];                
                end
            end
        end
    end
end

%%
% grayscale dialation to find maxima
se = strel3d(3);
dilatedBW = imdilate(FinalVol,se);
maxPos = find(FinalVol==dilatedBW & FinalVol>Th);
maxVals = FinalVol(maxPos);
[~,sortInd] = sort(maxVals,'descend');
maxPos = maxPos(sortInd);

fprintf(1,'numpeak = %d \n',length(maxPos));

maxXYZ = zeros(length(maxPos),3);
for i=1:length(maxPos)
    [xx,yy,zz] = ind2sub(size(FinalVol),maxPos(i));
    maxXYZ(i,:) = [xx yy zz];  
end

%%
Q = 0.5;
cropHalfSize = SearchRad;
[X,Y,Z] = ndgrid(-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize);
XYZdata.X = X;
XYZdata.Y = Y;
XYZdata.Z = Z;

Orders = fitCoeff(:,1:3);
PosArr = zeros(size(maxXYZ));
TotPosArr = zeros(size(maxXYZ));

% nonzero exit flag means the result is not good (fit failed, violation of
% minumum distance, etc).
exitFlagArr = zeros(1, size(maxXYZ,1));

CoeffArr = repmat(fitCoeff(:,4),[1 size(maxXYZ,1)]);

Alpha = 1;

minDist = 2 / Res;

% run fitting one by one
for i=1:size(maxXYZ,1)
    endFlag = 0;
    consecAccum = 0;
    iterNum = 0;
    while ~endFlag    
        iterNum = iterNum + 1;
        if iterNum>MaxIter
          exitFlagArr(i) = -4;
          endFlag = 1;
        end
        cropXind = maxXYZ(i,1) + (-cropHalfSize:cropHalfSize);
        cropYind = maxXYZ(i,2) + (-cropHalfSize:cropHalfSize);
        cropZind = maxXYZ(i,3) + (-cropHalfSize:cropHalfSize);

        cropVol = FinalVol(cropXind,cropYind,cropZind);

        Pos = PosArr(i,:);
        
        GaussWeight = exp(-1*Alpha*( (X-Pos(1)).^2 + (Y-Pos(2)).^2 + (Z-Pos(3)).^2 ) / cropHalfSize^2 );
        
        fun = @(p,xdata) calculate_3D_polynomial_Rogers(xdata.X,xdata.Y,xdata.Z,Pos,Orders,p).*GaussWeight;

        opts = optimset('Display','off');
        
        [p1,fminres1] = lsqcurvefit(fun,CoeffArr(:,i),XYZdata,cropVol.*GaussWeight,[],[],opts);
        CoeffArr(:,i) = p1;
        
        [dX, dY, dZ] = calc_dX_dY_dZ_Rogers(Orders,CoeffArr(:,i));
        if dX ==-100 && dY == -100 && dZ == -100
            exitFlagArr(i) = -1;
            endFlag = 1;
        else
            maxedShift = max([dX dY dZ],-1*[Q Q Q]);
            minedShift = min(maxedShift,[Q Q Q]);
            PosArr(i,:) = PosArr(i,:) + minedShift;
            if max(abs(PosArr(i,:))) > cropHalfSize
                exitFlagArr(i) = -2;
                endFlag = 1;
            elseif max(abs(minedShift)) < Q
                if consecAccum == CritIter-1
                    goodAtomTotPos = TotPosArr(1:i-1,:);
                    goodAtomTotPos = goodAtomTotPos(exitFlagArr(1:i-1)==0,:);
                    Dist = sqrt(sum((goodAtomTotPos - repmat(PosArr(i,:)+maxXYZ(i,:),[size(goodAtomTotPos,1) 1])).^2,2));
                    if min(Dist) < minDist
                      exitFlagArr(i) = -3;
                    else
                      TotPosArr(i,:) = PosArr(i,:) + maxXYZ(i,:);
                    end
                    endFlag = 1;
                else
                    consecAccum = consecAccum + 1;
                end
            else
                consecAccum = 0;
            end
        end
        
        

    end
    fprintf(1,'peak %d, flag %d \n',i,exitFlagArr(i));
end

%%
save(sprintf('%s_result.mat',ourputstring),'PosArr','TotPosArr','CoeffArr','Orders','exitFlagArr');


