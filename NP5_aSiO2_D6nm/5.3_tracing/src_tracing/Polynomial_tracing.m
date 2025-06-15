function [PosArr,TotPosArr,CoeffArr,Orders,exitFlagArr,maxXYZ] = ...
    Polynomial_tracing(FinalVol,MaxIter,CritIter,Th,SearchRad,Res,varargin)
if numel(varargin) > 0
    minDist_angs = varargin{1};
else
    minDist_angs = 1.5;
end
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

maxNumPeaks = 60000;
if length(maxPos) > maxNumPeaks
    maxPos = maxPos(1:60000);
end
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

minDist = minDist_angs / Res;

%% run fitting one by one
for i=1:size(maxXYZ,1)
    endFlag = 0;
    consecAccum = 0;
    iterNum = 0;
    while ~endFlag    
        iterNum = iterNum + 1;
        if iterNum>MaxIter
          exitFlagArr(i) = -4;
          endFlag = 1;
          break
        end
        cropXind = maxXYZ(i,1) + (-cropHalfSize:cropHalfSize);
        cropYind = maxXYZ(i,2) + (-cropHalfSize:cropHalfSize);
        cropZind = maxXYZ(i,3) + (-cropHalfSize:cropHalfSize);
        if min(cropXind) <= 0 || min(cropYind) <= 0 || min(cropZind) <= 0
            exitFlagArr(i) = -5;
            endFlag = 1;
            break
        end
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
end