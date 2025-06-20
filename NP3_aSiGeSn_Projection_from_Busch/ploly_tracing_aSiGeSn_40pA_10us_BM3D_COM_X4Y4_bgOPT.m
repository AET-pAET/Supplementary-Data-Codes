% Polynomial tracing based on reconstruction of original orientation
SearchRadAr = [3];
prefix = 'gd_cuda_SiGeSn_test';
Res = 0.195;

addpath('./src');
addpath('./src_tracing');

homeDir = './';
inputDir = './';


Recon_filename = [homeDir 'Resire_SiGeSn_40pA_10us_BM3D_COM_X4Y4_bgOPT.mat']
ourputstring_prefix = [inputDir 'polyn_tracing_RESIRE_SiGeSn_40pA']


classify_info = struct('Num_species', 3,  'halfSize',  3,  'plothalfSize',  1, ...
    'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  0,  'separate_part',  120);

%% load reconstructed 3D volume
disp(Recon_filename);
Dsetvol = importdata(Recon_filename);
Dsetvol = Dsetvol.reconstruction;
load('SiGeSn_MASK_VST_SIRT.mat')
Dsetvol = double(Dsetvol).*MASK*1500;

xx = (1:size(Dsetvol,1)) - round((size(Dsetvol,1)+1)/2);
yy = (1:size(Dsetvol,2)) - round((size(Dsetvol,2)+1)/2);
zz = (1:size(Dsetvol,3)) - round((size(Dsetvol,3)+1)/2);
xxi = ((3*xx(1)):(xx(end)*3))/3;
yyi = ((3*yy(1)):(yy(end)*3))/3;
zzi = ((3*zz(1)):(zz(end)*3))/3;
xxi = xxi(3:end); yyi = yyi(3:end); zzi = zzi(3:end);
[Y,X,Z] = meshgrid(yy,xx,zz);
[Yi,Xi,Zi] = meshgrid(yyi,xxi,zzi);
Dsetvol = interp3(Y,X,Z,Dsetvol,Yi,Xi,Zi,'spline',0);
clear Xi Yi Zi

FinalVol = My_paddzero(Dsetvol,size(Dsetvol)+20);
FinalVol_single = single(FinalVol);
save('Rec_SiGeSn_40pA_10us_BM3D_COM_X4Y4_bgOPT.mat','FinalVol_single','-v7.3');
% FinalVol_single = importdata(intpRec_filename);
FinalVol = double(FinalVol_single);

[ii,jj,kk] = meshgrid(0:4,0:4,0:4);
fitCoeff = [ii(:),jj(:),kk(:),0*kk(:)];
fitCoeff(sum(fitCoeff,2)>4,:) = [];
fitCoeff(max(fitCoeff,[],2) == 4,4) = -1;

parpool_size = 12;
if parpool_size~=0
    pjob = gcp('nocreate');
    if isempty(pjob)
        parpool(parpool_size)
    elseif pjob.NumWorkers ~= parpool_size
        delete(pjob)
        parpool(parpool_size)
    end
end

Th           = 1;
minDist      = 2 / Res * 3;
MaxIter      = 14;
CritIter     = 7;

%%
se           = strel3d(3);
dilatedBW    = imdilate(FinalVol, se);
maxPos       = find( FinalVol == dilatedBW & FinalVol > Th );
maxVals      = FinalVol(maxPos);
[~,sortInd]  = sort(maxVals,'descend');
maxNum       = min(50000, length(sortInd) );
maxPos       = maxPos(sortInd(1:maxNum));

fprintf(1,'numpeak = %d \n',length(maxPos));

maxXYZ = zeros(length(maxPos),3);
for i=1:length(maxPos)
    [xx,yy,zz] = ind2sub(size(FinalVol),maxPos(i));
    maxXYZ(i,:) = [xx yy zz];
end

for ind_a = 1:length(SearchRadAr)
    SearchRad    = SearchRadAr(ind_a);
    savestring = sprintf(ourputstring_prefix,SearchRad);
    disp(savestring)

    %%
    Q = 0.5;
    cropHalfSize = SearchRad;
    [X,Y,Z] = ndgrid(-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize);
    SphereInd = find( X.^2+Y.^2+Z.^2 <= (SearchRad+0.5)^2 );
    X_sph = X(SphereInd);
    Y_sph = Y(SphereInd);
    Z_sph = Z(SphereInd);

    XYZdata.X = X_sph;
    XYZdata.Y = Y_sph;
    XYZdata.Z = Z_sph;

    Orders = fitCoeff(:,1:3);
    PosArr = zeros(size(maxXYZ));
    TotPosArr = zeros(size(maxXYZ));

    exitFlagArr = zeros(1, size(maxXYZ,1));

    CoeffArr = repmat(fitCoeff(:,4),[1 size(maxXYZ,1)]);

    Alpha = 1;

    [box_inten] = get_box_intensity(FinalVol, maxXYZ', cropHalfSize, 1, 0, 'linear');
    parfor i=1:size(maxXYZ,1)
        %         tic
        endFlag = 0;
        consecAccum = 0;
        iterNum = 0;

        cropVol = reshape(box_inten(:,i),(cropHalfSize*2+1).*[1,1,1]);

        while ~endFlag
            iterNum = iterNum + 1;
            if iterNum>MaxIter
                exitFlagArr(i) = -4;
                endFlag = 1;
            end

            Pos = PosArr(i,:);
            GaussWeight = exp(-1*Alpha*( (X_sph-Pos(1)).^2 + (Y_sph-Pos(2)).^2 + (Z_sph-Pos(3)).^2 ) / cropHalfSize^2 );

            fun = @(p,xdata) calculate_3D_polynomial_Rogers(xdata.X,xdata.Y,xdata.Z,Pos,Orders,p).*GaussWeight;

            opts = optimset('Display','off');
            try
                [p1,fminres1] = lsqcurvefit(fun,CoeffArr(:,i),XYZdata,cropVol(SphereInd).*GaussWeight,[],[],opts);
            catch
                endFlag = 1;
                continue;
            end
            CoeffArr(:,i) = p1;

            [dX, dY, dZ] = calc_dX_dY_dZ_Rogers(Orders,CoeffArr(:,i));
            if dX ==-100 && dY == -100 && dZ == -100
                exitFlagArr(i) = -1;
                endFlag = 1;
            else
                confinedShift = min(max([dX dY dZ],-1*[Q Q Q]),[Q Q Q]);
                PosArr(i,:) = PosArr(i,:) + confinedShift;
                if max(abs(PosArr(i,:))) > cropHalfSize
                    exitFlagArr(i) = -2;
                    endFlag = 1;
                elseif max(abs(confinedShift)) < Q
                    if consecAccum == CritIter-1  %|| max(abs(confinedShift)) < 1e-5
                        %                         goodAtomTotPos = TotPosArr(1:i-1,:);
                        %                         goodAtomTotPos = goodAtomTotPos(exitFlagArr(1:i-1)==0,:);
                        %                         %                     Dist = sqrt(sum((goodAtomTotPos - repmat(PosArr(i,:)+maxXYZ(i,:),[size(goodAtomTotPos,1) 1])).^2,2));
                        %                         Dist = pdist2(PosArr(i,:)+maxXYZ(i,:),goodAtomTotPos);
                        %                         if min(Dist) < minDist
                        %                             exitFlagArr(i) = -3;
                        %                         else
                        TotPosArr(i,:) = PosArr(i,:) + maxXYZ(i,:);
                        %                         end
                        endFlag = 1;
                    else
                        consecAccum = consecAccum + 1;
                    end
                else
                    consecAccum = 0;
                end
            end
        end
        fprintf('peak %d, flag %d, consecAccum %d\n',i,exitFlagArr(i),consecAccum);
        %         toc
    end
    for i = 1:size(maxXYZ,1)
        if exitFlagArr(i) == 0
            goodAtomTotPos = TotPosArr(1:i-1,:);
            goodAtomTotPos = goodAtomTotPos(exitFlagArr(1:i-1)==0,:);
            %                     Dist = sqrt(sum((goodAtomTotPos - repmat(PosArr(i,:)+maxXYZ(i,:),[size(goodAtomTotPos,1) 1])).^2,2));
            Dist = pdist2(PosArr(i,:)+maxXYZ(i,:),goodAtomTotPos);
            if min(Dist) < minDist
                exitFlagArr(i) = -3;
            end
        end
    end
    save(sprintf('%s_result.mat',savestring),'PosArr','TotPosArr','CoeffArr','Orders','exitFlagArr');
    % load(sprintf('%s_result.mat',savestring),'PosArr','TotPosArr','CoeffArr','Orders','exitFlagArr');

    atom_pos = TotPosArr(exitFlagArr==0,:)';
    classify_info.halfSize = SearchRad;

    b1 = find(atom_pos(1,:)<15 | atom_pos(1,:)>size(FinalVol,1)-15);
    b2 = find(atom_pos(2,:)<15 | atom_pos(2,:)>size(FinalVol,2)-15);
    b3 = find(atom_pos(3,:)<15 | atom_pos(3,:)>size(FinalVol,3)-15);

    bT = union(union(b1,b2),b3);
    atom_pos(:,bT) = [];

    [temp_model, temp_atomtype] = initial_class_kmean(...
        FinalVol, atom_pos, classify_info);

    [peak_info,intensity_arr,intensity_plot_arr] = ...
        plot_class_hist(FinalVol, temp_model, temp_atomtype, classify_info);

    temp_model_o = temp_model / 3 - 2 ;

    save(sprintf('%s_result.mat',savestring),'temp_model_o',...
        'temp_atomtype','peak_info','intensity_arr','intensity_plot_arr',...
        'PosArr','TotPosArr','CoeffArr','Orders','exitFlagArr');
end