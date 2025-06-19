function Projs = calcLinearPjs_fit_para(para,xdata)

Angles   = xdata.angles;
model    = xdata.model;
atomtype = xdata.atoms;
Z_arr 	 = xdata.Z_arr;
Res      = xdata.Res;
volSize  = xdata.volSize;

CropHalfWidth = xdata.CropHalfWidth;

fixedfa = make_fixedfa_man(volSize, Res, Z_arr);

model = model ./ Res;
Heights = para(1,:);
FHeights = Heights;
FWidths = para(2,:) / pi^2 / Res^2;

if length(volSize) == 3
    x = (1:volSize(1)) - round((volSize(1)+1)/2);
    y = (1:volSize(2)) - round((volSize(2)+1)/2);
    z = (1:volSize(3)) - round((volSize(3)+1)/2);
elseif length(volSize) == 1
    x = (1:volSize(1)) - round((volSize(1)+1)/2);
    y = x;
    z = x;
else
    error('volSize should be either length 3 or length 1!')
end

sizeX = [length(x) length(y)];

inInd = find(model(1,:) >= min(x) & model(1,:) <= max(x) & ...
    model(2,:) >= min(y) & model(2,:) <= max(y) & ...
    model(3,:) >= min(z) & model(3,:) <= max(z));

calcModel = model(:,inInd);
calcAtomtype = atomtype(:,inInd);

cropIndRef = -CropHalfWidth:CropHalfWidth;
[cropX,cropY] = ndgrid(cropIndRef,cropIndRef);

Projs = zeros(sizeX(1),sizeX(2),size(Angles,1));
for i=1:size(Angles,1)
    disp(num2str(i))
    RM1 = MatrixQuaternionRot([0 0 1],Angles(i,1));
    RM2 = MatrixQuaternionRot([0 1 0],Angles(i,2));
    RM3 = MatrixQuaternionRot([1 0 0],Angles(i,3));

    calcModel_rot = (RM1*RM2*RM3)'*calcModel;

    finalvol_padded = zeros( [sizeX, length(Heights)]);
    cenPos = round((size(finalvol_padded)+1)/2);

    for j = 1:size(calcModel_rot,2)

        currPos = calcModel_rot(1:2,j) + cenPos(1:2)';
        currRndPos = round(currPos);

        cropInd1 = cropIndRef + currRndPos(1);
        cropInd2 = cropIndRef + currRndPos(2);

        CropVol = finalvol_padded(cropInd1,cropInd2,calcAtomtype(j));

        diffPos = currPos-currRndPos;
        diffPosZ = calcModel_rot(3,j)-round(calcModel_rot(3,j));

        gaussCalc = FHeights(calcAtomtype(j))*exp( -1*( (cropX-diffPos(1)).^2 + (cropY-diffPos(2)).^2 )/FWidths(calcAtomtype(j)) );
        gaussZcalc = exp(-1*(cropIndRef - diffPosZ).^2 / FWidths(calcAtomtype(j)) );

        finalvol_padded(cropInd1,cropInd2,calcAtomtype(j)) = CropVol + gaussCalc*sum(gaussZcalc);
    end

    finalvol_summed = zeros( sizeX);
    for j = 1:length(Heights)
        %fa = reshape(fatom_vector( sqrt(q2),AtomNumbers(j)),sizeX);
        CVol = My_stripzero(finalvol_padded(:,:,j),sizeX);
        FVol = My_FFTN(CVol);
        FVol = FVol .*  reshape(fixedfa,sizeX) ;
        finalvol_summed =finalvol_summed+FVol;
    end

    Vol = real(My_IFFTN(finalvol_summed));
    Projs(:,:,i) = sum(Vol,3);
end