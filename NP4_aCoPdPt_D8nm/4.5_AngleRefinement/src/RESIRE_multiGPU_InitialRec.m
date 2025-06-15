function OBJ = RESIRE_multiGPU_InitialRec(ParaInput)
%% read projections after BM3D and tilt angles and parameters
projections=ParaInput.projections;
[dimx, dimy, Num_pj] = size(projections);
angles=ParaInput.angles;
step_size      = ParaInput.step_size;  
iterations     = ParaInput.iterations;
dimz_scale     = ParaInput.dimz_scale;
InitialRec     = ParaInput.InitialRec;

plotflag     = ParaInput.plotflag;

dimz           = dimz_scale*dimx;
positivity     = true;
%% input
rotation       = 'ZYX';  % Euler angles setting ZYZ
dtype          = 'single';
projections_refined = cast(projections,dtype);
angles_refined      = cast(angles,dtype);

% compute normal vector of rotation matrix
matR = zeros(3,3);
if length(rotation)~=3
    disp('rotation not recognized. Set rotation = ZYX\n'); rotation = 'ZYX';
end
for i=1:3
    switch rotation(i)
        case 'X',   matR(:,i) = [1;0;0];
        case 'Y',   matR(:,i) = [0;1;0];
        case 'Z',   matR(:,i) = [0;0;1];
        otherwise,  matR = [0,0,1;
                0,1,0;
                1,0,0];
            disp('Rotation not recognized. Set rotation = ZYX');
            break
    end
end
vec1 = matR(:,1); vec2 = matR(:,2); vec3 = matR(:,3);
%% rotation matrix
Rs = zeros(3,3,Num_pj, dtype);
for k = 1:Num_pj
    phi   = angles_refined(k,1);
    theta = angles_refined(k,2);
    psi   = angles_refined(k,3);
    
    % compute rotation matrix R w.r.t euler angles {phi,theta,psi}
    rotmat1 = MatrixQuaternionRot(vec1,phi);
    rotmat2 = MatrixQuaternionRot(vec2,theta);
    rotmat3 = MatrixQuaternionRot(vec3,psi);
    R =  single(rotmat1*rotmat2*rotmat3)';
    Rs(:,:,k) = R;
end
%% Run RESIRE with multi-GPU
fprintf('\nResire code:(multi-GPU)\n');
dim_ext = [dimx,dimy,dimz];
tic
[rec] = RT3_film_multiGPU( (projections_refined), (Rs), dim_ext, ...
    iterations, (step_size) , (positivity),InitialRec);
toc
%% Cal. Projections
cal_projs = calculate3Dprojection_multiGPU(single(rec), Rs);
%% Cal. R-factor
for i=1:Num_pj
    pj = projections(:,:,i);
    resi_i=projections(:,:,i)-cal_projs(:,:,i);
    Rarr(i) = sum(abs(resi_i(:)))/ sum(abs(pj(:)));

    pj_cal = cal_projs(:,:,i);
    pj_cal = ImageNorm(pj_cal,sum(sum(pj)));
    resi_i_norm=pj-pj_cal;
    Rarr_norm(i) = sum(abs(resi_i_norm(:)))/ sum(abs(pj(:)));
    clear resi_i  pj  pj_cal resi_i_norm
end

R=mean(Rarr)
R_norm=mean(Rarr_norm)
%% Output
OBJ.ParaInput=ParaInput;
OBJ.reconstruction=rec;
OBJ.cal_projs=cal_projs;

OBJ.cal_projs;
OBJ.Rarr=Rarr;
OBJ.R1=R;
OBJ.Rarr_norm=Rarr_norm;
OBJ.R1_norm=R_norm;
%%
if plotflag ==1
figure();img(cal_projs,'colormap','gray')
figure();img(rec,'colormap','gray')
end

end
