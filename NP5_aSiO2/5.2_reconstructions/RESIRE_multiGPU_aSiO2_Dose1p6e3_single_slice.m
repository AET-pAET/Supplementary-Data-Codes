clear all; close all
%%
addpath('..\5.1_projections')
addpath('src\')
%% read projections and tilt angles
projections = importdata('phases_single_slice_1.6e3_SiO2_cov35pm.mat');
angles(1,:)=(-72:3:72);
angles(2,:)= (-72:3:72);
angles(3,:)=(-72:3:72);
angles(3,:)=0;
angles(1,:)=0;
angles=angles';
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

% extract size of projections & num of projections
[dimx, dimy, Num_pj] = size(projections_refined);
%% parameter
step_size      = 2;  
iterations     = 1000;
dimz           = dimx;
positivity     = true;
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
    iterations, (step_size) , (positivity) );
%% Cal. Projections and R-factor
cal_projs = calculate3Dprojection_multiGPU(single(rec), Rs);
figure(), img(cal_projs,'colormap','gray')

[N1, N2,Num_pj]=size(projections);

for i=1:Num_pj
    pj = projections(:,:,i);
    resi_i=projections(:,:,i)-cal_projs(:,:,i);
    Rarr(i) = sum(abs(resi_i(:)))/ sum(abs(pj(:)));
end
R1=mean(Rarr)
%% Output
OBJ.InputProjections=projections;
OBJ.InputAngles=angles;
OBJ.Calprojection=cal_projs;
OBJ.reconstruction=rec;
OBJ.step_size=step_size;
OBJ.Dimx = dimx;
OBJ.Dimy = dimy;
OBJ.Dimz = dimz;
OBJ.numIterations=iterations;
OBJ.Rarr=Rarr;
OBJ.R1=R1;

%%
save("Resire_aSiO2_Dose1p6e3_single_slice_cov35pm_S2p0N1000.mat","OBJ",'-v7.3')



