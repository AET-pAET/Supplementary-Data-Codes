clear all; close all
%%
addpath('src')
addpath('../4.1_projections')
%%
projections = importdata('Proj_CoPdPt_NoRadmDefocus_cov35pm_Dose1p7e4_BM3DOPT.mat');
angles_GT=load('CoPdPt_Angle_Dose1p7e4.mat');% Ground truth
angles_GT=angles_GT.angles_refined_out;
angles      = importdata('CoPdPt_ExperimentalAngles.mat');% Angle from experiment
angles_OPT      = importdata('CoPdPt_RefinedAngles.mat'); % Refined Angles
angles_OPT=angles_OPT.Angle_OPT;
%%
ParaInput.projections=projections;
ParaInput.angles=angles;% Run reconstruction with experimental measured angle
%ParaInput.angles=angles_OPT; Run reconstruction with optimized angle
ParaInput.step_size=1;
ParaInput.iterations=1000;
ParaInput.dimz_scale=1;
ParaInput.plotflag=1;

OBJ = RESIRE_multiGPU(ParaInput);
%% Output
%save("Resire_CoPdPt_NoRadmDefocus_BM3D_ExpAngle_S1p0N1000.mat","OBJ",'-v7.3') % Angle from experiment
%save("Resire_CoPdPt_NoRadmDefocus_BM3D_RefinedAngle_S1p0N1000.mat","OBJ",'-v7.3')% refined Angle


