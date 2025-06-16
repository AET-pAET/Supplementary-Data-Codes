% Polynomial tracing based on reconstruction of original orientation
clear all; close all
%%
addpath('./src_classification');
%% load rec data
RESIRE=load('3_Resire_Dose5p6e4_BM3D_After_Alignment.mat')
Dsetvol=RESIRE.OBJ.reconstruction*10*50;

N = size(Dsetvol)
fprintf('Recstruction Dimension = [%d,%d,%d];\n',N(1),N(2),N(3));
%return
Dsetvol0 = My_stripzero(Dsetvol,[N(1),N(2),N(3)]);
Dsetvol_1 = Dsetvol0(:,1:end,:);

FinalVol_tracing = (Dsetvol_1);
N1=size(FinalVol_tracing);
fprintf('Tracing Dimension = [%d,%d,%d];\n',N1(1),N1(2),N1(3));
%% oversamplling the reconstruction matrix by 3*3*3 by linear interpolation
ov=3;% oversamplling
xx = (1:size(FinalVol_tracing,1)) - round((size(FinalVol_tracing,1)+1)/2);
yy = (1:size(FinalVol_tracing,2)) - round((size(FinalVol_tracing,2)+1)/2);
zz = (1:size(FinalVol_tracing,3)) - round((size(FinalVol_tracing,3)+1)/2);

xxi = ((ov*xx(1)):(xx(end)*ov))/ov;
yyi = ((ov*yy(1)):(yy(end)*ov))/ov;
zzi = ((ov*zz(1)):(zz(end)*ov))/ov;

xxi = xxi(3:end); yyi = yyi(3:end); zzi = zzi(3:end);
[Y,X,Z]     = meshgrid(yy,xx,zz);
[Yi,Xi,Zi]  = meshgrid(yyi,xxi,zzi);

FinalVol_tracing=single(FinalVol_tracing);
FinalVol_tracing_interp     = interp3(Y,X,Z,FinalVol_tracing,Yi,Xi,Zi,'spline',0);
FinalVol    = My_paddzero(FinalVol_tracing_interp,size(FinalVol_tracing_interp)+20);% 

FinalVol=double(FinalVol);
N2=size(FinalVol);
fprintf('Tracing Dimension (oversamplling) = [%d,%d,%d];\n',N2(1),N2(2),N2(3));
%% load tracing data
load('3_polyn_tracing_BM3D_After_Alignment_result.mat') %
atom_pos0 = TotPosArr(exitFlagArr==0,:); %
atom_pos1=atom_pos0(1:end,:)-0;
%% initial speration of type-2 (Pd) and type-3 (Pt)

classify_info1 = struct('Num_species', 4,  'halfSize',  3,  'plothalfSize',  1, 'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  120, 'lnorm',2);

[temp_model1, temp_atomtype1] = initial_class_kmean_sub(FinalVol, atom_pos1', classify_info1);  
%% Classification of non-atoms and Co
ID_type4=find(temp_atomtype1==4);
ID_type3=find(temp_atomtype1==3);
ID_type2=find(temp_atomtype1==2);
ID_type2_type4=[ID_type4 ID_type3 ID_type2];

atom_pos2=atom_pos1(:,:);
atom_pos2(ID_type2_type4,:)=[];
atom_pos_type2_type4=atom_pos1(ID_type2_type4,:);

classify_info2 = struct('Num_species', 2,  'halfSize',  3,  'plothalfSize',  1, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  120, 'lnorm',2);

[temp_model2, temp_atomtype2,score] = initial_classification1_1types_fixnonatom(FinalVol, atom_pos2', classify_info2,0.00);
%% remove non-atoms
ID_type1=find(temp_atomtype2==1);
atom_pos_type1=atom_pos2(:,:);
atom_pos_type1(ID_type1,:)=[];
%% Classification of type-1 (Co), type-2 (Pd) type-3 (Pt)
Atom_pos_all1=[atom_pos_type2_type4; atom_pos_type1];

classify_info = struct('Num_species', 4,  'halfSize',  3,  'plothalfSize',  1, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  120, 'lnorm',2);

[temp_model0, temp_atomtype0] = initial_class_kmean_sub(FinalVol, Atom_pos_all1', classify_info);
%% merge type-3 and type-3 as Pt
ID4=find(temp_atomtype0==4);
temp_atomtype0(ID4)=3;
ID4_check=find(temp_atomtype0==4);
ID3=find(temp_atomtype0==3);
%% Output
temp_model=[temp_model0];
temp_atomtype=[temp_atomtype0];

%save('Model_CoPdPt_Dose5p6e4_BM3D_RESIRE.mat','temp_atomtype','temp_model');
%% Plot hist Extended Data Figure5g
classify_info0 = struct('Num_species', 3,  'halfSize',  3,  'plothalfSize',  3, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  150, 'lnorm',2);

[peak_info,intensity_arr,intensity_plot_arr] = plot_class_hist(FinalVol, temp_model, temp_atomtype, classify_info0);
%
figure(18); clf; hold on; 
set(gcf,'position',[250,250,750,500]);
bar(peak_info(1,:)/3.5/10/50,peak_info(3,:),'facecolor',[0,0.3,0.7],'linewidth',0.5);
bar(peak_info(1,:)/3.5/10/50,peak_info(4,:),'facecolor',[0,0.7,0.3],'linewidth',0.5);
bar(peak_info(1,:)/3.5/10/50,peak_info(5,:),'facecolor',[0.7,0.7,0],'linewidth',0.5);
box on; 

set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
ylabel(['Number of atoms'])
xlabel('Intensity (a.u.)')
xticks(0:500:2000)
