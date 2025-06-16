function Angles_opt = RESIRE_multiGPU_AngleRefine_GetOpAng(RESIRE,Para_Agnle_refine)
%%
addpath('src\')
%% 
kx=Para_Agnle_refine.kx
ky=Para_Agnle_refine.ky
kz=Para_Agnle_refine.kz

dx=Para_Agnle_refine.dx
dy=Para_Agnle_refine.dy
dz=Para_Agnle_refine.dz

cx=round(kx/2)
cy=round(ky/2)
cz=round(kz/2)

plotflag=Para_Agnle_refine.plotflag;
%%
angles_ref=RESIRE.OBJ.ParaInput.angles;
projections=RESIRE.OBJ.ParaInput.projections;
reconstruction=RESIRE.OBJ.reconstruction;

%%
%mkdir('Angles_temp')

folder = 'Angles_temp';
if ~exist(folder, 'dir')
    mkdir(folder);
end

Path = './Angles_temp/';

for j=1:kx
    for i=1:ky
        for k=1:kz
            angles(:,1)=angles_ref(:,1)+dx*(j-cx);
angles(:,2)=angles_ref(:,2)+dy*(i-cy);
angles(:,3)=angles_ref(:,3)+dz*(k-cz);
        eval(['angles_param' num2str(j) num2str(i) num2str(k) ' =angles;']);
        Matrixname=['angles_param',num2str(j) num2str(i) num2str(k)];
   eval(['save ',Path, Matrixname,' ',Matrixname]); 


Filename=[Path,'Resire_CoPdPt_AngleOPT',num2str(j) num2str(i) num2str(k)]

ParaInput_InitialRec_temp.projections=projections;
ParaInput_InitialRec_temp.angles=angles;
ParaInput_InitialRec_temp.step_size=0;
ParaInput_InitialRec_temp.iterations=1;
ParaInput_InitialRec_temp.dimz_scale=1;
ParaInput_InitialRec_temp.plotflag=0;

ParaInput_InitialRec_temp.InitialRec = reconstruction;

OBJ_InitialRec_temp=RESIRE_multiGPU_InitialRec(ParaInput_InitialRec_temp);

OBJ_angleRefine.refinedAngles=angles;
OBJ_angleRefine.Rarr=OBJ_InitialRec_temp.Rarr;
OBJ_angleRefine.R1=OBJ_InitialRec_temp.R1;

save(Filename, "OBJ_angleRefine", '-v7.3')

clear angles;
        end
    end
end
%%
addpath('./Angles_temp')

for j=1:kx
    for i=1:ky
        for k=1:kz

       dataname=['Resire_CoPdPt_AngleOPT' num2str(j) num2str(i) num2str(k) '.mat']
    Resire=load(dataname);

        R(:,j,i,k)=(Resire.OBJ_angleRefine.Rarr);

        Angle_Input_x(:,j,i,k)=(Resire.OBJ_angleRefine.refinedAngles(:,1));
        Angle_Input_y(:,j,i,k)=(Resire.OBJ_angleRefine.refinedAngles(:,2));
        Angle_Input_z(:,j,i,k)=(Resire.OBJ_angleRefine.refinedAngles(:,3));


        end
    end
end
%%
[nx ny nz nk]=size(R)

Angle_Input_x_1D_temp=reshape(Angle_Input_x(round(nx/2),:,:,:),ny*nz*nk,1);
Angle_Input_y_1D_temp=reshape(Angle_Input_y(round(nx/2),:,:,:),ny*nz*nk,1);
Angle_Input_z_1D_temp=reshape(Angle_Input_z(round(nx/2),:,:,:),ny*nz*nk,1);
R_1D=reshape(R(round(nx/2),:,:,:),ny*nz*nk,1);
%%
for i=1:nx
R_1D_temp=reshape(R(i,:,:,:),ny*nz*nk,1);
[R_min_temp,ID_temp]=min(R_1D_temp);
R_min(i)=R_min_temp;
ID_temp;

Angle_Input_x_1D_temp=reshape(Angle_Input_x(i,:,:,:),ny*nz*nk,1);
Angle_Input_y_1D_temp=reshape(Angle_Input_y(i,:,:,:),ny*nz*nk,1);
Angle_Input_z_1D_temp=reshape(Angle_Input_z(i,:,:,:),ny*nz*nk,1);

Angle_opt_x(i)=Angle_Input_x_1D_temp(ID_temp);
Angle_opt_y(i)=Angle_Input_y_1D_temp(ID_temp);
Angle_opt_z(i)=Angle_Input_z_1D_temp(ID_temp);
end

%%
Angles_opt=[Angle_opt_x',Angle_opt_y',Angle_opt_z'];
%%
if plotflag==1
figure
scatter3(Angle_Input_x_1D_temp,Angle_Input_y_1D_temp,Angle_Input_z_1D_temp,R_1D.^4*10^6,'filled')
%
figure
plot(R_min)
%
fig=figure
plot(Angle_opt_z-angles_ref(:,3)','-o')
hold on
plot(Angle_opt_y-angles_ref(:,2)','-o')
plot(Angle_opt_x-angles_ref(:,1)','-o')
legend('z','y','x')
print(fig,"Angle_diff","-dpng")
end

%folder = 'Angles_temp';;
delete(fullfile(folder, '*.mat'));

end
