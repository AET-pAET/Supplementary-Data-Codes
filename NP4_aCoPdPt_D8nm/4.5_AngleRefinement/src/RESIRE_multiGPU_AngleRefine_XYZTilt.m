function RESIRE_OPT = RESIRE_multiGPU_AngleRefine_XYZTilt(RESIRE_inital,kx,ky,kz,dx,dy,dz,steps)
%%
RESIRE_ref=RESIRE_inital;%load('Resire_CoPdPt_Dose1p7e4_bg1_corp_S2N300_ExpAng_DAngle0p1.mat');
RESIRE_ref.OBJ.ParaInput.angles

plotflag=0;

Para_Agnle_refine.kx = kx
Para_Agnle_refine.ky = ky
Para_Agnle_refine.kz = kz

Para_Agnle_refine.dx = dx
Para_Agnle_refine.dy = dy
Para_Agnle_refine.dz = dz

Para_Agnle_refine.plotflag=1;
RESIRE=RESIRE_ref;
%return
%%
for inter=1:steps
inter
    if mod(inter-1,1) == 0
        InitialRec=RESIRE.OBJ.reconstruction;
    end

R_inter(inter)=RESIRE.OBJ.R1;

Angles_opt = RESIRE_multiGPU_AngleRefine_GetOpAng(RESIRE,Para_Agnle_refine)

ParaInput.projections=RESIRE.OBJ.ParaInput.projections;
ParaInput.angles=Angles_opt;
ParaInput.step_size=2;
ParaInput.iterations=30;
ParaInput.dimz_scale=1;
ParaInput.plotflag=0;

OBJ = RESIRE_multiGPU(ParaInput);
RESIRE.OBJ=OBJ;


if inter==inter
    RESIRE_OPT=RESIRE;
    RESIRE_OPT.R_inter=R_inter;
end


fig_R=figure
plot(R_inter)
print(fig_R,"R_inter","-dpng")

clear Angles_opt;
end

end
