clear all; clear all
%%
RMSD_D1p7e4_cov35pm_WithoutRadmDefocus=load('RMSD_Model_CoPdPt_NoRadmDefocus_cov35pmDose1p7e4_S1N1000.mat');
RMSD_D1p7e4_cov35pm_WithoutRadmDefocus_Co=RMSD_D1p7e4_cov35pm_WithoutRadmDefocus.dif_Co;
RMSD_D1p7e4_cov35pm_WithoutRadmDefocus_Pd=RMSD_D1p7e4_cov35pm_WithoutRadmDefocus.dif_Pd;
RMSD_D1p7e4_cov35pm_WithoutRadmDefocus_Pt=RMSD_D1p7e4_cov35pm_WithoutRadmDefocus.dif_Pt;
%%
RMSD_D1p7e4_cov35pm_WithRadmDefocus=load('RMSD_Model_CoPdPt_RadmDefocus_cov35pmDose1p7e4_S2N1000.mat');
RMSD_D1p7e4_cov35pm_WithRadmDefocus_Co=RMSD_D1p7e4_cov35pm_WithRadmDefocus.dif_Co;
RMSD_D1p7e4_cov35pm_WithRadmDefocus_Pd=RMSD_D1p7e4_cov35pm_WithRadmDefocus.dif_Pd;
RMSD_D1p7e4_cov35pm_WithRadmDefocus_Pt=RMSD_D1p7e4_cov35pm_WithRadmDefocus.dif_Pt;
%% with Random defocus
figure; clf; set(gcf,'Position',[0,0,210*4,210]);
subplot(1,3,2);
hold on
hh=histogram(pdist2(double(RMSD_D1p7e4_cov35pm_WithRadmDefocus_Pd'),[0,0,0]), 'BinEdges', (0:0.03:1.0)*0.5*0.8, 'EdgeColor', 'black', ...
'FaceColor', [220 111 0]/255, 'FaceAlpha', 0.5,'linewidth',0.6);
xlim([0,0.4])
ylim([0,700])
xlabel(['Deviation (',char(197),')'])
ylabel('Number of atoms')
set(gca,'FontSize',11,'FontName', 'Arial','linewidth',1.0);
box on

subplot(1,3,1);
hold on
hh=histogram(pdist2(double(RMSD_D1p7e4_cov35pm_WithRadmDefocus_Co'),[0,0,0]), 'BinEdges', (0:0.03:1.0)*0.5*0.8/0.4*0.8, 'EdgeColor', 'black', ...
'FaceColor', [220 111 0]/255, 'FaceAlpha', 0.5,'linewidth',0.6);

xlim([0,0.8])
xlabel(['Deviation (',char(197),')'])
ylabel('Number of atoms')
set(gca,'FontSize',11,'FontName', 'Arial','linewidth',1.0);
ylim([0,600])
set(gca,'FontSize',11,'FontName', 'Arial','linewidth',1.0);
box on


subplot(1,3,3);
hold on
hh=histogram(pdist2(double(RMSD_D1p7e4_cov35pm_WithRadmDefocus_Pt'),[0,0,0]), 'BinEdges', (0:0.015:0.5)*0.4, 'EdgeColor', 'black', ...
'FaceColor', [220 111 0]/255, 'FaceAlpha', 0.5,'linewidth',0.6);
xlim([0,0.2])
xlabel(['Deviation (',char(197),')'])
ylabel('Number of atoms')
ylim([0,300])
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
box on
%% without Random defocus
figure; clf; set(gcf,'Position',[0,0,210*4,210]);
subplot(1,3,2);
hh=histogram(pdist2(double(RMSD_D1p7e4_cov35pm_WithoutRadmDefocus_Pd'),[0,0,0]), 'BinEdges', (0:0.03:1.0)*0.5*0.8, 'EdgeColor', 'black', ...
'FaceColor', [0 0 160]/255, 'FaceAlpha', 0.5,'linewidth',0.6);
%xlim([0,1.2])
xlabel(['Deviation (',char(197),')'])
ylabel('Number of atoms')
set(gca,'FontSize',11,'FontName', 'Arial','linewidth',1.0);

subplot(1,3,1);
hh=histogram(pdist2(double(RMSD_D1p7e4_cov35pm_WithoutRadmDefocus_Co'),[0,0,0]), 'BinEdges', (0:0.03:1.0)*0.5*0.8/0.4*0.8, 'EdgeColor', 'black', ...
'FaceColor', [0 0 160]/255, 'FaceAlpha', 0.5,'linewidth',0.6);

xlim([0,0.8])
xlabel(['Deviation (',char(197),')'])
ylabel('Number of atoms')
set(gca,'FontSize',11,'FontName', 'Arial','linewidth',1.0);
ylim([0,900])

subplot(1,3,3);
hh=histogram(pdist2(double(RMSD_D1p7e4_cov35pm_WithoutRadmDefocus_Pt'),[0,0,0]), 'BinEdges', (0:0.015:0.5)*0.4, 'EdgeColor', 'black', ...
'FaceColor', [0 0 160]/255, 'FaceAlpha', 0.5,'linewidth',0.6);
ylabel('Number of atoms')
xlim([0,0.2])
xlabel(['Deviation (',char(197),')'])
ylabel('Number of atoms')
ylim([0,400])

set(gca,'FontSize',11,'FontName', 'Arial','linewidth',1.0);

print('-dbmp', '-r1200', 'Fig_3e.bmp');
print('-djpeg', '-r1200', 'Fig_3e.jpg');