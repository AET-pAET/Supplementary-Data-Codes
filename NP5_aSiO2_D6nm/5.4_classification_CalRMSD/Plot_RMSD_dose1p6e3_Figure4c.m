clear all; close all;
%%
RMSD_D1p6e3_cov35pm_multi=load('RMSD_Model_aSiO2_Dose1p6e3_multislice_cov35pm_S1p0N1000.mat');
RMSD_D1p6e3_cov35pm_multi_Si=RMSD_D1p6e3_cov35pm_multi.dif_Si;
RMSD_D1p6e3_cov35pm_multi_O=RMSD_D1p6e3_cov35pm_multi.dif_O;
RMSD_D1p6e3_cov35pm_multi_All=RMSD_D1p6e3_cov35pm_multi.dif_all;
%%
RMSD_D1p6e3_cov35pm_single=load('RMSD_Model_aSiO2_Dose1p6e3_single_slice_cov35pm_S2p0N1000.mat');
RMSD_D1p6e3_cov35pm_single_Si=RMSD_D1p6e3_cov35pm_single.dif_Si;
RMSD_D1p6e3_cov35pm_single_O=RMSD_D1p6e3_cov35pm_single.dif_O;
RMSD_D1p6e3_cov35pm_single_All=RMSD_D1p6e3_cov35pm_single.dif_all;
%%

figure(11); clf; set(gcf,'Position',[0,0,230*3,230]); 

subplot(1,3,2);
hh=histogram(pdist2(double(RMSD_D1p6e3_cov35pm_multi_O'),[0,0,0]), 'BinEdges', (0:0.03:1.0)*0.6, 'EdgeColor', 'black', ...
'FaceColor', [0 0 160]/255, 'FaceAlpha', 0.5,'linewidth',0.6);
ylabel('Number of atoms')

hold on
hh=histogram(pdist2(double(RMSD_D1p6e3_cov35pm_single_O'),[0,0,0]), 'BinEdges', (0:0.03:1.0)*0.6, 'EdgeColor', 'black', ...
'FaceColor', [220 111 0]/255, 'FaceAlpha', 0.5,'linewidth',0.6);
xlim([0,0.6])
ylim([0,700])
xlabel(['Deviation (',char(197),')'])
ylabel('Number of atoms')
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);


subplot(1,3,1);
hh=histogram(pdist2(double(RMSD_D1p6e3_cov35pm_multi_Si'),[0,0,0]), 'BinEdges', (0:0.03:1.0)*0.8/0.4*0.2, 'EdgeColor', 'black', ...
'FaceColor', [0 0 160]/255, 'FaceAlpha', 0.5,'linewidth',0.6);
ylabel('Number of atoms')

hold on
hh=histogram(pdist2(double(RMSD_D1p6e3_cov35pm_single_Si'),[0,0,0]), 'BinEdges', (0:0.03:1.0)*0.8/0.4*0.2, 'EdgeColor', 'black', ...
'FaceColor', [220 111 0]/255, 'FaceAlpha', 0.5,'linewidth',0.6);
xlim([0,0.4])
xlabel(['Deviation (',char(197),')'])
ylabel('Number of atoms')
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
ylim([0,500])
lgd=legend('Multislice','Single-slice');
legend('boxoff') 
lgd.ItemTokenSize = [20, 6]
legend('boxoff') 
set(legend,'FontSize',8,'linewidth',1.0)
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
set(legend,'Location','NorthWest') 

subplot(1,3,3);
hh=histogram(pdist2(double(RMSD_D1p6e3_cov35pm_multi_All'),[0,0,0]), 'BinEdges', (0:0.03:1)*0.6, 'EdgeColor', 'black', ...
'FaceColor', [0 0 160]/255, 'FaceAlpha', 0.5,'linewidth',0.6);
ylabel('Number of atoms')

hold on
hh=histogram(pdist2(double(RMSD_D1p6e3_cov35pm_single_All'),[0,0,0]), 'BinEdges', (0:0.03:1)*0.6, 'EdgeColor', 'black', ...
'FaceColor', [220 111 0]/255, 'FaceAlpha', 0.5,'linewidth',0.6);
xlim([0,0.6])
xlabel(['Deviation (',char(197),')'])
ylabel('Number of atoms')
ylim([0,1300])

set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);


% print('-dbmp', '-r1200', 'Fig_4c.bmp');
% print('-djpeg', '-r1200', 'Fig_4c.jpg');












