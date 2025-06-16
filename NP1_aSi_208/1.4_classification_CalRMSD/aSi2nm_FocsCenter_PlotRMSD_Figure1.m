clear all; close all
%% load displacement from the GT
dif_Dose8p0e3=load('dif_Dose8p0e3_AtomTypeAll.mat');
dif_Dose8p0e3=dif_Dose8p0e3.dif;
dif_Dose1p6e4=load('dif_Dose1p6e4_AtomTypeAll.mat');
dif_Dose1p6e4=dif_Dose1p6e4.dif;
dif_Dose1p6e5=load('dif_Dose1p6e5_AtomTypeAll.mat');
dif_Dose1p6e5=dif_Dose1p6e5.dif;
%% plot RMSD
figure; clf; set(gcf,'Position',[200,200,240,180]);
set(gca, 'Color', 'none');
set(gcf, 'Color', 'white');
hh=histogram(pdist2(double(dif_Dose1p6e5'),[0,0,0]), 'BinEdges', (0:0.015:0.5), 'EdgeColor', 'black', ...
'FaceColor', [0 0 160]/255, 'FaceAlpha', 0.5,'linewidth',0.6);
hold on
hh=histogram(pdist2(double(dif_Dose1p6e4'),[0,0,0]), 'BinEdges', (0:0.015:0.5), 'EdgeColor', 'black', ...
'FaceColor', [0 102 102]/255, 'FaceAlpha', 0.5,'linewidth',0.6);
hold on
hh=histogram(pdist2(double(dif_Dose8p0e3'),[0,0,0]), 'BinEdges', (0:0.015:0.5), 'EdgeColor', 'black', ...
'FaceColor', [220 111 0]/255, 'FaceAlpha', 0.5,'linewidth',0.6);

ylim([0,80])
xlim([0,0.5])
xlabel(['Deviation (',char(197),')'])
ylabel('Number of atoms')
lgd = legend('1.6×10^5 e^-/Å^2','1.6×10^4 e^-/Å^2','8.0×10^3 e^-/Å^2');
lgd.Box = 'off'; 
lgd.ItemTokenSize = [15, 5]
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
copygraphics(gcf, 'Resolution', 1200); 
xticks(0:0.1:0.4)
%print('-djpeg', '-r1200', 'Fig_1d.jpg');
print('-dbmp', '-r1200', 'Fig_1d.bmp');
%%





