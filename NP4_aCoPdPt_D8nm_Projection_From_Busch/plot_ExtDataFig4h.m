close all; clear all;
%%
load('peak_info_SiGeSn_40pA_After_BM3D_Alignment.mat')
figure(15); clf; hold on; set(gcf,'position',[250,250,434/2*1.01,418/2*0.92]);
bar(peak_info(1,:)/100000*4.5,peak_info(3,:),'facecolor',[1,0.3,0.5],'BarWidth',1,'EdgeColor','none');
bar(peak_info(1,:)/100000*4.5,peak_info(4,:),'facecolor',[0,0.7,0.3],'BarWidth',1,'EdgeColor','none');
bar(peak_info(1,:)/100000*4.5,peak_info(5,:),'facecolor',[0,0.6,0.9],'BarWidth',1,'EdgeColor','none');
xlim([0,16]); xticks(0:4:12); box on;
set(gca,'FontSize',10,'FontName', 'Arial','linewidth',1.0);
set(gca, 'Layer', 'top')
%%
print('-r1200','-djpeg','ExtData_Fig4h.jpg');