function plotCommonLine(Image,Axis,varargin)
% Axis = 1 : X-Axis
% Axis = 2 : Y-Axis
figure; set(gcf,'position',[600,200,700,500]);
for i=1:size(Image,3)
%     plot(refXProfile,'g-','linewidth',1);
    if Axis ~= 1 && Axis ~= 2 
        fprintf('Type the right commonline axis! \n');
        break
    end
    hold all;
    if Axis == 1
        plot(1:size(Image,2),sum(Image(:,:,i),1));
    else
        plot(1:size(Image,1),sum(Image(:,:,i),2));
    end
    
    if numel(varargin) > 0
        tPause = varargin{1};
        title(sprintf('ProjNum: %d',i));
        if tPause > 0
            pause(tPause);
        else
            pause();
        end
    end
end
hold off; box on;
set(gca,'linewidth',1.5,'fontsize',12);
end
