clear;
figure;hold on;
%set
markerSize = 10;
files = dir('*.col');
for i=1:length(files)
    colorMap = importdata([files(i).name]);
    for j=1:size(colorMap,1)
        plot(i,colorMap(j,1),'ko','markersize',markerSize,'linewidth',2,'markerfacecolor',[colorMap(j,2),colorMap(j,3),colorMap(j,4)],'color',[colorMap(j,2),colorMap(j,3),colorMap(j,4)])
    end
    legendName{i} = files(i).name(1:end-4);
end

set(gcf, 'Position',  [200, 200, 800, 400]);
legend(legendName);
xlim([0,20]);