clear;
close all;

close all; 
clear;

addpath('../../Plot');

cases= ["villasor","ky1","ferto","ky2","sanchegy","ky3","buk","ky4","lovo","ky5","nagycenk","ky6","vashegy","ky7","varis","ky8",...
    "becsidomb","ky9","tomalom","ky10","szakov","ky11","kohegy","ky12","harka","ky13","pozsonyiut","ky14"];
 
idx = 1:28;
idx_iso = 1:27;
idx_net = 1:27;

%blackBody, blackBodyExt, cividis, coolWarmBent, coolWarmSmooth, inferno, jet, kindlmann, kindlmannExt, magma, plasma, viridis
%discrete: lines, prism
% colorMapName = 'grayscale'; 
colorMapName = 'plasma';
colorMap = importdata(['../../Plot/ColorMaps/',colorMapName,'.col']);

fig = figure('Position',[200 200 900 600]);
hold on; grid on;
% moderately vulnerable
plot([1e-3,1e-3],[1e-1,1e7],'--k','linewidth',1.5,'HandleVisibility','off')
annotation('textarrow',[0.647,0.717],[0.71,0.71],'FontSize',16,'Linewidth',1.5)
text(1e-2,3e6,'Moderately','FontSize',14,'HorizontalAlignment','center');
text(1e-2,1.2e6,'vulnerable','FontSize',14,'HorizontalAlignment','center');
% text(1e-2,3e6,'Mérsékelten','FontSize',14,'HorizontalAlignment','center');
% text(1e-2,1.2e6,'sebezhető','FontSize',14,'HorizontalAlignment','center');
% highly vulnerable
plot([1e-2,1e-2],[1e-1,1e5],'--k','linewidth',1.5,'HandleVisibility','off')
annotation('textarrow',[0.732,0.802],[0.53,0.53],'FontSize',16,'Linewidth',1.5)
text(1e-1,3e4,'Highly','FontSize',14,'HorizontalAlignment','center');
text(1e-1,1.2e4,'vulnerable','FontSize',14,'HorizontalAlignment','center');
% text(1e-1,3e4,'Kritikusan','FontSize',14,'HorizontalAlignment','center');
% text(1e-1,1.2e4,'sebezhető','FontSize',14,'HorizontalAlignment','center');
for I=idx
    if(mod(I,2)==0)
        data = importdata(join(['Network Data/',cases(I),'/vulner_N.txt'],''));
        col = 0.25;
    else
        data = importdata(join(['Network Data/',cases(I),'/vulner_N.txt'],''));
        col = 0.75;
    end

    set(gca,'XScale','log','YScale','log');

    gammaBar = data;
    gammaBarNZ = gammaBar(gammaBar~=0);
    n = length(gammaBarNZ);
    if(n>100)
        r = round(log(n)/log(2)+1);
    else
        r = round(sqrt(n));
    end
    b = prctile(gammaBarNZ,(0:1/r:1)*100);
    x = (b(1:end-1)+b(2:end))/2;

    f = zeros(r,1);
    for j=1:r
        f(j) = length(gammaBarNZ(gammaBarNZ>=b(j) & gammaBarNZ<b(j+1)));
    end
    f = f/sum(f); % relative frequency

    y = f./diff(b');

    colr = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),col,'spline')));
    colg = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),col,'spline')));
    colb = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),col,'spline')));

    if(mod(I,2) == 0)
        marker = '+';
    else
        marker = 'o';
    end
    plot(x,y,marker,'color',[colr,colg,colb],'linewidth',1.5);
end
xlabel('Local vulnerability [-]','fontsize',14);
ylabel('Probability density [-]','fontsize',14);
% xlabel('Lokális sebezhetőség [-]','fontsize',14);
% ylabel('Valószínűségi sűrűség [-]','fontsize',14);
% legend(cases(idx));
xlim([1e-9,1e0])
ylim([1e-1,1e8])
% set(gca,'Position',[100,100,600,400])
% legend('Eredeti','N szabĂˇly','N-1 szabĂˇly');
legend('Real-life WDNs: SN*','Artificial WDNs: ky*','location','southwest');
% set(gca,'FontSize',14);
% ColumnLegend(3,num2str(idx(:)));
% ColumnLegend(3,num2str(idx(:)),'location','southwest_sf');
% legend('Original','N rule','N-1 rule');
% rectangle('Position',[1.3e-9 1.9e-1 7.5e-7 5e3],'FaceColor',[1 1 1])
% txt1 = text(2e-8,3e3,'Ivóvízhálózat','FontSize',14,'HorizontalAlignment','center');
% txt2 = text(2e-8,1.2e3,'sorszám','FontSize',14,'HorizontalAlignment','center');
% txt1 = text(2e-8,3e3,'Water distribution','FontSize',14,'HorizontalAlignment','center');
% txt2 = text(2e-8,1.2e3,'network','FontSize',14,'HorizontalAlignment','center');
% title('Failure rate: relative pipe length, ky: N-1 rule');
set(gca,'FontSize',14);
% title('Failure rate: pipe material');
saveas(gca,'Plots/Vulner_all_N.fig','fig');
saveas(gca,'Plots/Vulner_all_N.png','png');
saveas(gca,'Plots/Vulner_all_N.eps','epsc');


