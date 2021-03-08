close all; 
clear;

addpath('../../Plot');

% cases= ["villasor","ferto","sanchegy","buk","lovo","nagycenk","vashegy","varis","becsidomb","tomalom",...
%     "szakov","kohegy","harka","pozsonyiut","sopronkovesd","dudlesz","ivan","agyagosszergeny","kofejto","simasag",...
%     "acsad","csaford","nagylozs","balf","csapod","und","rojtokmuzsaj","brennberg","pusztacsalad","kutyahegy",...
%     "nyarliget","meszlen","fertoujlak","gorbehalom","tozeggyarmajor","ebergoc","csillahegy","jerevan","gloriette",...
%     "ohermes","ujhermes","linear_330"];

cases = ["ky1","ky2","ky3","ky4","ky5","ky6","ky7","ky8","ky9","ky10","ky11","ky12","ky13","ky14"];

idx = 1:14;
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
annotation('textarrow',[0.647,0.717],[0.71,0.71],'FontSize',20,'Linewidth',1.5)
text(1e-2,3e6,'Moderately','FontSize',14,'HorizontalAlignment','center');
text(1e-2,1.2e6,'vulnerable','FontSize',14,'HorizontalAlignment','center');
% text(1e-2,3e6,'Mérsékelten','FontSize',14,'HorizontalAlignment','center');
% text(1e-2,1.2e6,'sebezhetõ','FontSize',14,'HorizontalAlignment','center');
% highly vulnerable
plot([1e-2,1e-2],[1e-1,1e5],'--k','linewidth',1.5,'HandleVisibility','off')
annotation('textarrow',[0.732,0.802],[0.53,0.53],'FontSize',20,'Linewidth',1.5)
text(1e-1,3e4,'Highly','FontSize',14,'HorizontalAlignment','center');
text(1e-1,1.2e4,'vulnerable','FontSize',14,'HorizontalAlignment','center');
% text(1e-1,3e4,'Kritikusan','FontSize',14,'HorizontalAlignment','center');
% text(1e-1,1.2e4,'sebezhetõ','FontSize',14,'HorizontalAlignment','center');
for I=idx
%     data{1} = importdata(join(['Network Data/',cases(I),'/vulner_orig_rand.txt'],''));
%     data{1} = importdata(join(['Network Data/',cases(I),'/vulner_orig.txt'],''));
    data{1} = importdata(join(['Network Data/',cases(I),'/vulner_N.txt'],''));
%     data{1} = importdata(join(['Network Data/',cases(I),'/vulner_Nm1.txt'],''));
    set(gca,'XScale','log','YScale','log');
    for i=1:size(data,2)
        gammaBar = data{i};
        gammaBarNZ = gammaBar(gammaBar~=0);
        n = length(gammaBarNZ);
        if(n>100)
            r = round(log(n)/log(2)+1);
        else
            r = round(sqrt(n));
        end
        b = prctile(gammaBarNZ,(0:1/r:1)*100);
        x = (b(1:end-1)+b(2:end))/2;
    %     x = b(1:end-1);

        f = zeros(r,1);
        % frequency
        for j=1:r
            f(j) = length(gammaBarNZ(gammaBarNZ>=b(j) & gammaBarNZ<b(j+1)));
        end
        f = f/sum(f); % relative frequency

        y = f./diff(b');

        colr = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),(I-1)/size(idx,2),'spline')));
        colg = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),(I-1)/size(idx,2),'spline')));
        colb = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),(I-1)/size(idx,2),'spline')));

        if(mod(I,5) == 0)
            marker = 'x';
        elseif(mod(I,5) == 1)
            marker = 'o';
        elseif(mod(I,5) == 2)
            marker = '+';
        elseif(mod(I,5) == 3)
            marker = '*';
        elseif(mod(I,5) == 4)
            marker = 'p';
        end
        plot(x,y,marker,'color',[colr,colg,colb],'linewidth',1.5);

    end
end
ax=gca;
ax.FontSize = 18;
xlabel('Local vulnerability [-]','fontsize',24);
ylabel('Probability density [-]','fontsize',24);
% xlabel('Lokális sebezhetõség [-]','fontsize',14);
% ylabel('Valószínûségi sûrûség [-]','fontsize',14);
% legend(cases(idx));
xlim([1e-9,1e0])
ylim([1e-1,1e8])
% set(gca,'Position',[100,100,600,400])
% legend('Eredeti','N szabÃ¡ly','N-1 szabÃ¡ly');
% legend(cases,'location','southwest','FontSize',15);
ColumnLegend(2,cases(:),'location','southwest','FontSize',16);
rectangle('Position',[1.3e-9 7e-1 7.5e-7 5e3],'FaceColor',[1 1 1])
% set(gca,'FontSize',14);
% ColumnLegend(3,num2str(idx(:)));
% ColumnLegend(3,num2str(idx(:)),'location','southwest_sf');
% legend('Original','N rule','N-1 rule');
% rectangle('Position',[1.3e-9 1.9e-1 7.5e-7 5e3],'FaceColor',[1 1 1])
% txt1 = text(2e-8,3e3,'Ivóvízhálózat','FontSize',14,'HorizontalAlignment','center');
% txt2 = text(2e-8,1.2e3,'sorszám','FontSize',14,'HorizontalAlignment','center');
% txt1 = text(2e-8,3e3,'Water distribution','FontSize',14,'HorizontalAlignment','center');
% txt2 = text(2e-8,1.2e3,'network','FontSize',14,'HorizontalAlignment','center');
title('Failure rate: relative pipe length, N rule','FontSize',24);

% set(gca,'FontSize',16);
% title('Failure rate: pipe material');
saveas(gca,'Plots/Vulner_ky_N.fig','fig');
saveas(gca,'Plots/Vulner_ky_N.png','png');
saveas(gca,'Plots/Vulner_ky_N.eps','epsc');

% number of iso valves, segments
% iso_valves = zeros(length(idx_iso),3);
% no_segments = zeros(length(idx_iso),3);
% for i=idx_iso
%     data = importdata(join(['Network Data/',cases(i),'/number_of_valves.txt'],''));
%     iso_valves(i,:) = data(:,1);
%     no_segments(i,:) = data(:,2);
% end

% fig2 = figure('Position',[250 200 1500 500]);
% bar(iso_valves);
% title('Number of ISO valves');
% legend('Original','N rule','N-1 rule');
% 
% fig3 = figure('Position',[300 200 1500 500]);
% bar(no_segments);
% title('Number of segments');
% legend('Original','N rule','N-1 rule');
% 
% % network vulnerability
% 
% network_vulner=zeros(length(idx_net),3);
% for i=idx_net
%     data2{1} = importdata(join(['Network Data/',cases(i),'/vulner_orig.txt'],''));
%     data2{2} = importdata(join(['Network Data/',cases(i),'/vulner_N.txt'],''));
%     data2{3} = importdata(join(['Network Data/',cases(i),'/vulner_Nm1.txt'],''));
%     network_vulner(i,1) = sum(data2{1});
%     network_vulner(i,2) = sum(data2{2});
%     network_vulner(i,3) = sum(data2{3});
% end
% 
% fig4 = figure('Position',[350 200 1500 500]);
% bar(network_vulner);
% title('Network vulnerability');
% legend('Original','N rule','N-1 rule');







