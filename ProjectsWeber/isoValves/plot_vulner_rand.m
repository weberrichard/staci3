close all; 
clear;

addpath('../../Plot');

cases= ["villasor","ferto","sanchegy","buk","lovo","nagycenk","vashegy","varis","becsidomb","tomalom",...
    "szakov","kohegy","harka","pozsonyiut","sopronkovesd","dudlesz","ivan","agyagosszergeny","kofejto","simasag",...
    "acsad","csaford","nagylozs","balf","csapod","und","rojtokmuzsaj","brennberg","pusztacsalad","kutyahegy",...
    "nyarliget","meszlen","fertoujlak","gorbehalom","tozeggyarmajor","ebergoc","csillahegy","jerevan","gloriette",...
    "ohermes","ujhermes"];

% cases = ["ky1","ky2","ky3","ky4","ky5","ky6","ky7","ky8","ky9","ky10","ky11","ky12","ky13","ky14"];

idx = 1:27;

colors = [0.25, 0.75];

%blackBody, blackBodyExt, cividis, coolWarmBent, coolWarmSmooth, inferno, jet, kindlmann, kindlmannExt, magma, plasma, viridis
%discrete: lines, prism
% colorMapName = 'grayscale'; 
colorMapName = 'plasma';
colorMap = importdata(['../../Plot/ColorMaps/',colorMapName,'.col']);

fig = figure('Position',[200 200 900 600]);
hold on; grid on;
for I=idx
    data{1} = importdata(join(['Network Data/',cases(I),'/vulner_orig.txt'],''));
    data{2} = importdata(join(['Network Data/',cases(I),'/vulner_orig_rand.txt'],''));
%     data{1} = importdata(join(['Network Data/',cases(I),'/vulner_N.txt'],''));
%     data{2} = importdata(join(['Network Data/',cases(I),'/vulner_Nm1.txt'],''));
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

        colr = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),colors(i),'spline')));
        colg = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),colors(i),'spline')));
        colb = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),colors(i),'spline')));

        if(mod(i,2) == 0)
            marker = 'x';
        else
            marker = 'o';
        end
        plot(x,y,marker,'color',[colr,colg,colb],'linewidth',1.5);

    end
end

xlabel('Local vulnerability','fontsize',14);
ylabel('Probability density','fontsize',14);
legend('Original','Randomized');
set(gca,'FontSize',14);
% xlabel('Szegmens sebezhetőség [-]','interpreter','latex','fontsize',14);
% ylabel('Valószínűség sűrűség [-]','interpreter','latex','fontsize',14);
% xlim([1e-9,1e-1])
% ylim([1e-1,1e8])
% set(gca,'Position',[100,100,600,400])
% legend(cases(idx),'location','west');
% set(gca,'FontSize',14);
% ColumnLegend(3,num2str(idx(:)));
% ColumnLegend(3,num2str(idx(:)),'location','southwest_sf');
% legend('Original','N rule','N-1 rule');
% rectangle('Position',[1.3e-9 1.9e-1 3.5e-7 5e3],'FaceColor',[1 1 1])
% saveas(gca,'Plots/GraphAbs.fig','fig');
% saveas(gca,'Plots/GraphAbs.png','png');
% saveas(gca,'Plots/GraphAbs.eps','epsc');





