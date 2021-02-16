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
for I=idx
%     data{1} = importdata(join(['Network Data/',cases(I),'/vulner_orig.txt'],''));
    data{1} = importdata(join(['Network Data/',cases(I),'/vulner_N.txt'],''));
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

xlabel('Local vulnerability','fontsize',14);
ylabel('Probability density','fontsize',14);
legend(cases(idx));
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







