close all; 
clear;

addpath('../../Plot');

%blackBody, blackBodyExt, cividis, coolWarmBent, coolWarmSmooth, inferno, jet, kindlmann, kindlmannExt, magma, plasma, viridis
%discrete: lines, prism
% colorMapName = 'grayscale'; 
colorMapName = 'plasma';
colorMap = importdata(['../../Plot/ColorMaps/',colorMapName,'.col']);

n_lin = 100; % number of segments in linear wdn
n_grid = 10*10; % number of segments in grid wdn

m_pir = 4; % number of layers in pyramid wds
l_pir = 2; % number of new edges
n_pir = (1-l_pir^(m_pir+1))/(1-l_pir);

% vulnerability of artificial wdns
data{1} = (1-linspace(1,n_lin,n_lin)/n_lin)/n_lin;
data{2} = ones(1,n_grid)/n_grid^2;

v = zeros(n_pir,1);
v(1) = 1;
for i=1:m_pir
    li = 2^i;
    ui = 2^(i+1)-1;
    v(li:ui) = (v(i)-1/n_pir)/2;
end
%adding some rand
v = v.*(randn(n_pir,1)/100+1);
data{3} = v;


fig = figure('Position',[200 200 900 600]);
hold on; grid on;
set(gca,'XScale','log','YScale','log');
for i=3
    gammaBar = data{i};
    gammaBarNZ = gammaBar(gammaBar~=0);
    n = length(gammaBarNZ);
    if(n>100)
        r = round(log(n)/log(2)+1);
    else
        r = 2*round(sqrt(n));
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

    colr = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),(i-1)/size(data,2),'spline')));
    colg = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),(i-1)/size(data,2),'spline')));
    colb = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),(i-1)/size(data,2),'spline')));

    if(mod(i,5) == 0)
        marker = 'x';
    elseif(mod(i,5) == 1)
        marker = 'o';
    elseif(mod(i,5) == 2)
        marker = '+';
    elseif(mod(i,5) == 3)
        marker = '*';
    elseif(mod(i,5) == 4)
        marker = 'p';
    end
    plot(x,y,marker,'color',[colr,colg,colb],'linewidth',1.5);
end
xlabel('Local vulnerability','fontsize',14);
ylabel('Probability density','fontsize',14);
set(gca,'FontSize',14);

% saveas(gca,'Plots/GraphAbs.fig','fig');
% saveas(gca,'Plots/GraphAbs.png','png');
% saveas(gca,'Plots/GraphAbs.eps','epsc');




