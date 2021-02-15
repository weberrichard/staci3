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
for I=idx
    if(mod(I,2)==0)
        data = importdata(join(['Network Data/',cases(I),'/vulner_Nm1.txt'],''));
        col = 0.25;
    else
        data = importdata(join(['Network Data/',cases(I),'/vulner_orig.txt'],''));
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

xlabel('Local vulnerability','fontsize',14);
ylabel('Probability density','fontsize',14);
legend('Sopron Networks','ky*');
set(gca,'FontSize',14);

