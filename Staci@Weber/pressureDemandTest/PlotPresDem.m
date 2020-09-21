clear;

addpath('../../Plot');
% blackBody, blackBodyExt, cividis, coolWarmBent, coolWarmSmooth, inferno, jet, kindlmann, 
% kindlmannExt, magma, plasma, viridis
colorMapName = 'plasma'; 

colorMap = importdata(['../../Plot/ColorMaps/',colorMapName,'.txt']);

pmin = 15;
pdes = 25;
m = [1.5,2.0,2.5,3.,3.5];
dreq = 1;

p = pmin:0.01:pdes;
d = zeros(size(p));
dd = zeros(size(p));

fig1 = figure;
hold on; grid on;
for j=1:length(m)
    for i=1:length(p)
        if(p(i)<=pmin)
              d(i) = 0;
              dd(i) = 0;
        elseif(p(i)>pmin && p(i)<pdes)
              d(i) =  dreq*((p(i)-pmin)/(pdes-pmin))^(1/m(j));
              dd(i) = dreq*(1/m(j))*((p(i)-pmin)/(pdes-pmin))^(1/m(j)-1)/(pdes-pmin);
        else
              d(i) =  dreq;
              dd(i) = 0;
        end
    end
    colr = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),(j-1)/size(m,2),'spline')));
    colg = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),(j-1)/size(m,2),'spline')));
    colb = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),(j-1)/size(m,2),'spline')));
    
    if(m(j) == 2.5)
        lw = 2.5;
        colr_width = colr;
        colg_width = colg;
        colb_width = colb;
    else
        lw = 1.1;
    end
    plot(p,d,'linewidth',lw,'color',[colr,colg,colb]);
end
plot([10,15],[0,0],'linewidth',lw,'color',[colr_width,colg_width,colb_width]);
plot([25,30],[1,1],'linewidth',lw,'color',[colr_width,colg_width,colb_width]);
xlabel('Pressure [m]' ,'fontsize',12);
ylabel('Relative consumption (c/d) [-]' ,'fontsize',12);
xticks([pmin,pdes])
xticklabels({'p_{min}','p_{des}'})
lgstr = join([repmat("m = ",length(m),1),num2str(m')]);
legend(lgstr);
saveas(fig1,'PresDem.png','png');
saveas(fig1,'PresDem.fig','fig');
saveas(fig1,'PresDem.eps','epsc');

