clear;

% eps = logspace(-6,log10(0.05),11);
eps = [1e-6,5e-6,1e-5,5e-5,1e-4,2e-4,5e-4,1e-3,2e-3,5e-3,1e-2,2e-2,3e-2,4e-2,5e-2];
Re = logspace(log10(2000),7,200);

colorMapName = 'plasma'; 
colorMap = importdata(['../../Plot/ColorMaps/',colorMapName,'.txt']);

f_imp = zeros(length(eps),length(Re));
f_epa1 = f_imp;

% staci now
counter = 0;
for i = 1:length(eps)
    for j = 1:length(Re)
        f_new = 10;
        f_old = 0.03;
        while(abs((f_new - f_old)/f_old) > 1e-4)
            if(f_new ~= 10)
                f_old = f_new;
            end
            tmp = -2*log10(eps(i)/3.71 + 2.51/Re(j)/sqrt(f_old));
            f_new = 1 / tmp / tmp;
            counter = counter+1;
        end
        f_imp(i,j) = f_new;
    end
end
m = counter/length(eps)/length(Re);

% epanet
Re1 = Re(Re<4000);
Re1 = [Re1,4000];
Re2 = Re(Re>4000);
Re2 = [4000,Re2];

Y3 = -0.86859*log(eps/3.7+5.744/4000^.9);
Y2 = eps/3.7 + 5.74*Re1'.^0.9;
FA = Y3.^-2;
FB = FA.*(2-0.00514215./Y2./Y3);
R = Re1'/2000;
X4 = R.*(0.032-3*FA+0.5*FB);
X3 = -0.128 + 13*FA-2*FB;
X2 = 0.128-17*FA+2.5*FB;
X1 = 7*FA-FB;
f_epa1 = (X1 + R.*(X2+R.*(X3+X4)));

f_epa2 = zeros(length(eps),length(Re2));
for i = 1:length(eps)
    for j = 1:length(Re2)
        f_epa2(i,j) = 0.25 / ( log10(eps(i)/3.7 + 5.74/Re2(j)^0.9))^2;    
    end
end


% laminar
Re_lam = logspace(log10(200),log10(2000),5);
f_lam = 64./Re_lam;

fig = figure('position',[100 100 900 600]);
hold on; grid on;
plot(1e5,1e-1,'w')
plot(1e5,1e-1,'w')
set(gca,'xscale','log','yscale','log');
for i=length(eps):-1:1
    r = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),i/1.3/(length(eps)-1),'spline')));
    g = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),i/1.3/(length(eps)-1),'spline')));
    b = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),i/1.3/(length(eps)-1),'spline'))); 
%     plot(Re1,f_epa1(:,i),'color',[r,g,b],'linewidth',1.75);
    plot(Re2,f_epa2(i,:),'color',[r,g,b],'linewidth',1.75);
end

for i=length(eps):-1:1
    r = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),i/1.3/(length(eps)-1),'spline')));
    g = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),i/1.3/(length(eps)-1),'spline')));
    b = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),i/1.3/(length(eps)-1),'spline'))); 
    plot(Re1,f_epa1(:,i),'color',[r,g,b],'linewidth',1.75);
%     plot(Re2,f_epa2(i,:),'color',[r,g,b],'linewidth',1.75);
end

plot(Re_lam,f_lam,'color',[r,g,b],'linewidth',1.75);

plot([2000,2000],[4e-3,1],'--k','linewidth',1.1)
plot([4000,4000],[4e-3,1],'--k','linewidth',1.1)
xlabel('Reynolds number [-]','fontsize',14)
ylabel('Friction coefficient [-]','fontsize',14)
xlim([5e1,1e7])
ylim([7e-3,3e-1])

legend(["Relative";"roughness [-]";num2str(flip(eps)')],'location','southwest','fontsize',10)

saveas(fig,'darcy-weisbach.eps','epsc')
saveas(fig,'darcy-weisbach.fig','fig')




