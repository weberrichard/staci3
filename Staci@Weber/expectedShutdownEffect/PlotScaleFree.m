close all; 
clear;

addpath('../../Plot');

cases= ["villasor","ferto","sanchegy","buk","lovo","nagycenk","vashegy","varis","becsidomb","tomalom",...
    "szakov","kohegy","harka","pozsonyiut","sopronkovesd","dudlesz","ivan","agyagosszergeny","kofejto","simasag",...
    "acsad","csaford","nagylozs","balf","csapod","und","rojtokmuzsaj","brennberg","pusztacsalad","kutyahegy",...
    "nyarliget","meszlen","fertoujlak","gorbehalom","tozeggyarmajor","ebergoc","csillahegy","jerevan","gloriette",...
    "ohermes","ujhermes"];

idx = 1:27;

%blackBody, blackBodyExt, cividis, coolWarmBent, coolWarmSmooth, inferno, jet, kindlmann, kindlmannExt, magma, plasma, viridis
%discrete: lines, prism
% colorMapName = 'grayscale'; 
colorMapName = 'plasma';

colorMap = importdata(['../../Plot/ColorMaps/',colorMapName,'.txt']);

% DEALING WITH ALFA*LENGTH STUFF
data = importdata('gammaBarRelLength.txt');
data = data(idx,:);
figure('Position',[100 100 900 600]);
hold on; grid on;
set(gca,'XScale','log','YScale','log');
everyX=[];
everyY=[];
everyData=[];
for i=1:size(data,1)
    gammaBar{i} = str2num(data{i});
    everyData = [everyData, gammaBar{i}];
    gammaBarNZ{i} = gammaBar{i}(gammaBar{i}~=0);
    n = length(gammaBarNZ{i});
    if(n>100)
        r = round(log(n)/log(2)+1);
    else
        r = round(sqrt(n));
    end
    b = prctile(gammaBarNZ{i},(0:1/r:1)*100);
    x = (b(1:end-1)+b(2:end))/2;
%     x = b(1:end-1);
    
    everyX = [everyX,x];
    
    f = zeros(r,1);
    % frequency
    for j=1:r
        f(j) = length(gammaBarNZ{i}(gammaBarNZ{i}>=b(j) & gammaBarNZ{i}<b(j+1)));
    end
    f = f/sum(f); % relative frequency
    
    y = f./diff(b');
    everyY = [everyY,y'];
    
    colr = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),(i-1)/size(data,1),'spline')));
    colg = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),(i-1)/size(data,1),'spline')));
    colb = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),(i-1)/size(data,1),'spline')));
    
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

    % drawing columns
%     cx = reshape([b;b], [], 1)';
%     cx = cx(2:end-1);
%     cy = reshape([y';y'], [], 1)';
%     plot(cx,cy,'color',[colr,colg,colb],'linewidth',1.1);
end
xlabel('Local vulnerability [-]','fontsize',14);
ylabel('Probability density [-]','fontsize',14);
% xlabel('Szegmens sebezhetõség [-]','interpreter','latex','fontsize',14);
% ylabel('Valószínûség sûrûség [-]','interpreter','latex','fontsize',14);
xlim([1e-9,1e-1])
ylim([1e-1,1e8])
% set(gca,'Position',[100,100,600,400])
% legend(cases(idx),'location','west');
set(gca,'FontSize',14);
ColumnLegend(3,num2str(idx(:)));
% ColumnLegend(3,cases(idx),'location','southwest_sf');
rectangle('Position',[3e-4 1.5e4 8e-2 6.5e7],'FaceColor',[1 1 1])
saveas(gca,'Plots/GammaDensity.fig','fig');
saveas(gca,'Plots/GammaDensity.png','png');
saveas(gca,'Plots/GammaDensity.eps','epsc');


% % CREATING THE OVERALL DENSITY FUNCTION
everyDataSort = sort(everyData(everyData~=0));
n = length(everyDataSort);
if(n>100)
    r = round(log(n)/log(2)+1);
else
    r = round(sqrt(n));
end
b = prctile(everyDataSort,(0:1/r:1)*100);
xo = (b(1:end-1)+b(2:end))/2;

f = zeros(r,1);
% frequency
for j=1:r
    f(j) = length(everyDataSort(everyDataSort>=b(j) & everyDataSort<b(j+1)));
end
f = f/sum(f); % relative frequency
yo = (f./diff(b'))';
% plot(xo,yo,'o','color',[0,0,0],'linewidth',1.5,'markersize',8,'markerfacecolor',[0,0,0]);


% lognormal FIT BRUTE FORCE YEEAAAAAHH
% dmin = min(everyDataSort);
dmin = 1e-7;
deltae = 1-length(everyData(everyData>dmin))/length(everyData);
xe = everyX(everyX>dmin);
meanY = mean(everyY(everyX>dmin));

muV = -logspace(0,2,100);
sigmaV = (logspace(0,2,100))';
CeV = (1-deltae) * sqrt(2./pi./sigmaV.^2) .* ( erf(-muV./sqrt(2)./sigmaV) - erf((log(dmin)-muV)./sqrt(2)./sigmaV) ).^-1;

ff = zeros(length(muV),length(sigmaV));
for i=1:length(muV)
    for j=1:length(sigmaV)
       ye = CeV(j,i) .* (1./xe .* exp(-(log(xe)-muV(i)).^2/2./sigmaV(j).^2));
       ff(i,j) = sum(log(ye./everyY(everyX>dmin)).^2);
    end
end
minMatrix = min(ff(:));
[row,col] = find(ff==minMatrix);
mu = muV(row);
sigma = sigmaV(col);

% BRUTE FORCE FOR R2log
% mu = -1000;
% sigma = 1;

% Maximum Likelihood Estimator for mu and sigma
% mu = sum(log(everyDataSort)) / length(everyDataSort);
% sigma = sum((log(everyDataSort) - mu).^2) / length(everyDataSort);

% Method of Moments estimator
% mu = -log(sum(everyDataSort.^2))/2 + 2*log(sum(everyDataSort)) - 3/2*log(length(everyDataSort));
% sigma = log(sum(everyDataSort.^2)) - 2*log(sum(everyDataSort)) + log(length(everyDataSort));

% R2 calc
Ce = (1-deltae) * sqrt(2/pi/sigma^2) * ( erf(-mu/sqrt(2)/sigma) - erf((log(dmin)-mu)/sqrt(2)/sigma) )^-1;
ye = Ce .* (1./xe .* exp(-(log(xe)-mu).^2/2./sigma.^2));
R2log = 1 - sum(log(ye./everyY(everyX>dmin)).^2) / sum(log(meanY./everyY(everyX>dmin)).^2);

figure('Position',[100 100 900 600]);
grid on; hold on;
set(gca,'XScale','log','YScale','log')
plot(everyX,everyY,'ko','markerfacecolor',[0,0,0],'markersize',3);
% xx = [min(everyX)/2,max(everyX)*2];
% yy = Ce .* (1./xx .* exp(-(log(xx)-mu).^2/2./sigma.^2));
[xes,si] = sort(xe);
plot(xes,ye(si),'r','linewidth',2.0);
plot(1e-7,1e3,'w','linewidth',1.5);
legend('Networks', 'Log-normal fit',['R^2_{log} = ',num2str(R2log)]);
xlabel('$\overline{\gamma}$','interpreter','latex','fontsize',14);
ylabel('$f(\overline{\gamma})$','interpreter','latex','fontsize',14);
% saveas(gca,'Plots/GammaDensityFit.fig','fig');
% saveas(gca,'Plots/GammaDensityFit.png','png');

% CHI SQUARE TEST 
p = 0.95;
chi2test = zeros(length(cases(idx)),1);
for i = 1:length(cases(idx))
    n = length(gammaBarNZ{i});
    if(n>100)
        r = round(log(n)/log(2)+1);
    else
        r = round(sqrt(n));
    end
    
    b = prctile(gammaBarNZ{i},(0:1/r:1)*100);
    
    f = zeros(r,1);
    ppi = zeros(r,1);
    % frequency
    for j=1:r
        f(j) = length(gammaBarNZ{i}(gammaBarNZ{i}>=b(j) & gammaBarNZ{i}<b(j+1)));
        pu = .5*(1+erf((log(b(j+1))-mu)/sigma));
        pl = .5*(1+erf((log(b(j))-mu)/sigma));
        ppi(j) = pu-pl;
    end
    
    chi2 = (f-n*ppi).^2./(n.*ppi);
    
    chi2obs = sum(chi2);
    chi2crit = chi2inv(p,r-1);
    if(chi2obs < chi2crit)
        chi2test(i) = 1;
    end
end
% legendText = [cases(idx),'Illesztett fv.'];
% ColumnLegend(2,legendText);
% % rectangle('Position',[4.6 0.43 3.35 0.55],'FaceColor',[1 1 1])
% saveas(gca,'RankDistFit.fig','fig');
% saveas(gca,'RankDistFit.png','png');

nOk = sum(chi2test);
nCases = length(cases(idx));
disp(['chi2 test: ', num2str(nOk),' of ', num2str(nCases), ' (',num2str(nOk/nCases*100),' %)']);


% for checking the territoriry under the function
Area = Ce*sqrt(pi/2)*sigma * ( erf((log(1)-mu)/sqrt(2)/sigma) - erf((log(dmin)-mu)/sqrt(2)/sigma) ) + deltae

% DRAWING THE CUMULATIVE DISTRIBUTION FUNCTION
% figure(17);cla;
% hold on; grid on;
% for i=1:size(data,1)
%     set(gca,'XScale','log','YScale','linear');
%     
%     nz = gammaBar{i}(gammaBar{i}~=0);
%     x = sort(nz);
%     x = reshape([x;x],[1,2*length(x)]);
%     x = [min(x)/10,x,max(x)*10];
%     
%     delta = sum(gammaBar{i}(gammaBar{i}==0))/length(gammaBar{i});
%     y = delta:1/length(nz):1;
%     y = reshape([y;y],[1,2*length(y)]);
%     
%     colr = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),(i-1)/size(data,1),'spline')));
%     colg = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),(i-1)/size(data,1),'spline')));
%     colb = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),(i-1)/size(data,1),'spline')));
%     plot(x,y,'color',[colr,colg,colb],'linewidth',1.5);
% end
% plot(logspace(-15,0,10),.5*(1+erf((log(logspace(-15,0,10))-mu)/sigma)));

% saveas(gca,'Plots/GammaDist.fig','fig');
% saveas(gca,'Plots/GammaDist.png','png');

% % USING CURVE FITTING STUFF
% % EXPONENT   f(x) = C*x^g
% [xData, yData] = prepareCurveData( everyX, everyY );
% % Set up fittype and options.
% x_min = min(everyX);
% x_max = max(everyX);
% ft = fittype(['(g+1)/(',num2str(1),'^(g+1)-',num2str(x_min),'^(g+1) * x^g)'],...
%     'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.StartPoint = -1.059;
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% 
% a_fit = fitresult.g;
% C_fit = (a_fit+1)/(1^(a_fit+1)-x_min^(a_fit+1));
% x_fit = min(everyX):(max(everyX)-min(everyX))/10:max(everyX);
% y_fit = C_fit*(x_fit.^a_fit);
% F_fit = sum((C_fit*(everyX.^a_fit) - everyY).^2);
% 
% plot(x_fit,y_fit,'Color',[0.0,0.7,0.0],'linewidth',1.5);
% plot(x_fit(1),y_fit(1),'w','linewidth',1.5);
% % coeff of determination
% R21 = 1 - F_fit / sum((mean(everyY)-everyY).^2);
% leg = legend([cases(idx),[num2str(C_fit),'*(x^{',num2str(a_fit),'})'],...
%     ['R^2 = ',num2str(R21)]]);
% 
% % probability
% x2 = x_min;
% x1 = 1;
% pi = C_fit * ( x_max^(a_fit+1) / (a_fit+1) - x_min^(a_fit+1) / (a_fit+1) ); 
% 
% % LOGARITHMIC FIT
% logx = log(everyX);
% logy = log(everyY);
% [xData, yData] = prepareCurveData( logx, logy );
% x_min = min(everyX);
% x_max = max(everyX);
% ft = fittype(['log((g+1)/(',num2str(x_max),'^(g+1)-',num2str(x_min),'^(g+1))) + g*x'],'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.StartPoint = -1.01;
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% 
% g_fit = fitresult.g;
% % g_fit = -0.9;
% x_fit = min(logx):(max(logx)-min(logx))/10:max(logx);
% C_fit = log((g_fit+1)/(1-x_min^(g_fit+1)));
% y_fit = C_fit + g_fit*x_fit;
% F_fit = sum((C_fit + g_fit*logx - logy).^2);
% plot(exp(x_fit),exp(y_fit),'Color',[1.0,0.0,0.0],'linewidth',1.5);
% plot(x_fit(1),y_fit(1),'w','linewidth',1.5);
% % coeff of determination
% R22 = 1 - F_fit / sum((mean(logy)-logy).^2);
% legend([cases(idx),[num2str(C_fit),'*(x^{-',num2str(g_fit),'})'], ['R^2 = ',num2str(R21)], ...
%     [num2str(C_fit),' + ',num2str(g_fit), '*x'],['R^2 = ',num2str(R22)]]);
% 
% 
% % CREATING THE OVERALL DENSITY FUNCTION
% 
% n = length(everyData);
% if(n>100)
%     r = round(log(n)/log(2)+1);
% else
%     r = round(sqrt(n));
% end
% b = prctile(everyData,(0:1/r:1)*100);
% x = (b(1:end-1)+b(2:end))/2;
% 
% f = zeros(r,1);
% % frequency
% for j=1:r
%     f(j) = length(everyData(everyData>=b(j) & everyData<b(j+1)));
% end
% f = f/sum(f); % relative frequency
% y = f./diff(b');
% 
% plot(x,y,'o','color',[0,0,0],'linewidth',1.5,'markersize',8,'markerfacecolor',[0,0,0]);

% drawing columns
% cx = reshape([b;b], [], 1)';
% cx = cx(2:end-1);
% cy = reshape([y';y'], [], 1)';
% plot(cx,cy,'color',[colr,colg,colb],'linewidth',1.1);





%%%%%%%%%%%%%%%%%%%%%%
% FITTING VILLAGE
%%%%%%%%%%%%%%%%%%%%%%

% USING CURVE FITTING STUFF
% EXPONENT   f(x) = a/x + a*log(x_min)/(1-x_min)
% [xData, yData] = prepareCurveData( everyX, everyY );
% % Set up fittype and options.
% x_min = min(everyX);
% ft = fittype(['a/x + a*log(',num2str(x_min),')/(1-',num2str(x_min),')'],...
%     'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.StartPoint = 0.0589;
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% 
% a_fit = fitresult.a;
% x_fit = min(everyX):(max(everyX)-min(everyX))/10:max(everyX);
% y_fit = a_fit./x_fit + a_fit*log(x_min)/(1-x_min);
% y_fit = y_fit + abs(min(y_fit));
% F_fit = sum((a_fit./everyX + a_fit*log(x_min)/(1-x_min) - everyY).^2);
% plot(x_fit,y_fit,'Color',[1.0,0.0,0.0],'linewidth',1.5);
% plot(x_fit(1),y_fit(1),'w','linewidth',1.5);
% % coeff of determination
% R22 = 1 - F_fit / sum((mean(everyY)-everyY).^2);
% legend([cases(idx),[num2str(C_fit),'*(x^{-',num2str(gamma_fit),'}-1)'], ['R^2 = ',num2str(R21)], ...
%     [num2str(a_fit),'/x+',num2str(a_fit*log(x_min)/(1-x_min))],['R^2 = ',num2str(R22)]]);


% % USING CURVE FITTING STUFF
% % EXPONENT   f(x) = a/x + a*log(x_min) + x_min
% [xData, yData] = prepareCurveData( everyX, everyY );
% % Set up fittype and options.
% x_min = min(everyX);
% ft = fittype(['a/x+a*log(',num2str(x_min),')+',num2str(x_min)],'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.StartPoint = 0.0589;
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% 
% a_fit = fitresult.a;
% % a_fit = 0.1;
% x_fit = min(everyX):(max(everyX)-min(everyX))/10:max(everyX);
% y_fit = a_fit./x_fit + a_fit*log(x_min) + x_min;
% F_fit = sum((a_fit./everyX + a_fit*log(x_min) + x_min - everyY).^2);
% plot(x_fit,y_fit,'Color',[0.7,0.0,0.0],'linewidth',1.5);
% plot(x_fit(1),y_fit(1),'w','linewidth',1.5);
% % coeff of determination
% R22 = 1 - F_fit / sum((mean(everyY)-everyY).^2);
% legend([cases(idx),[num2str(C_fit),'*(x^{-',num2str(gamma_fit),'}-1)'], ['R^2 = ',num2str(R21)], ...
%     [num2str(a_fit),'/x+',num2str(a_fit*log(x_min)+x_min)],['R^2 = ',num2str(R22)]]);


% % USING CURVE FITTING STUFF
% % EXPONENT   f(x) = D * x^g
% [xData, yData] = prepareCurveData( everyX, everyY );
% % Set up fittype and options.
% x_min = min(everyX);
% ft = fittype(['(g+1)/(1-',num2str(x_min),'^(g+1)) * x^g'],'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.StartPoint = -1.05;
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% 
% g_fit = fitresult.g;
% D_fit = (g_fit+1)/(1-x_min^(g_fit+1));
% x_fit = min(everyX):(max(everyX)-min(everyX))/10:max(everyX);
% y_fit = D_fit*x_fit.^g_fit;
% F_fit = sum((D_fit*everyX.^g_fit - everyY).^2);
% plot(x_fit,y_fit,'Color',[0.7,0.0,0.0],'linewidth',1.5);
% plot(x_fit(1),y_fit(1),'w','linewidth',1.5);
% % coeff of determination
% R22 = 1 - F_fit / sum((mean(everyY)-everyY).^2);
% legend([cases(idx),[num2str(C_fit),'*(x^{-',num2str(gamma_fit),'}-1)'], ['R^2 = ',num2str(R21)], ...
%      [num2str(D_fit),'*x^{-',num2str(g_fit),'}'],['R^2 = ', num2str(R22)]]);


% USING CURVE FITTING STUFF
% DOUBLE EXPONENT   f(x) = C * (x^-gamma - x^-delta)
% [xData, yData] = prepareCurveData( everyX, everyY );
% % Set up fittype and options.
% x_min = min(everyX);
% ft = fittype(['1/ ( 1/(g+1)*(1-',num2str(x_min),'^(g+1)) - 1/(d+1)*(1-'...
%     ,num2str(x_min),'^(d+1)) ) * ( x^g - x^d )'],'independent', 'x',...
%     'dependent', 'y' );
% 
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.StartPoint = [-1.1; -1.001];
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% 
% gf = fitresult.g;
% df = fitresult.d;
% Cf = 1 / ( 1/(gf+1) * (1-x_min^(gf+1)) - 1/(df+1) * (1-x_min^(df+1)) );
% x_fit = min(everyX):(max(everyX)-min(everyX))/10:max(everyX);
% y_fit = Cf*(x_fit.^gf-x_fit.^df);
% F_fit = sum((Cf*(everyX.^gf-everyX.^df) - everyY).^2);
% 
% plot(x_fit,y_fit,'Color',[0.7,0.0,0.0],'linewidth',1.5);
% plot(x_fit(1),y_fit(1),'w','linewidth',1.5);
% % coeff of determination
% R22 = 1 - F_fit / sum((mean(everyY)-everyY).^2);
% legend([cases(idx),[num2str(C_fit),'*(x^{-',num2str(gamma_fit),'}-1)'],...
%     ['R^2 = ',num2str(R21)],[num2str(Cf),...
%     '*(x^{',num2str(gf),'}-x^{',num2str(df),'})'],['R^2 = ',num2str(R22)]]);


% USING CURVE FITTING STUFF
% % EXPONENT   f(x) = (gamma-1) * x_min^(gamma-1) * x^-gamma
% [xData, yData] = prepareCurveData( everyX, everyY );
% % Set up fittype and options.
% x_min = min(everyX);
% ft = fittype(['(a-1)*',num2str(x_min),'^(a-1)*x^-a'], 'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.StartPoint = 1.059;
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% 
% gamma_fit = fitresult.a; % (1.059, 1.06)
% % gamma_fit = 1.001;
% C_fit = (gamma_fit-1)*min(everyX)^(gamma_fit-1);
% x_fit = min(everyX):(max(everyX)-min(everyX))/10:max(everyX);
% y_fit = C_fit*x_fit.^-gamma_fit;
% F_fit = sum((C_fit*everyX.^-gamma_fit - everyY).^2);
% 
% plot(x_fit,y_fit,'Color',[0.0,0.7,0.0],'linewidth',1.5);
% plot(x_fit(1),y_fit(1),'w','linewidth',1.5);
% % coeff of determination
% R21 = 1 - F_fit / sum((mean(everyY)-everyY).^2);
% legend([cases(idx),[num2str(C_fit),'*x^{-',num2str(gamma_fit),'}'],['R^2 = ',num2str(R21)]]);


% power law with cutoff
% dmin = 1e-7;
% deltae = 1-length(everyData(everyData>dmin))/length(everyData);
% alfae = 0.75;
% lambdae = 7e1;
% % Ce = 0.16;
% Ce = lambdae^(1-alfae)*(1-deltae)/(gammainc(lambdae*dmin,1-alfae,'upper')-gammainc(lambdae,1-alfae,'upper'));
% xe = everyX;
% ye = Ce * xe.^-alfae .* exp(-lambdae*xe);
% plot(xe,ye,'xb');
% 
% pi = Ce*(lambdae^(alfae-1)*(gammainc(lambdae*dmin,1-alfae,'upper')-gammainc(lambdae,1-alfae,'upper'))) + deltae

% ESTIMATION OF THE EXPONENT ACCORDING TO AARON CLAUSET
% dataMin = 1e-8;
% everyDataSort = sort(everyData(everyData>dataMin));
% % dataMin = everyDataSort(1);
% alfa = 1 + length(everyDataSort) * sum(log(everyDataSort/dataMin))^-1
% delta = 1-length(everyDataSort)/length(everyData);
% Cf = (1-delta)*(1-alfa)/(1-dataMin^(1-alfa));
% ff = plot(everyX,Cf*everyX.^-alfa,'k');

% dataMin = 1e-8;
% A = dataMin;
% logB = sum(log(everyDataSort));
% n = length(everyDataSort);
% alfa2 = (n*log(A)-logB-n)/(n*log(A)-logB);
% 
% x = 1.001:0.001:2;
% C = 1-n/length(everyData);
% out = n*(1./(1-x) + A.^(1-x).*log(A)./(1-A.^(1-x))) + logB;
% out2 = n .* ( -(1-C)./(1-A.^(1-x)) - A.^(1-x).*(1-C).*(1-x).*log(A)...
%     ./ (1-A.^(1-x)).^2 ) - ((1-C).*(1-x)./(1-A.^(1-x))).*logB;


% lognormal
% dmin = 1e-8;
% deltae = 1-length(everyData(everyData>dmin))/length(everyData);
% mu = -100;
% sigma = 100;
% xe = everyX;
% Ce = (1-deltae) * sqrt(2/pi/sigma^2) * ( erf(-mu/sqrt(2)/sigma) - erf((log(dmin)-mu)/sqrt(2)/sigma) )^-1;
% % Ce = 0.1;
% ye = Ce * (1./xe .* exp((-log(xe)-mu).^2/2/sigma^2));
% plot(xe,ye,'bs');
