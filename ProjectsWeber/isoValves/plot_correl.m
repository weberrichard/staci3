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
colorMap = importdata(['../../Plot/ColorMaps/',colorMapName,'.col']);

case_type = 'orig';

a = 1;
b = 1;
Lambda_f1 = zeros(length(idx),1);
Lambda_f2 = zeros(length(idx),1);
Lambda = zeros(length(idx),1);
Gamma = zeros(length(idx),1);

for i=idx
    data = importdata(join(['Network Data/',cases(i),'/Lambda_',case_type,'.txt'],''));
    Lambda_f1(i) = data(1);
    Lambda_f2(i) = data(2);
    Lambda(i) = a*Lambda_f1(i) + b*Lambda_f2(i);
    Gamma(i) = importdata(join(['Network Data/',cases(i),'/network_vulner_',case_type,'.txt'],''));
end

figure();
plot(Lambda,Gamma,'x','linewidth',1.5,'markersize',8);
xlabel('Lambda [-]');
ylabel('Gamma [-]');
saveas(gca,'Plots/Lambda_Gamma_correl.png','png');
corr(Lambda,Gamma,'Type','Pearson')
corr(Lambda,Gamma,'Type','Spearman')
