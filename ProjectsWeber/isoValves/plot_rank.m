clear;
close all;

addpath('../../Plot');

cases= ["villasor","ferto","sanchegy","buk","lovo","nagycenk","vashegy","varis","becsidomb","tomalom",...
    "szakov","kohegy","harka","pozsonyiut","sopronkovesd","dudlesz","ivan","agyagosszergeny","kofejto","simasag",...
    "acsad","csaford","nagylozs","balf","csapod","und","rojtokmuzsaj","brennberg","pusztacsalad","kutyahegy",...
    "nyarliget","meszlen","fertoujlak","gorbehalom","tozeggyarmajor","ebergoc","csillahegy","jerevan","gloriette",...
    "ohermes","ujhermes"];

idx = 1:5;

colours = [0.2,0.5,0.8];

%blackBody, blackBodyExt, cividis, coolWarmBent, coolWarmSmooth, inferno, jet, kindlmann, kindlmannExt, magma, plasma, viridis
%discrete: lines, prism
% colorMapName = 'grayscale'; 
colorMapName = 'plasma';
colorMap = importdata(['../../Plot/ColorMaps/',colorMapName,'.col']);

figure()
hold on; grid on;

for i=idx
    data{1} = importdata(join(['Network Data/',cases(i),'/rank_orig.txt'],''));
    data{2} = importdata(join(['Network Data/',cases(i),'/rank_N.txt'],''));
    data{3} = importdata(join(['Network Data/',cases(i),'/rank_Nm1.txt'],''));

    max_rank = max([data{1}',data{2}',data{3}']);

    rank_orig = zeros(max_rank+1,1);
    rank_N = zeros(max_rank+1,1);
    rank_Nm1 = zeros(max_rank+1,1);

    for j=2:max_rank+1
        rank_orig(j) = sum(data{1}==j-1);
        rank_N(j) = sum(data{2}==j-1);
        rank_Nm1(j) = sum(data{3}==j-1);
    end

    rank_orig = rank_orig/length(data{1});
    rank_N = rank_N/length(data{2});
    rank_Nm1 = rank_Nm1/length(data{3});

    colr = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),colours(1),'spline')));
    colg = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),colours(1),'spline')));
    colb = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),colours(1),'spline')));
    plot(0:max_rank,rank_orig,'x','color',[colr,colg,colb],'linewidth',1.5,'markersize',10);

    colr = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),colours(2),'spline')));
    colg = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),colours(2),'spline')));
    colb = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),colours(2),'spline')));
    plot(0:max_rank,rank_N,'o','color',[colr,colg,colb],'linewidth',1.5,'markersize',10);

    colr = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),colours(3),'spline')));
    colg = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),colours(3),'spline')));
    colb = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),colours(3),'spline')));
    plot(0:max_rank,rank_Nm1,'+','color',[colr,colg,colb],'linewidth',1.5,'markersize',10);
end
legend('Original','N rule','N-1 rule');





