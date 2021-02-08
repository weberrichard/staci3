clear;
close all;

cases= ["villasor","ferto","sanchegy","buk","lovo","nagycenk","vashegy","varis","becsidomb","tomalom",...
    "szakov","kohegy","harka","pozsonyiut","sopronkovesd","dudlesz","ivan","agyagosszergeny","kofejto","simasag",...
    "acsad","csaford","nagylozs","balf","csapod","und","rojtokmuzsaj","brennberg","pusztacsalad","kutyahegy",...
    "nyarliget","meszlen","fertoujlak","gorbehalom","tozeggyarmajor","ebergoc","csillahegy","jerevan","gloriette",...
    "alomhegy","ohermes","ujhermes"];

%   i | n_pres | n_pool |  n_seg |  n_iso | nodes |   n1/n3+ |   L1/L3+ |   V1/V3+ |     alfa |    alfa2 |     beta |    gamma |   l [km] |  c [cmh] |      apl |     diam |  cluster |    artic
% --------------------------------------------------------------------------------------------------------------------------------------------------------------------
M1=[
  1 ,      1 ,      2 ,    717 ,    882 ,  7188 ,   0.2934 ,   0.2001 ,   0.0889 ,   0.1131 ,   0.1131 , 3.2356e-04 , 4.5798e-01 ,  97.2850 , 539.1875 ,  16.7692 ,       42 ,   0.0156 ,      105 ;
  2 ,      5 ,      1 ,    461 ,    533 ,  6908 ,   0.5587 ,   0.5360 ,   0.2964 ,   0.1458 ,   0.1458 , 1.9487e-04 , 1.3628e-02 , 111.9518 ,  51.9506 ,  26.6182 ,       73 ,   0.0097 ,      130 ;
  3 ,      1 ,      0 ,    222 ,    268 ,  2940 ,   0.5638 ,   0.2309 ,   0.0837 ,   0.0918 ,   0.0918 , 1.9277e-04 , 1.4598e-02 ,  36.2403 ,  50.4736 ,  13.1093 ,       34 ,   0.0183 ,       56 ;
  4 ,      2 ,      2 ,    143 ,    178 ,  2912 ,   0.5738 ,   0.2948 ,   0.1647 ,   0.1418 ,   0.1418 , 2.3471e-03 , 4.3418e-01 ,  49.9463 , 310.9869 ,  10.7592 ,       27 ,   0.0347 ,       49 ;
  5 ,      5 ,      2 ,    173 ,    193 ,  2459 ,   0.6923 ,   0.4625 ,   0.2330 ,   0.3572 ,   0.3572 , 3.0344e-04 , 3.5940e-02 ,  45.2231 ,  30.5721 ,  21.2484 ,       60 ,   0.0099 ,       71 ;
  6 ,      2 ,      2 ,    142 ,    159 ,  2058 ,   0.7000 ,   0.6368 ,   0.4181 ,   0.8265 ,   0.8265 , 9.2214e-05 , 3.2658e-02 ,  33.1569 ,  10.3756 ,  17.3137 ,       43 ,   0.0000 ,       61 ;
  7 ,      1 ,      4 ,    179 ,    201 ,  1991 ,   0.6190 ,   0.3359 ,   0.1424 ,   0.6276 ,   0.6276 , 4.3397e-04 , 6.5391e-02 ,  31.9167 ,  36.0580 ,  15.1267 ,       38 ,   0.0181 ,       63 ;
  8 ,      1 ,      1 ,    226 ,    256 ,  1736 ,   0.5172 ,   0.2979 ,   0.0905 ,   0.3060 ,   0.3060 , 2.1618e-03 , 9.9482e-02 ,  27.1178 ,  78.2735 ,  14.8443 ,       43 ,   0.0000 ,       79 ;
  9 ,      1 ,      1 ,    195 ,    231 ,  1590 ,   0.3412 ,   0.3751 ,   0.1345 ,   0.2216 ,   0.2216 , 8.3940e-05 , 1.5546e-02 ,  22.6747 ,  19.4385 ,  10.4627 ,       24 ,   0.0000 ,       40 ;
 10 ,      1 ,      0 ,     65 ,     99 ,  1294 ,   0.2571 ,   0.0926 ,   0.0644 ,   0.1312 ,   0.1312 , 4.7355e-05 , 2.5074e-02 ,  17.6230 ,   6.1531 ,   5.2861 ,       11 ,   0.0270 ,        9 ;
 11 ,      1 ,      1 ,     85 ,     92 ,  1083 ,   0.7778 ,   0.6029 ,   0.3823 ,   0.7267 ,   0.7267 , 1.9506e-04 , 8.1120e-02 ,  19.4874 ,   8.6565 ,  13.0067 ,       34 ,   0.0000 ,       36 ;
 12 ,      1 ,      0 ,     57 ,     61 ,   976 ,   0.8571 ,   0.5389 ,   0.4497 ,   0.3461 ,   0.3461 , 1.8308e-03 , 1.3859e-01 ,  13.0838 ,  20.1315 ,   6.5802 ,       15 ,   0.0000 ,       26 ;
 13 ,      1 ,      0 ,     45 ,     54 ,   808 ,   0.5556 ,   0.4524 ,   0.3507 ,   0.1532 ,   0.1532 , 2.5927e-04 , 4.9883e-02 ,  10.7363 ,  10.2996 ,   5.2283 ,       13 ,   0.0309 ,       11 ;
 14 ,      2 ,      1 ,     74 ,     79 ,   765 ,   0.9600 ,   0.8259 ,   0.2698 ,   0.1466 ,   0.1466 , 7.0947e-04 , 4.8501e-03 ,  22.4582 , 526.6205 ,   8.5128 ,       22 ,   0.0000 ,       31 ;
 15 ,      1 ,      1 ,     64 ,     72 ,   696 ,   0.5238 ,   0.8457 ,   0.5944 ,   0.4235 ,   0.4235 , 6.9715e-05 , 4.0027e-02 ,  10.0819 ,   6.5767 ,   8.1290 ,       19 ,   0.0571 ,       15 ;
 16 ,      1 ,      0 ,     47 ,     53 ,   615 ,   1.0667 ,   0.6459 ,   0.3474 ,   0.1352 ,   0.1352 , 1.4544e-04 , 1.0166e-01 ,  17.9590 ,   5.1270 ,   4.9315 ,       11 ,   0.0288 ,       17 ;
 17 ,      1 ,      1 ,     44 ,     48 ,   596 ,   0.6875 ,   3.0793 ,   1.4043 ,   0.4201 ,   0.4201 , 1.2603e-04 , 9.6051e-02 ,   7.9404 ,   4.7235 ,   6.4461 ,       15 ,   0.0405 ,       15 ;
 18 ,      0 ,      1 ,     63 ,     73 ,   595 ,   0.4167 ,   1.0457 ,   0.8303 ,   0.3431 ,   0.3431 , 9.7343e-06 , 8.0142e-02 ,   8.8966 ,   0.4050 ,   6.7921 ,       17 ,   0.0000 ,        9 ;
 19 ,      1 ,      0 ,     45 ,     52 ,   549 ,   0.5294 ,   0.3131 ,   0.4515 ,   0.2682 ,   0.2682 , 2.7542e-04 , 1.8220e-01 ,   7.2223 ,   5.3955 ,   5.0899 ,       11 ,   0.0000 ,       11 ;
 20 ,      1 ,      1 ,     41 ,     45 ,   570 ,   0.8000 ,   0.4983 ,   0.4448 ,   0.4583 ,   0.4583 , 2.4103e-04 , 2.2054e-01 ,  10.0889 ,   3.9345 ,   6.7610 ,       16 ,   0.0405 ,       17 ;
 21 ,      1 ,      0 ,     27 ,     29 ,   488 ,   0.7500 ,   0.4571 ,   0.3040 ,   0.4655 ,   0.4655 , 1.7357e-03 , 4.5573e-01 ,   8.9835 ,  14.4520 ,   5.2934 ,       13 ,   0.0000 ,       12 ;
 22 ,      1 ,      2 ,     38 ,     39 ,   459 ,   1.0833 ,   1.0065 ,   0.5621 ,  38.2099 ,  38.2099 , 2.0815e-05 , 1.2184e-01 ,   9.9325 ,   0.6150 ,   6.8890 ,       18 ,   0.0000 ,       24 ;
 23 ,      2 ,      1 ,     27 ,     32 ,   437 ,   0.7273 ,   0.5028 ,   0.4297 ,   2.9292 ,   2.9292 , 2.5804e-05 , 5.9468e-02 ,   6.8783 ,   1.6920 ,   3.6809 ,        7 ,   0.0508 ,        8 ;
 24 ,      0 ,      1 ,     12 ,     12 ,   379 ,   1.0000 ,   1.3322 ,   0.7971 ,   0.3722 ,   0.3722 , 7.3161e-04 , 1.1791e-01 ,   5.0957 ,   8.8990 ,   2.8788 ,        6 ,   0.0000 ,        5 ;
 25 ,      1 ,      1 ,     24 ,     24 ,   378 ,   1.1250 ,   1.0184 ,   1.0304 ,   0.5093 ,   0.5093 , 1.6897e-04 , 2.1535e-01 ,   6.3901 ,   2.8275 ,   4.6486 ,       11 ,   0.0000 ,       12 ;
 26 ,      1 ,      1 ,     29 ,     31 ,   316 ,   0.7500 ,   2.2615 ,   2.3633 ,   4.0408 ,   4.0408 , 3.2704e-05 , 7.6207e-02 ,   7.6735 ,   2.0869 ,   5.4335 ,       13 ,   0.0000 ,       12 ;
 27 ,      1 ,      1 ,     28 ,     28 ,   301 ,   1.0000 ,   7.1298 ,   4.5626 ,   1.0074 ,   1.0074 , 1.0106e-04 , 1.9707e-01 ,   4.5924 ,   1.8780 ,   5.3122 ,       13 ,   0.0000 ,       14 ;
 28 ,      0 ,      1 ,     26 ,     28 ,   283 ,   1.7143 ,   1.4811 ,   1.2605 ,   0.1680 ,   0.1680 , 1.4212e-04 , 5.7722e-01 ,   6.9343 ,   3.8300 ,   3.8215 ,        8 ,   0.0462 ,        7 ;
 29 ,      2 ,      1 ,     23 ,     25 ,   233 ,   0.8000 ,   0.4451 ,   0.2990 ,   0.4575 ,   0.4575 , 3.5425e-05 , 1.0984e-01 ,   4.2050 ,   1.1610 ,   3.7352 ,        8 ,   0.0000 ,        8 ;
 30 ,      1 ,      0 ,      8 ,      8 ,   161 ,   2.0000 ,   0.8046 ,   0.8049 ,   0.3684 ,   0.3684 , 1.4492e-04 , 1.5537e-01 ,   2.1014 ,   1.4385 ,   2.1071 ,        4 ,   0.2000 ,        3 ;
 31 ,      1 ,      0 ,      8 ,      9 ,   162 ,   0.5000 ,   0.0911 ,   0.0911 ,   0.3815 ,   0.3815 , 1.1163e-04 , 2.7942e-01 ,   2.4904 ,   1.3749 ,   1.9286 ,        3 ,   0.0000 ,        1 ;
 32 ,      1 ,      1 ,      8 ,      8 ,   147 ,   1.0000 ,   4.3745 ,   4.3745 ,   7.2737 ,   7.2737 , 3.1317e-05 , 1.3791e-01 ,   3.8759 ,   0.8175 ,   2.3571 ,        5 ,   0.0000 ,        3 ;
 33 ,      1 ,      0 ,     11 ,     10 ,   133 ,   2.3333 ,   3.1756 ,   2.4704 ,   0.2708 ,   0.2708 , 5.1336e-05 , 2.4952e-01 ,   2.8356 ,   0.6825 ,   2.5818 ,        5 ,   0.0000 ,        4 ;
 34 ,      0 ,      1 ,     22 ,     24 ,   151 ,   1.0000 ,   0.6836 ,   0.6627 ,   0.2085 ,   0.2085 , 2.2341e-04 , 7.0112e-01 ,   5.6910 ,   4.1400 ,   3.4935 ,        8 ,   0.0000 ,        8 ;
 35 ,      1 ,      1 ,      7 ,      6 ,   121 ,   2.0000 ,   2.3395 ,   2.8032 ,   4.6151 ,   4.6151 , 7.3514e-05 , 4.5007e-01 ,   1.8588 ,   0.5880 ,   2.1905 ,        4 ,   0.0000 ,        3 ;
 36 ,      1 ,      0 ,      7 ,      6 ,   110 ,   2.0000 ,   4.5196 ,   4.3933 ,   0.3835 ,   0.3835 , 7.0946e-05 , 3.6992e-01 ,   1.4528 ,   0.6660 ,   2.1905 ,        4 ,   0.0000 ,        3 ;
 37 ,      1 ,      0 ,     13 ,     13 ,    93 ,   1.6667 ,   0.9330 ,   0.8926 ,   0.3943 ,   0.3943 , 7.6157e-05 , 3.4893e-01 ,   1.6092 ,   0.7853 ,   3.0256 ,        6 ,   0.0000 ,        7 ;
 38 ,      1 ,      0 ,     18 ,     20 ,    95 ,   0.5000 ,   0.1478 ,   0.1229 ,   0.3922 ,   0.3922 , 4.8556e-04 , 2.4937e-01 ,   1.9119 ,   5.9599 ,   3.3072 ,        7 ,   0.0000 ,        3 ;
 39 ,      1 ,      0 ,      6 ,      5 ,    61 ,   4.0000 ,  18.2604 ,  18.2833 ,   0.7127 ,   0.7127 , 1.2140e-04 , 7.0921e-01 ,   1.0122 ,   0.6132 ,   1.8667 ,        3 ,   0.0000 ,        2 ;
 40 ,      1 ,      0 ,      1 ,      0 ,    28 ,     -nan ,     -nan ,     -nan ,   1.0000 ,   1.0000 , 5.0987e-05 , 0.0000e+00 ,   0.5113 ,   0.1836 ,     -nan ,        0 ,   0.0000 ,        0 ;
 41 ,      1 ,      2 ,      6 ,      5 ,    32 ,   3.0000 ,  69.8286 ,  70.9304 , 203.7408 , 203.7408 , 3.0104e-06 , 9.5204e-01 ,   1.5090 ,   0.3465 ,   2.0667 ,        4 ,   0.0000 ,        3 ;
 42 ,      1 ,      0 ,      5 ,      4 ,    28 ,   3.0000 ,   2.6270 ,   2.6270 ,   0.8539 ,   0.8539 , 7.9343e-05 , 8.4941e-01 ,   0.3867 ,   0.3345 ,   1.8000 ,        3 ,   0.0000 ,        2 ;];

% SET INDEX
% idx = M1(:,1); %ALL
idx = [1:27]; %BEST
% idx = [1:39,41,42]; %LARGE

M1 = M1(idx,:);
cases = cases(idx);

nSeg = M1(:,4);
nISO = M1(:,5);
nNodes = M1(:,6);
n1n3p = M1(:,7);
L1L3p = M1(:,8);
V1V3p = M1(:,9);
gamma = M1(:,13);
lPipes = M1(:,14);
cons = M1(:,15);
apl = M1(:,16);
diam = M1(:,17);
cc = M1(:,18);
ac = M1(:,19);

% LOADING THE SEGMENT EDGE VECTOR
sg = zeros(length(idx),1); % spectral gap
ld = zeros(length(idx),1); % link density
and = zeros(length(idx),1); % average node degree
mc = zeros(length(idx),1); % meshdness coeff
for i=1:size(idx,2)
%     i=42;
    ev = importdata(join(["Network Data/",cases(i),"/EdgeVector.txt"],''));
    ev = ev+1;
    G = graph(ev(:,1),ev(:,2));
    Adj = adjacency(G);
    eigval = eig(Adj);
    eigval = sort(eigval);
    sg(i) = eigval(end)-eigval(end-1); 
%     lap = laplacian(G);
%     ac = sort(eig(lap));
%     ac = ac(end-1);

end
ld = 2*nISO./nSeg./(nSeg-1);
and = 2.*nISO./nSeg;
mc = (nISO + 1 - nSeg)./(2*nSeg-5);


disp(' ');
disp('Correlations');
disp(' gamma');
disp([' n1n3p : ', num2str(corr(n1n3p,gamma)), '    ', num2str(corr(n1n3p,gamma,'Type','Spearman'))]);
disp([' L1L3p : ', num2str(corr(L1L3p,gamma)), '    ', num2str(corr(L1L3p,gamma,'Type','Spearman'))]);
disp([' V1V3p : ', num2str(corr(V1V3p,gamma)), '    ', num2str(corr(V1V3p,gamma,'Type','Spearman'))]);
disp([' apl   : ', num2str(corr(apl,gamma)), '    ', num2str(corr(apl,gamma,'Type','Spearman'))]);
disp([' diam  : ', num2str(corr(diam,gamma)), '    ', num2str(corr(diam,gamma,'Type','Spearman'))]);
disp([' cc    : ', num2str(corr(cc,gamma)), '    ', num2str(corr(cc,gamma,'Type','Spearman'))]);
disp([' sg    : ', num2str(corr(sg,gamma)), '    ', num2str(corr(sg,gamma,'Type','Spearman'))]);
disp([' ld    : ', num2str(corr(ld,gamma)), '    ', num2str(corr(ld,gamma,'Type','Spearman'))]);
disp([' and   : ', num2str(corr(and,gamma)), '    ', num2str(corr(and,gamma,'Type','Spearman'))]);
disp([' mc    : ', num2str(corr(mc,gamma)), '    ', num2str(corr(mc,gamma,'Type','Spearman'))]);
disp([' ac    : ', num2str(corr(ac,gamma)), '    ', num2str(corr(ac,gamma,'Type','Spearman'))]);
disp([' ac/n  : ', num2str(corr(ac./nSeg,gamma)), '    ', num2str(corr(ac./nSeg,gamma,'Type','Spearman'))]);
disp(' ');

figure();grid on; hold on;
plot(gamma,n1n3p,'x');
plot(gamma,L1L3p,'o');
plot(gamma,V1V3p,'p');
plot(gamma,apl,'d');
plot(gamma,diam,'*');
plot(gamma,cc,'x');
plot(gamma,sg,'o');
plot(gamma,ld,'p');
plot(gamma,and,'d');
plot(gamma,mc,'*');
plot(gamma,ac,'x');
plot(gamma,ac./nSeg,'o');
set(gca,'XScale','log','YScale','log')

% BOX PLOT
figure(2156);cla;
set(gcf,'units','points','position',[100,100,600,200])
boxplot(gamma,'orientation','horizontal')
xlim([0,0.5])
a = gca;
a.Children.Children(1).MarkerEdgeColor = [0.87 0.4 0.38]; %outlier
a.Children.Children(1).LineWidth =1.5;
a.Children.Children(2).Color = [0.87 0.4 0.38]; %median
a.Children.Children(2).LineWidth =1.5;
a.Children.Children(3).Color = [0,0,0]; %box
a.Children.Children(3).LineWidth =1.5;
a.Children.Children(4).Color = [0.468 0.0 0.662]; %bottom
a.Children.Children(4).LineWidth =1.5;
a.Children.Children(5).Color = [0.468 0.0 0.662]; %top
a.Children.Children(5).LineWidth =1.5;
a.Children.Children(6).Color = [0.468 0.0 0.662]; %bottom
a.Children.Children(6).LineWidth =1.5;
a.Children.Children(7).Color = [0.468 0.0 0.662]; %top
a.Children.Children(7).LineWidth =1.5;

% a.Children.Children(1).MarkerEdgeColor = [0.2 0.2 0.2]; %outlier
% a.Children.Children(1).LineWidth =1.5;
% a.Children.Children(2).Color = [0.2 0.2 0.2]; %median
% a.Children.Children(2).LineWidth =1.5;
% a.Children.Children(3).Color = [0,0,0]; %box
% a.Children.Children(3).LineWidth =1.5;
% a.Children.Children(4).Color = [0.4 0.4 0.4]; %bottom
% a.Children.Children(4).LineWidth =1.5;
% a.Children.Children(5).Color = [0.4 0.4 0.4]; %top
% a.Children.Children(5).LineWidth =1.5;
% a.Children.Children(6).Color = [0.4 0.4 0.4]; %bottom
% a.Children.Children(6).LineWidth =1.5;
% a.Children.Children(7).Color = [0.4 0.4 0.4]; %top
% a.Children.Children(7).LineWidth =1.5;

set(gca,'ytick',[])
ylim([0.8,1.2])
xlabel('Network vulnerability [-]')
% xlabel('Hálózati sebezhetõség [-]')
set(gca,'FontSize',14);
saveas(gca,'Plots/GammaBoxPlot.fig','fig');
saveas(gca,'Plots/GammaBoxPlot.png','png');
saveas(gca,'Plots/GammaBoxPlot.eps','epsc');