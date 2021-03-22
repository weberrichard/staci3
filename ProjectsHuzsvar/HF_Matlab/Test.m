T = importdata('vulnerability.csv')
GG=T(:,5);
averagedemand=T(:,6);
figure;
GG=GG.';
averagedemand=averagedemand.';
scatter(averagedemand,GG);
 
