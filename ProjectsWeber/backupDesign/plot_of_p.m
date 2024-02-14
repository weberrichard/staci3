clear;

p = linspace(0,60,1000);
pref = 20;
p0 = 20;
alpha = log(p0)/pref;

of_p = exp(-alpha'*(p-pref));

figure;
plot(p,of_p);
% set(gca,'xscale','log','yscale','log');
