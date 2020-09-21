clear;

eps = logspace(-6,log10(0.05),11);
% eps = 0.01;
Re = logspace(log10(2000),7,20);
% Re = 3565.07;

f_imp = zeros(length(eps),length(Re));
% f_wiki = f_imp;
f_epa1 = f_imp;

f_wiki2 = f_imp;

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

% wiki 1
a = 2.51 ./ Re';
b = eps ./ 3.7;
f_wiki = 1 ./ ((2.*lambertw(log(10)./2./a.*10.^(b./2./a)))./log(10) - b./a).^2;

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


% wiki 2
Rew = [logspace(2,log10(2000),5),Re];
f_wiki2 = zeros(length(eps),length(Rew));
for i = 1:length(eps)
    for j = 1:length(Rew)
        a = 1 ./ (1+(Rew(j)'/2712).^8.4);
        b = 1 ./ (1+(Rew(j)'/150*eps(i)).^1.8);
        f_wiki2(i,j) = (64/Rew(j)')^a * (0.75*log(Rew(j)'/5.37))^(2*(a-1)*b) * (0.88*log(6.82/eps(i)))^(2*(a-1)*(1-b));
    end
end

% laminar
Re_lam = logspace(2,log10(2300),5);
f_lam = 64./Re_lam;

figure();
hold on; grid on;
set(gca,'xscale','log','yscale','log');
plot(1e4,1e-1,'k','linewidth',1.5);
plot(1e4,1e-1,'b','linewidth',1.5);
plot(1e4,1e-1,'r','linewidth',1.5);
plot(1e4,1e-1,'g','linewidth',1.5);
plot(1e4,1e-1,'m','linewidth',1.5);
plot(Re_lam,f_lam,'k','linewidth',1.5);
plot(Re,f_imp,'b','linewidth',1.5);
% plot(Re,f_wiki,'r','linewidth',1.5);
plot(Re1,f_epa1,'g','linewidth',1.5);
plot(Re2,f_epa2,'g','linewidth',1.5);
% plot(Rew,f_wiki2,'m','linewidth',1.5);
legend('laminar','imp','wiki 1', 'epa', 'wiki 2')







