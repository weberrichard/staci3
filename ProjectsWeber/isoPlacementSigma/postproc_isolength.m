clear; close all;

cases = ["tomalom","tomalom_sigma_opt","szakov","szakov_sigma_opt"];

for i=[1,3]
   l_orig = importdata("Network Data/" + cases(i) + "/isoLength.txt");
   l_opt = importdata("Network Data/" + cases(i+1) + "/isoLength.txt");
   
   figure();
   histogram(l_orig,10);
   hold on;
   histogram(l_opt,10);
   legend('original','optimized');
   title(cases(i));
end
