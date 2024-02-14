clear;

network_name = 'balf_mat_year_simp';

% cost_rel,of_p1,of_p2,pref,diam[i]
d1 = importdata(network_name+"_of_p1_log.txt");
d2 = importdata(network_name+"_of_p2_log.txt");

cost_1 = d1(:,1);
of_p1_1 = d1(:,2);
of_p2_1 = d1(:,3);

cost_2 = d2(:,1);
of_p1_2 = d2(:,2);
of_p2_2 = d2(:,3);

% figure;
% plot(cost_1,of_p1_1,'x');
% hold on; grid on;
% plot(cost_2,of_p1_2,'o');
% legend(["Classic","Backup"]);
% xlabel("Relative Cost [-]");
% ylabel("Objective Function 1 [-]");
% 
% figure;
% plot(cost_1,of_p2_1,'x');
% hold on; grid on;
% plot(cost_2,of_p2_2,'o');
% legend(["Classic","Backup"]);
% xlabel("Relative Cost [-]");
% ylabel("Objective Function 1 [-]");

% loading pareto fronts
pf1 = importdata(network_name+"_of_p1_pareto.txt");
pf2 = importdata(network_name+"_of_p2_pareto.txt");


