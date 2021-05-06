clear;

cases = importdata('gamma_cases.txt');
gamma_orig = importdata('gamma_orig.txt');
gamma_orig_std = importdata('gamma_orig_std.txt');
gamma_opt = importdata('gamma_opt.txt');
gamma_opt_std = importdata('gamma_opt_std.txt');

A = [cases,gamma_orig];

