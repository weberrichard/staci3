clear;

A1 = [4,1,1;2,3,1;3,3,1];
[l1,u1,p1] = lu(A1);

A2 = [21,3,1,3;4,18,3,1;2,4,9,5;5,2,7,11];
[l2,u2,p2] = lu(A2);

A3 = [5 0 0 0 0 0 0 0 -1
0 6 0 0 0 1 0 0 -1
0 0 7 0 0 0 -1 1 0
0 0 0 8 0 0 1 0 -1
0 0 0 0 9 -1 0 1 0
0 1 0 0 -1 14 0 0 0
0 0 -1 1 0 0 13 0 0
0 0 1 0 1 0 0 12 0
-1 -1 0 -1 0 0 0 0 11];

lu_matrix(A3);
[l3,u3,p3] = lu(A3);
[ll3,uu3,pp3] = ldl(A3);

A3 = sparse(A3);
[row, col, v] = find(A3);
dlmwrite('jacobi2_grid.txt',[col row v], 'delimiter', ',');

J11 = A3(1:5,1:5);
J12 = A3(1:5,6:9);
J21 = A3(6:9,1:5);
J22 = A3(6:9,6:9);

lb3 = [eye(5),zeros(5,4);J21*inv(J11),eye(4)];
db3 = [J11,zeros(5,4);zeros(4,5),J22-J21*inv(J11)*J12];
ub3 = [eye(5),inv(J11)*J12;zeros(4,5),eye(4)];

A4 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1
0 0.72 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 -1
0 0 0.72 0 0 0 0 0 0 0 0 0 0 -1 1 0 0 0 0 0 0 0
0 0 0 0.72 0 0 0 0 0 0 0 0 0 0 0 -1 1 0 0 0 0 0
0 0 0 0 0.72 0 0 0 0 0 0 0 0 0 0 0 -1 1 0 0 0 0
0 0 0 0 0 0.72 0 0 0 0 0 0 0 0 0 0 0 0 -1 1 0 0
0 0 0 0 0 0 0.72 0 0 0 0 0 0 0 0 0 0 0 0 -1 1 0
0 0 0 0 0 0 0 0.72 0 0 0 0 0 0 0 1 0 0 0 0 0 -1
0 0 0 0 0 0 0 0 0.72 0 0 0 0 -1 0 0 1 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0.72 0 0 0 0 -1 0 0 1 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0.72 0 0 0 0 -1 0 0 1 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0.72 0 0 0 0 -1 0 0 1 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0.72 0 0 0 0 -1 0 0 1 0
0 1 -1 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 1 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 -1 0 0 0 1 0 0 -1 0 0 0 0 0 0 0 0 0 0 0
0 0 0 1 -1 0 0 0 1 0 0 -1 0 0 0 0 0 0 0 0 0 0
0 0 0 0 1 0 0 0 0 1 0 0 -1 0 0 0 0 0 0 0 0 0
0 0 0 0 0 -1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 1 -1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0
-1 -1 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
d1 = zeros(size(A4));
d2 = randi(32,1,length(A4));
for i=1:length(A4)
   d1(i,i) = d2(i); 
end
A4 = A4 + d1;
[l4,u4,p4] = lu(A4);

A4 = sparse(A4);
[row, col, v] = find(A4);
dlmwrite('jacobi2_ohermes.txt',[col row v], 'delimiter', ',');

[L1,U1] = lu_matrix(A1);
[L2,U2] = lu_matrix(A2);

Ai = [3,3,3]; % non-zero elements in each row
Aj = [1,2,3,1,2,3,1,2,3]; % column index of each elements
Ax = [4,1,1,2,3,1,3,3,1]; % element values
n = 3;

A5 = readmatrix('jacobi_ferto.txt');
A5 = A5(:,1:end-1);
A5 = sparse(A5);
[row, col, v] = find(A5);
dlmwrite('jacobi2_ferto.txt',[col row v], 'delimiter', ',');

tic;
[l5,u5,p5] = lu(A5);
toc;
tic;
[l6,u6,p6] = ldl(A5);
toc;
figure;
spy(A5,'x');
hold on
spy(l5,'o');
spy(u5,'+');

function [L,U] = lu_matrix(A)
    n = length(A);

    U = zeros(n);
    L = eye(n);
    for i=1:n
        %f3m
        for j=i:n
            U(i,j) = A(i,j) - sum(L(i,1:i-1)'.*U(1:i-1,j));
        end

        %a3m
        for j=i+1:n
            L(j,i) = (A(j,i) - sum(L(j,1:i-1)'.*U(1:i-1,i)))/U(i,i);
        end
    end
end