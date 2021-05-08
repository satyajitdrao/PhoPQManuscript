%     Q Qp  P  Pp QpP QPp
A1 =[-1  1  0  0  0  0
      0 -1 -1  0  1  0
      1  0  0  1 -1  0
     -1  0  0 -1  0  1
      1  0  1  0  0 -1]; % without X
A = [A1 zeros(5,5)];
B =[0 0  0  0 0 0 0  1 -1  0  0
    0 0  1  0 0 0 0  0  1 -1  0
    0 0  0 -1 0 0 0 -1  0  1  0
    0 0  0  1 0 0 0  1  0  0 -1
    0 0 -1  0 0 0 0 -1  0  0  1
    1 0  0  0 0 0 1 -1  0  0  0];
B1 = [A; B]; % with X
rank(B1)
rank(A1)
C1 = B1; C1(:,7) = []; % pseudo first order
rank(C1) 