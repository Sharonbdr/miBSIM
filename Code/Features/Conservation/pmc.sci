m=1;
L=9;
k=4;
starts = [ 2 8 9 1 ];
tmp = [
1,2,7.485000e-02;
1,8,3.742500e-01;
1,9,4.311400e-01;
1,1,1.197600e-01;
2,2,1.613787e-01;
2,3,2.962976e-01;
2,9,4.470864e-01;
2,1,9.523724e-02;
3,2,2.172200e-01;
3,4,3.933200e-01;
3,9,1.812300e-01;
3,1,2.082300e-01;
4,5,2.172200e-01;
4,8,3.933200e-01;
4,9,1.812300e-01;
4,1,2.082300e-01;
5,2,1.613787e-01;
5,3,2.962976e-01;
5,6,4.470864e-01;
5,1,9.523724e-02;
6,2,1.759682e-01;
6,9,3.490665e-01;
6,1,1.373386e-01;
7,2,2.172200e-01;
7,8,3.933200e-01;
7,9,1.812300e-01;
7,1,2.082300e-01;
8,2,2.172200e-01;
8,8,3.933200e-01;
8,9,1.812300e-01;
8,1,2.082300e-01;
9,2,1.759682e-01;
9,8,3.376266e-01;
9,9,3.490665e-01;
9,1,1.373386e-01;
];
P=sparse(tmp(:,1:2),tmp(:,3),[9 9]);
tmp = [
6,7,3.376266e-01;
];
Q=sparse(tmp(:,1:2),tmp(:,3),[9 9]);
index=1:L; final=index(full(sum(Q,1))>0);
