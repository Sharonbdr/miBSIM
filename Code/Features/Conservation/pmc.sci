m=1;
L=9;
k=4;
starts = [ 8 2 9 1 ];
tmp = [
1,8,1.760000e-01;
1,2,1.968000e-01;
1,9,3.104000e-01;
1,1,3.168000e-01;
2,8,3.973092e-01;
2,3,2.076296e-01;
2,9,8.754182e-02;
2,1,3.075194e-01;
3,4,3.973092e-01;
3,3,2.076296e-01;
3,9,8.754182e-02;
3,1,3.075194e-01;
4,5,3.231400e-01;
4,2,1.535700e-01;
4,9,2.743800e-01;
4,1,2.489100e-01;
5,6,3.231400e-01;
5,2,1.535700e-01;
5,9,2.743800e-01;
5,1,2.489100e-01;
6,8,3.231400e-01;
6,2,1.535700e-01;
6,1,2.489100e-01;
7,8,3.153200e-01;
7,2,2.205500e-01;
7,9,2.542100e-01;
7,1,2.099200e-01;
8,8,3.231400e-01;
8,2,1.535700e-01;
8,9,2.743800e-01;
8,1,2.489100e-01;
9,8,3.153200e-01;
9,2,2.205500e-01;
9,9,2.542100e-01;
9,1,2.099200e-01;
];
P=sparse(tmp(:,1:2),tmp(:,3),[9 9]);
tmp = [
6,7,2.743800e-01;
];
Q=sparse(tmp(:,1:2),tmp(:,3),[9 9]);
index=1:L; final=index(full(sum(Q,1))>0);
