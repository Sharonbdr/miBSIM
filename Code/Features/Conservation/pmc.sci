m=1;
L=10;
k=4;
starts = [ 2 9 10 1 ];
tmp = [
1,2,1.554300e-01;
1,9,2.878700e-01;
1,10,3.299000e-01;
1,1,2.268000e-01;
2,2,2.466100e-01;
2,3,2.088400e-01;
2,10,3.485400e-01;
2,1,1.960100e-01;
3,2,3.190200e-01;
3,4,3.272000e-01;
3,10,5.317000e-02;
3,1,3.006100e-01;
4,5,3.190200e-01;
4,9,3.272000e-01;
4,10,5.317000e-02;
4,1,3.006100e-01;
5,2,2.466100e-01;
5,3,2.088400e-01;
5,6,3.485400e-01;
5,1,1.960100e-01;
6,2,2.887600e-01;
6,7,2.439400e-01;
6,10,2.770000e-01;
6,1,1.903000e-01;
7,9,3.272000e-01;
7,10,5.317000e-02;
7,1,3.006100e-01;
8,2,2.466100e-01;
8,3,2.088400e-01;
8,10,3.485400e-01;
8,1,1.960100e-01;
9,2,3.190200e-01;
9,9,3.272000e-01;
9,10,5.317000e-02;
9,1,3.006100e-01;
10,2,2.887600e-01;
10,9,2.439400e-01;
10,10,2.770000e-01;
10,1,1.903000e-01;
];
P=sparse(tmp(:,1:2),tmp(:,3),[10 10]);
tmp = [
7,8,3.190200e-01;
];
Q=sparse(tmp(:,1:2),tmp(:,3),[10 10]);
index=1:L; final=index(full(sum(Q,1))>0);
