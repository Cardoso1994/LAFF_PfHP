% number of repeats:% 3
% enter first, last, inc:% 48 480 48 
data = [
%  n          reference      |         current implementation 
%        time       GFLOPS   |    time       GFLOPS     diff 
   480 3.0944e-02 7.1478e+00    6.1850e-02 3.5761e+00 3.5527e-13
   432 7.4088e-03 2.1764e+01    4.5445e-02 3.5481e+00 3.1264e-13
   384 5.3038e-03 2.1352e+01    2.8724e-02 3.9426e+00 1.9895e-13
   336 3.3634e-03 2.2557e+01    1.8685e-02 4.0602e+00 1.7053e-13
   288 1.9186e-03 2.4902e+01    1.0910e-02 4.3792e+00 1.1369e-13
   240 1.0832e-03 2.5524e+01    6.0479e-03 4.5715e+00 4.2633e-14
   192 5.9291e-04 2.3875e+01    3.0313e-03 4.6698e+00 2.8422e-14
   144 2.3632e-04 2.5270e+01    1.2278e-03 4.8641e+00 2.8422e-14
    96 8.2048e-05 2.1566e+01    3.6070e-04 4.9057e+00 1.4211e-14
    48 1.7610e-05 1.2560e+01    4.0597e-05 5.4483e+00 7.1054e-15
];

% Maximum difference between reference and your implementation: 3.552714e-13.
