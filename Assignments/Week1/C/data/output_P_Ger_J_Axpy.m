% number of repeats:% 3
% enter first, last, inc:% 48 480 48 
data = [
%  n          reference      |         current implementation 
%        time       GFLOPS   |    time       GFLOPS     diff 
   480 1.0281e-02 2.1513e+01    6.7356e-02 3.2838e+00 3.5527e-13
   432 7.4396e-03 2.1674e+01    4.6055e-02 3.5011e+00 3.1264e-13
   384 4.7460e-03 2.3862e+01    3.2073e-02 3.5309e+00 1.9895e-13
   336 3.3586e-03 2.2588e+01    1.9357e-02 3.9194e+00 1.7053e-13
   288 2.0316e-03 2.3517e+01    1.1465e-02 4.1670e+00 1.1369e-13
   240 1.1237e-03 2.4604e+01    6.6143e-03 4.1801e+00 4.2633e-14
   192 5.7124e-04 2.4781e+01    3.1729e-03 4.4615e+00 2.8422e-14
   144 2.4792e-04 2.4088e+01    1.1870e-03 5.0309e+00 2.8422e-14
    96 8.8587e-05 1.9974e+01    3.3743e-04 5.2440e+00 1.4211e-14
    48 1.8643e-05 1.1864e+01    4.2752e-05 5.1737e+00 7.1054e-15
];

% Maximum difference between reference and your implementation: 3.552714e-13.
