% number of repeats:% 3
% enter first, last, inc:% 48 480 48 
data = [
%  n          reference      |         current implementation 
%        time       GFLOPS   |    time       GFLOPS     diff 
   480 1.1098e-02 1.9931e+01    6.7673e-02 3.2684e+00 3.5527e-13
   432 6.9743e-03 2.3120e+01    4.1025e-02 3.9304e+00 3.1264e-13
   384 4.5344e-03 2.4975e+01    3.0139e-02 3.7575e+00 1.9895e-13
   336 4.5503e-03 1.6673e+01    1.8560e-02 4.0877e+00 1.7053e-13
   288 1.9124e-03 2.4983e+01    1.1302e-02 4.2273e+00 1.1369e-13
   240 1.0420e-03 2.6532e+01    6.2401e-03 4.4307e+00 4.2633e-14
   192 5.2628e-04 2.6898e+01    3.0328e-03 4.6676e+00 2.8422e-14
   144 2.3516e-04 2.5396e+01    1.1244e-03 5.3112e+00 2.8422e-14
    96 8.1579e-05 2.1690e+01    3.0479e-04 5.8056e+00 1.4211e-14
    48 1.6911e-05 1.3079e+01    3.8878e-05 5.6892e+00 7.1054e-15
];

% Maximum difference between reference and your implementation: 3.552714e-13.
