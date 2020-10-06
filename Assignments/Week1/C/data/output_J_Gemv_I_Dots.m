% number of repeats:% 3
% enter first, last, inc:% 48 480 48 
data = [
%  n          reference      |         current implementation 
%        time       GFLOPS   |    time       GFLOPS     diff 
   480 9.9663e-03 2.2193e+01    2.5982e-01 8.5129e-01 3.5527e-13
   432 6.6198e-03 2.4358e+01    1.4702e-01 1.0967e+00 3.1264e-13
   384 4.3535e-03 2.6013e+01    1.1630e-01 9.7376e-01 1.9895e-13
   336 3.0367e-03 2.4983e+01    6.1924e-02 1.2252e+00 1.7053e-13
   288 1.7365e-03 2.7513e+01    3.8518e-02 1.2403e+00 1.1369e-13
   240 8.5813e-04 3.2219e+01    2.1922e-02 1.2612e+00 4.2633e-14
   192 4.4689e-04 3.1676e+01    1.1463e-02 1.2350e+00 2.8422e-14
   144 1.9759e-04 3.0224e+01    4.4574e-03 1.3398e+00 2.8422e-14
    96 6.9479e-05 2.5468e+01    1.2979e-03 1.3633e+00 1.4211e-14
    48 1.5081e-05 1.4666e+01    1.2717e-04 1.7393e+00 7.1054e-15
];

% Maximum difference between reference and your implementation: 3.552714e-13.