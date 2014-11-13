
K = 10^6;

C = [rand(K, 1) rand(K, 1)];
D = min(C, [], 2);

mean(C(:, 1))
mean(C(:, 2))
mean(D)