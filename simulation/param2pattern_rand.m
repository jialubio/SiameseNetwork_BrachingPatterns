function [W, D] = param2pattern_rand(N)

% linear
N1 = 2;   W1 = 1;   D1 = 0.2;
N2 = 8;   W2 = 5;   D2 = 0.2;

W = ((W2 - W1) * N + W1 * N2 - W2 * N1) / (N2 - N1);
D = ((D2 - D1) * N + D1 * N2 - D2 * N1) / (N2 - N1);

