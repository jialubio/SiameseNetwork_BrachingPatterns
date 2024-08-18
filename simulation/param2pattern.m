function [W, D] = param2pattern(N)
% load optimal width and density
load([pwd '/results.mat'])
% linear
N_range = N0V;
W_range = optimWidth;
D_range = optimDensity;

W = zeros(size(N));
D = zeros(size(N));

for i = 1 :size(N,1),
    for j = 1 :size(N,2),
        % find N range
        diff = N_range - N(i,j);
        idx = min(find(diff >= 0));
        
        N1 = N_range(idx-1); W1 = W_range(idx-1); D1 = D_range(idx-1);
        N2 = N_range(idx);   W2 = W_range(idx);   D2 = D_range(idx);
        
        W(i,j) = ((W2 - W1) * N(i,j) + W1 * N2 - W2 * N1) / (N2 - N1);
        D(i,j) = ((D2 - D1) * N(i,j) + D1 * N2 - D2 * N1) / (N2 - N1);
    end
end

