
sigma = 0.1; %variance
d = 8; % dimensions

%expected in-sample error
function E_IN = calc_e_in(sigma, d, N)
	E_IN = sigma ^2 * (1 - (d + 1)/N);
end

Ens = [10 25 100 500 1000];
out = zeros(length(Ens), 1);
E_IN_BOUNDARY = 0.008;

for i=1:length(Ens)
	out(i) = calc_e_in(sigma, d, Ens(i));
endfor

out > E_IN_BOUNDARY