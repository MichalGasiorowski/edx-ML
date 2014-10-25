
%calculate min N such that probability bound is satisfied
% P[ |E_in(g) - E_out(g)| > eps ] <= 2Me^(-2eps^2N)
EMS = [1 10 100];
eps = 0.05;
p_bound = 0.03;
function N = hoeff_calc(M, eps, p_bound)
	N = log(p_bound/(2*M))/(-2*eps*eps);
end



hoeff_calc(EMS(1), eps, p_bound)
hoeff_calc(EMS(2), eps, p_bound)
hoeff_calc(EMS(3), eps, p_bound)
