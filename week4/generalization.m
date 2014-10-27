clear all;

% P( |E_in - E_out| > eps ) <= 4 * m_h(2N) * exp(-1/8*eps*eps*N)
% M_H(N) = N**d_vc
% M_H(2N) = (2N)**d_vc = 
d_vc = 10;
delta = 0.95;
eps = 0.05;

Ns = [10**5; 400000; 420000; 440000; 460000; 480000];

%function delta = confidence_bound(N, eps, d_vc)
	%delta = 4 * (realpow((2*N), d_vc)) * exp(-1/8 * eps * eps * N);
%end

function epsilon = VC_bound(N, delta, d_vc)
	%epsilon = sqrt(8/N * log(4 * ((2* N) ** d_vc) / delta ) ); 
	epsilon = sqrt(8/N * ( log(4)  -log(delta) + 2*d_vc*log(N) ) );
end

function epsilon = RademacherPenalty_bound(N, delta, d_vc)
	epsilon = sqrt(2/N * ( log(2) + log(N) + d_vc*log(N) )  ) + sqrt(2/N * log(1/delta)) + 1/N; 
end


% eps <= sqrt ( 1/N * (2*eps + log(6 * m_h(2N) / delta ) )
function epsilon = PVdB_bound(N, delta, d_vc, eps)
	%epsilon = sqrt(1/N * ( 2*eps + log(6 * ( (2*N)**d_vc ) / delta   ) ) ) ;
	epsilon = sqrt(1/N * ( 2*eps + log(6) + 2*d_vc*log(N) - log(delta) ) ) ;
end

function epsilon = devroye_bound(N, delta, d_vc, eps)
	epsilon = sqrt( 1/(2*N) * ( 4*eps*(1+eps) + log(4) - log(delta) + 2*d_vc*log(N) ) ) ;
end

function epsilon = calc_iter(estimate_f, N, delta, d_vc, init_eps, iter)
	old_eps = epsilon = init_eps;
	for i=1:iter
		epsilon = estimate_f(N, delta, d_vc, epsilon);
		if(abs(old_eps - epsilon) < 0.000001)
			break;
		endif
	endfor
end

epsis = arrayfun(@(x) VC_bound(x, delta, d_vc), Ns);

rr = [Ns epsis epsis <= eps];

delta = 0.05;
d_vc= 50;

Ns = 100:10000;
%Ns = 1:10;

%PVdB_bound(10000, delta, d_vc, 1);
%calc_PVdB(10000, delta, d_vc, 1, 20);
%error('I dont want to play anymore')


y_vc = arrayfun(@(x) VC_bound(x, delta, d_vc), Ns);
y_rp_vc = arrayfun(@(x) RademacherPenalty_bound(x, delta, d_vc), Ns);
y_PVdB_vc = arrayfun(@(x) calc_iter(@PVdB_bound, x, delta, d_vc, 1, 20), Ns);
y_Devroye_vc = arrayfun(@(x) calc_iter(@devroye_bound, x, delta, d_vc, 1, 20), Ns);

plot(Ns, y_vc, 'r');
hold on;
plot(Ns, y_rp_vc, 'b');
plot(Ns, y_PVdB_vc, 'g');
plot(Ns, y_Devroye_vc, 'm');



xlabel('N');
ylabel('generalization error');
legend('VC Bound', 'Rademacher Penalty Bound', 'Parrondo and Van den Broek Bound', 'Devroye Bound');

hold off;
