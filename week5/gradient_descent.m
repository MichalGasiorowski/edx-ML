

eta = 0.1;
start_point = [1 1];

function [J, grad] = computeCost(p)
	u = p(1);
	v = p(2);
	grad = zeros(1, 2);
	J = (u*e^v - 2*v*e^(-u))^2;
	grad(1) = 2*e^(-2*u) * (u*e^(u+v) - 2*v) * (e^(u+v) + 2*v);
	grad(2) = 2*e^(-2*u) * (u*e^(u+v) - 2) * (u*e^(u+v) - 2*v);	
end

function new_p = gradient_step(p, eta, grad)
	%[J, grad] = computeCost(p);
	new_p = p .+ (-eta .* grad);
end



E_UV_MIN = 10^(-14);
new_p = start_point;
for iter = 0:100
	[J, grad] = computeCost(new_p);
	if (J < E_UV_MIN)
		break
	endif
	new_p = gradient_step(new_p, eta, grad);	
endfor

format 'long'
disp('iter='), disp(iter)
disp('J='), disp(J)
format

E_UV_MIN = 10^(-14);
new_p = start_point;
Js = zeros(1, 30);
M = [0 1; 1 0];
for iter = 1:length(Js);
	[J, grad] = computeCost(new_p);
	Js(iter) = J;
	t = M(mod(iter, 2) + 1, :);
	new_p = gradient_step(new_p, eta, t .* grad);	
endfor

format 'long'
disp('iter='), disp(iter);
disp('Js='), disp(Js);
format