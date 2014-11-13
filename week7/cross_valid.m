
A = [-1 0];
%X = [p 1];
B = [1 0];


function [a, b] = fit_constant_model(p1, p2)
	a = 0;
	b = (p1(2) + p2(2))/2;
end

function [a, b] = fit_linear_model(p1, p2)
	a = (p1(2) - p2(2)) / (p1(1) - p2(1));
	b = p1(2) - a * p1(1);
end

function err = square_error(p, a, b)
	err = (p(1)*a + b - p(2) ) .^2;
end

function [cross_error_const, cross_error_linear] = cross_valid_error(points)
	p1 = points(1, :);
	p2 = points(2, :);
	p3 = points(3, :);
	
	cross_error_const_tab = zeros(3, 1);
	cross_error_linear_tab = zeros(3, 1);
	
	[a, b] = fit_constant_model(p1, p2);
	cross_error_const_tab(1) = square_error(p3, a, b);
	[a, b] = fit_constant_model(p1, p3);
	cross_error_const_tab(2) = square_error(p2, a, b);
	[a, b] = fit_constant_model(p2, p3);
	cross_error_const_tab(3) = square_error(p1, a, b);
	
	cross_error_const = mean(cross_error_const_tab);
	
	[a, b] = fit_linear_model(p1, p2);
	cross_error_linear_tab(1) = square_error(p3, a, b);
	[a, b] = fit_linear_model(p1, p3);
	cross_error_linear_tab(2) = square_error(p2, a, b);
	[a, b] = fit_linear_model(p2, p3);
	cross_error_linear_tab(3) = square_error(p1, a, b);
	
	cross_error_linear = mean(cross_error_linear_tab);
	
end

ps = [sqrt(sqrt(3) + 4) sqrt(sqrt(3) - 1) sqrt(9 +  4*sqrt(6)) sqrt(9 - sqrt(6))];

cross_valid_error_tab = [];

for i=1:length(ps)
	X = [ps(i) 1];
	[cross_error_const, cross_error_linear] = cross_valid_error([A; B; X]);
	cross_valid_error_tab = [cross_valid_error_tab; [X cross_error_const cross_error_linear] ];
endfor
cross_valid_error_tab








