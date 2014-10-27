clear all; close all;
N = 10000;

function points = draw_points(N)
	points = ([rand(N, 1)] .- 0.5) .* 2;
	points = [points sin(pi*points(:, 1) )];
	points = [points ([rand(N, 1)] .- 0.5) .* 2];
	points = [points sin( pi*points(:, 3) ) ];
	
	%near = find(abs(points(:, 3) - points(:,1)) < 0.000001);
	%points(near, :) = 
end

function a = find_best_ax(p1, p2)
	a = ( p1(1) * p1(2) + p2(1) * p2(2) ) / ( p1(1)**2 + p2(1)**2 ) ;	
end

function a = find_best_ax_row(row)
	a = ( row(1) * row(2) + row(3) * row(4) ) / ( row(1)**2 + row(3)**2 ) ;	
end

function b = find_best_b_row(row)
	b = mean(row([2 4]));
end

function out = find_best_ax_b_row(row)
	x1 = row(1); y1 = row(2); x2 = row(3); y2 = row(4);

	a = (y2 - y1) / (x2 - x1);
	b = y1 - a*x1;
	out = [a b];
end

function out = find_best_axx_b_row(row)
	x1 = row(1); y1 = row(2); x2 = row(3); y2 = row(4);
	
	a = (y2 - y1) / (x2**2 - x1**2);
	b = y1 - a*x1**2;
	out = [a b];
end

function a = find_best_axx_row(row)
	x1 = row(1); y1 = row(2); x2 = row(3); y2 = row(4);
	
	a = (2*(y1*x1**2 + y2*x2**2))/(x1**4 + x2**4);
end



function out = calc_out_error(learn_f, N)
	points = draw_points(N);
	
	applyToGivenRow = @(func, matrix) @(row) func(matrix(row, :));
	applyToRows  = @(func, matrix) arrayfun(applyToGivenRow(func, matrix), 1:size(matrix,1))';
	
	gs = applyToRows(learn_f, points);
	g_bar = mean(gs);
	
	points = draw_points(N);
	
	bias = sum((points(:,2) - g_bar*points(:,1)) .**2)/N ; 
	variance = mean((gs .* points(:,1) - g_bar*points(:,1) ) .**2);
	out_error = bias + variance;
	out = [bias, variance, out_error];
end

function out = calc_out_error2(learn_f, N) % 2 parameters
	points = draw_points(N);
	
	gs = zeros(N, 2);
	
	for i=1:N
		gs(i,:) = learn_f(points(i, :));
	endfor
	
	g_bar = mean(gs);
	
	points = draw_points(N);
	
	bias = sum((points(:,2) - g_bar(1)*points(:,1) - g_bar(2)) .**2)/N ;
	variance = mean((gs(:,1) .* points(:,1) + gs(:,2) - g_bar(1)*points(:,1) - g_bar(2)) .**2);
	
	out_error = bias + variance;
	out = [bias, variance, out_error];
end



%points = ([rand(N, 1) rand(N, 1) rand(N, 1) rand(N, 1)] .- 0.5) .* 2;


%applyToGivenRow = @(func, matrix) @(row) func(matrix(row, :));
%applyToRows  = @(func, matrix) arrayfun(applyToGivenRow(func, matrix), 1:size(matrix,1), 'UniformOutput', false)';


f = @(x) (sin(pi*x) - g_bar*x)**2;

%quadv(f, -1, 1)

%gs = arrayfun(@(p1, p2) 

out_b  = calc_out_error(@find_best_b_row, N)
out_a  = calc_out_error(@find_best_ax_row, N)
out_ax_b  = calc_out_error2(@find_best_ax_b_row, N)
out_ax2  = calc_out_error(@find_best_axx_row, N)
out_ax2_b  = calc_out_error2(@find_best_axx_b_row, N)




