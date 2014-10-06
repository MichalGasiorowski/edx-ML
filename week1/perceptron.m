clear all;
N = 100;

function A = generateRandomLine()
% generate random line x2 = a*x1 + b
	points = (rand(2, 2) .- 0.5) .* 2;
	p1 = points(1, :);
	p2 = points(2, :);
	if(abs(p1(1) - p2(1)) < 0.0001)
		p1(1) = p1(1) /2;
	endif
	a = ( p2(2) - p1(2) ) / ( p2(1) - p1(1) );
	b = p1(2) - a * p1(1);
	A = [a; b];
end

function P = generatePoints(N) % generate N random points uniformly from [-1 1] X [1- 1] 
	P = (rand(N, 2) .- 0.5) .* 2;
end

function labels = getLabels(points, L)
	y = L(1) .* points(:, 1) .+ L(2);
	positive_side = (points(:, 2) >= y);
	negative_side = (points(:, 2) < y) .* -1;
	labels = positive_side .+ negative_side;
end

function labels = getLabelsW(points, weight)
	%weight
	%points
	y = weight * points';
	positive_side = (y >= 0)';
	negative_side = ((y < 0) .* -1)';
	labels = positive_side .+ negative_side;
end

function w = computeNewWeight(oldW, point, label)
	w = oldW + point .* label;
end

function i_pos = computeMissLabel(points, labels, weight)
	wLabels = getLabelsW(points, weight);
	misLabel = (labels != wLabels);
	[i_pos j] = find(misLabel != 0);
end

function [it, points, labels, L, weight, out_error] = runExperiment(N, OUT_ERROR_N)
	L = generateRandomLine();
	points = generatePoints(N);
	pointsEx = [ones(N, 1) points];
	labels = getLabels(points, L);
	
	y = L(1) .* points(:, 1) .+ L(2);
	
	MAX_ITER = 1000;
	weight = [0 0 0];
	
	for it = 1:MAX_ITER
		missL = computeMissLabel(pointsEx, labels, weight);
		misses_count(it) = length(missL);
		if(misses_count(it) == 0)
			break;
		endif
		ind = ceil(rand(1) * length(missL) ) ;
		p = pointsEx(missL(ind), :);
		
		weight = computeNewWeight(weight, p, labels(missL(ind)) );
	endfor
	
	% calculate error rate
	% P[f(x) != g(x)]
	
	out_points = generatePoints(OUT_ERROR_N);
	out_pointsEx = [ones(OUT_ERROR_N, 1) out_points];
	
	labels2 = getLabels(out_points, L);
	wLabels2 = getLabelsW(out_pointsEx, weight);
	misLabel = (labels2 != wLabels2);
	out_error = sum(misLabel) / OUT_ERROR_N;
	it = it - 1;
end

function plotData(points, labels, L, weight)
	y = L(1) .* points(:, 1) .+ L(2);
	[i_pos j] = find(labels == 1);
	[i_neg j] = find(labels == -1);
	
	plot( points(i_pos, 1) , points(i_pos, 2) , 'bx' );
	axis([-1 1 -1 1]);
	xlabel ('x1');
	ylabel ('x2');
	title ('Perceptron plot');
	hold on;
	plot( points(i_neg, 1) , points(i_neg, 2) , '@11' );
	plot( points(:, 1) , y , 'g' ) ;
	
	y2 = (-weight(2)/weight(3)) .* points(:, 1) .+ (-weight(1) / weight(3) );
	plot( points(:, 1) , y2 , 'k' ) ;
	legend ({'+1', '-1', 'target function', 'perceptron-learned function'}, 'location', 'northeastoutside');
	hold off;
end

function [iterations, errors] = runIterations(N, run_num)
	iterations = zeros(1, run_num);
	errors = zeros(1, run_num);
	for i = 1:run_num
		[it, points, labels, L, weight, out_error] = runExperiment(N, 10000);
		iterations(i) = it;
		errors(i) = out_error;
	end
end

RUNS = 1000;

%[iterations, errors] = runIterations(N, RUNS);
%mean_iter = mean(iterations)
%mean_error = mean(errors)


[it, points, labels, L, weight, out_error] = runExperiment(N, 10000);
plotData(points, labels, L, weight);




