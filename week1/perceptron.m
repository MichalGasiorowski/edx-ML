
N = 100;

function A = generateRandomLine()
% generate random line x2 = a*x1 + b
	points = (rand(2, 2) .- 0.5) .* 2;
	p1 = points(1, :);
	p2 = points(2, :);
	a = ( p2(2) - p1(2) ) / ( p2(1) - p1(1) );
	b = p1(2) - a * p1(1);
	A = [a; b];
end

function P = generatePoints(N)
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

function [it, L, weight, out_error] = runExperiment(N, OUT_ERROR_N)
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
end

RUNS = 1000;
iterations = zeros(1, RUNS);
errors = zeros(1, RUNS);

for i = 1:RUNS
	[it, L, weight, out_error] = runExperiment(N, 10000);
	iterations(i) = it;
	errors(i) = out_error;
end

mean(iterations)
mean(errors)

#{
L = generateRandomLine();
points = generatePoints(N);
pointsEx = [ones(N, 1) points];
labels = getLabels(points, L);

y = L(1) .* points(:, 1) .+ L(2);

[i_pos j] = find(labels == 1);
[i_neg j] = find(labels == -1);

clf;

plot( points(i_pos, 1) , points(i_pos, 2) , 'x' )
axis([-1 1 -1 1]);
hold on;

%plot( x_axis_1, y_axis_1 , "@12"  ) 
plot( points(:, 1) , y , 'g' ) ;

plot( points(i_neg, 1) , points(i_neg, 2) , '@11' )


ITER = 100;
weight = [0 0 0];
misses_count = zeros(1, ITER);
for i = 1:ITER
	missL = computeMissLabel(pointsEx, labels, weight);
	misses_count(i) = length(missL);
	if(misses_count(i) == 0)
		break;
	endif
	ind = ceil(rand(1) * length(missL) ) ;
	p = pointsEx(missL(ind), :);
	
	[p, labels(missL(ind)), weight * p'];
	weight = computeNewWeight(weight, p, labels(missL(ind)) );
	
endfor

weight;
L;
misses_count;
y2 = (-weight(2)/weight(3)) .* points(:, 1) .+ (-weight(1) / weight(3) );
plot( points(:, 1) , y2 , 'k' ) ;
hold off;
#}

