clear all;
N = 10;


function A = generateRandomLine()
% generate random line x2 = a*x1 + b
	points = (rand(2, 2) .- 0.5) .* 2;
	p1 = points(1, :);
	p2 = points(2, :);
	if(abs(p1(1) - p2(1)) < 0.001)
		p1(1) = p1(1) /2;
	endif
	a = ( p2(2) - p1(2) ) / ( p2(1) - p1(1) );
	b = p1(2) - a * p1(1);
	A = [b  a -1];
end

function labels = getLabels(points, fun)
	y = fun * points';
	positive_side = (y >= 0)';
	negative_side = ((y < 0) .* -1)';
	labels = positive_side .+ negative_side;
end

function i_pos = computeMissLabel(points, labels, weight)
	wLabels = target_fun(points, weight);
	misLabel = (labels != wLabels);
	[i_pos j] = find(misLabel != 0);
end

function P = generatePoints(N) % generate N random points uniformly from [-1 1] X [-1 1] 
	P = [ones(N, 1) (rand(N, 2) .- 0.5) .* 2 ];
end