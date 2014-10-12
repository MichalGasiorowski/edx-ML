clear all;
N = 10;
OUT_ERROR_N = 10000;
EXP_NUMBER = 50;

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
	A = [a; -1; b];
end

function P = generatePoints(N) % generate N random points uniformly from [-1 1] X [1- 1] 
	P = [(rand(N, 2) .- 0.5) .* 2 ones(N, 1) ];
end


function labels = getLabelsW(points, weight)
	%weight
	%points
	y = weight' * points';
	positive_side = (y >= 0)';
	negative_side = ((y < 0) .* -1)';
	labels = positive_side .+ negative_side;
end

function i_pos = computeMissLabel(points, labels, weight)
	labels;
	wLabels = getLabelsW(points, weight);
	misLabel = (labels != wLabels);
	[i_pos j] = find(misLabel != 0);
end

function plotData(points, labels, L, weight)
	y = L(1) .* points(:, 1) .+ L(3);
	[i_pos j] = find(labels == 1);
	[i_neg j] = find(labels == -1);
	
	plot( points(i_pos, 1) , points(i_pos, 2) , 'bx' );
	axis([-1 1 -1 1]);
	xlabel ('x1');
	ylabel ('x2');
	title ('linreg plot');
	hold on;
	plot( points(i_neg, 1) , points(i_neg, 2) , '@11' );
	plot( points(:, 1) , y , 'g' ) ;
	
	y2 = (-weight(2)/weight(3)) .* points(:, 1) .+ (-weight(1) / weight(3) );
	plot( points(:, 1) , y2 , 'k' ) ;
	legend ({'+1', '-1', 'target function', 'linreg-learned function'}, 'location', 'northeastoutside');
	hold off;
end

function data = linearRegression(N, OUT_ERROR_N)
	L = generateRandomLine();
	X = generatePoints(N);
	labels = getLabelsW(X, L);
	w = pinv(X'*X)*X'*labels;
	ERRS = [0 0];
	if (OUT_ERROR_N > 1)
		ERRS = calculateErrors(L, w, X, OUT_ERROR_N);
	else
		ERRS = [0 0];
	endif
	
	data = {L, X, w, ERRS};
	%data = [L ;X ;w ;calculateErrors(L, w, X, OUT_ERROR_N)];
end

function ERRS = calculateErrors(L, w, points, OUT_ERROR_N)
	
	in_hypo_labels = getLabelsW(points, w);
	correct_in_labels = getLabelsW(points, L);
	miss_in_label = (correct_in_labels != in_hypo_labels);
	E_IN = sum(miss_in_label) / size(points, 1);

	out_points = generatePoints(OUT_ERROR_N);
	
	correct_out_labels = getLabelsW(out_points, L);
	hipo_out_labels = getLabelsW(out_points, w);
	
	mis_out_label = (correct_out_labels != hipo_out_labels);
	
	E_OUT = sum(mis_out_label) / OUT_ERROR_N;
	ERRS = [E_IN E_OUT];
end

function data = runExperiment(EXP_NUMBER, N, OUT_ERROR_N)
	data = zeros(EXP_NUMBER, 2);
	for i=1:EXP_NUMBER
		LR = linearRegression(N, OUT_ERROR_N);
		data(i, :) = LR{4};
	endfor
end

function w = computeNewWeight(oldW, point, label)
	w = oldW + point' .* label;
end

function data = runPerceptron(points, L, startW)
	disp("runPerceptron startW=: " ), disp(startW)
	
	weight = startW;
	labels = getLabelsW(points, L)
	MAX_ITER = 1000;
	
	misses_count = zeros(MAX_ITER, 1);
	for it=1:MAX_ITER
		
		missL = computeMissLabel(points, labels, weight);
		
		misses_count(it) = length(missL);
		if(misses_count(it) == 0)
			disp("miss is ZERO")
			break;
		endif
		disp(length(missL));
		ind = ceil(rand(1) * length(missL) ) ;
		p = points(missL(ind), :);
		weight = computeNewWeight(weight, p, labels(missL(ind)) )
	endfor
	
	it = it - 1;
	data = it;
end

function data = runLinPerceptron(EXP_NUMBER, N)
	data = zeros(EXP_NUMBER, 1);
	for i=1:EXP_NUMBER
		LR = linearRegression(N, 0);
		points = LR{2};
		L = LR{1};
		startW = LR{3};
		data(i) = runPerceptron(points, L, startW);
	endfor
end

%data = runExperiment(EXP_NUMBER, N, OUT_ERROR_N);
data = runLinPerceptron(EXP_NUMBER, N);
%L = generateRandomLine();
%X = generatePoints(N);
%labels = getLabelsW(X, L);
%w = pinv(X'*X)*X'*labels;

%plotData(X, labels, L, w);



