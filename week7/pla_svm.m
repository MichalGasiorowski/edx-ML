%clear all;
N = 100;
COMPARE_RUNS = 20;

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

function labels = getLabels(features, fun)
	y = fun * features';
	positive_side = (y >= 0)';
	negative_side = ((y < 0) .* -1)';
	labels = positive_side .+ negative_side;
end

function i_pos = computeMissLabel(points, labels, weight)
	wLabels = getLabels(points, weight);
	misLabel = (labels != wLabels);
	[i_pos j] = find(misLabel != 0);
end

function P = generatePoints(N) % generate N random points uniformly from [-1 1] X [-1 1] 
	P = [ones(N, 1) (rand(N, 2) .- 0.5) .* 2 ];
end

function w = computeNewWeight(oldW, point, label)
	w = oldW .+ point .* label;
end

function [weight, it, in_error] = runPLA(points, labels, MAX_ITER)
	weight = [0 0 0];
	misses_count = 0;
	for it = 1:MAX_ITER
		missL = computeMissLabel(points, labels, weight);
		misses_count = length(missL);
		if (misses_count == 0)
			break;
		endif
		rand_ind = ceil(rand(1) * length(missL) ) ;
		p = points(missL(rand_ind), :);
		weight = computeNewWeight(weight, p, labels(missL(rand_ind)) );
	endfor
	in_error = misses_count / size(points, 1);
	it = it - 1;
end

function [weight, in_error, support_vector_count] = runSVM(points, labels)
	weight = 0;
	in_error = 0;
	support_vector_count = 0;

	X = points(:, 2:end);
	%X = points;
	Y = labels;
	
	N = size(X, 1);
	H = (Y*Y') .* (X*X');
	L = ones(N ,1) .* -1;
	A = Y';
	q = -1 * ones(N, 1);
	b = 0;
	
	lb = zeros(N, 1);
	ub = ones(N, 1) .* 10^22;
	opt = optimset ("MaxIter", 500);
	
	[alpha, obj, info, lambda] = qp([],H,L,A,b,lb,ub,opt);
	sup_vect = (alpha >= 10^-5) .+ (alpha <= (-1)*10^-5);
	support_vector_count = sum( sup_vect );
	sup_vect_ind= find(sup_vect == 1);
	
	%A = [b  a -1];
	wt = (alpha .* Y)' * X;
	X1 = X(sup_vect_ind(1), :);
	bb = 1/Y(sup_vect_ind(1)) - wt * X1';
	%Y(sup_vect_ind(1)) * wt * X1'
	
	weight = [bb (alpha .* Y)' * X ];
	in_error = compute_error_on_points(weight, labels, points);
end

function err = compute_error(hypo_fun, model_fun, N_SIZE)
	out_points = generatePoints(N_SIZE); 
	true_labels = getLabels(out_points, model_fun);
	hypo_labels = getLabels(out_points, hypo_fun);
	
	misLabel = (true_labels != hypo_labels);
	err = sum(misLabel) / N_SIZE;
end

function err = compute_error_on_points(hypo_fun, true_labels, points)
	N_SIZE = size(true_labels, 1);
	hypo_labels = getLabels(points, hypo_fun);
	
	misLabel = (true_labels != hypo_labels);
	err = sum(misLabel) / N_SIZE;
end

function plotPoints(points, labels)
	[i_pos j] = find(labels == 1);
	[i_neg j] = find(labels == -1);
	
	plot( points(i_pos, 2) , points(i_pos, 3) , 'bx' );
	axis([-1 1 -1 1]);
	xlabel ('x1');
	ylabel ('x2');
	title ('PLA plot');
	hold on;
	plot( points(i_neg, 2) , points(i_neg, 3) , '@11' );
	legend ({'+1', '-1'}, 'location', 'northeastoutside');
	
end

function plotY(points, L, color)
	y = (-L(2)/L(3)) .* points(:, 2) .+ (-L(1) / L(3) );
	plot( points(:, 2) , y , color );
end

function [hypo_funs, test_errors, out_errors ] = runPLAExperiment(N, runs)
	test_errors = zeros(runs, 1);
	out_errors = zeros(runs, 1);
	hypo_funs = [];
	points_all_on_same_side = true;
	points = 0;
	labels = 0;
	for i=1:runs
		model_fun = generateRandomLine();
		points_all_on_same_side = true;
		while (points_all_on_same_side)
			points = generatePoints(N);
			labels = getLabels(points, model_fun);
			points_all_on_same_side = (abs(sum(labels)) == N);
		endwhile
		[hypo_fun, it, in_error] = runPLA(points, labels, 1000);
		out_error_it = compute_error(hypo_fun, model_fun, 10000);
		out_errors(i) = out_error_it;
		test_errors(i) = in_error;
		hypo_funs = [hypo_funs; [hypo_fun]];
	endfor
	
end

function [pla_test_errors, svm_test_errors, pla_out_errors, svm_out_errors, support_vector_counts ] = compare_PLA_SVM(N, runs)
	pla_test_errors = zeros(runs, 1);
	pla_out_errors = zeros(runs, 1);
	svm_test_errors = zeros(runs, 1);
	svm_out_errors = zeros(runs, 1);
	support_vector_counts = zeros(runs, 1);
	ten_percent = ceil(runs / 10);
	

	disp('compare_PLA_SVM start')
	for i=1:runs
		if(mod(i, ten_percent) == 0)
			disp(i);
			fflush(stdout);
		endif
		model_fun = generateRandomLine();
		points_all_on_same_side = true;
		while (points_all_on_same_side)
			points = generatePoints(N);
			labels = getLabels(points, model_fun);
			points_all_on_same_side = (abs(sum(labels)) == N);
		endwhile
		
		[pla_hypo_fun, it, pla_in_error] = runPLA(points, labels, 100);
		[svm_weight, svm_in_error, support_vector_count] = runSVM(points, labels);
		
		pla_test_errors(i) = pla_in_error;
		svm_test_errors(i) = svm_in_error;
		
		out_points = generatePoints(1000); 
		true_labels = getLabels(out_points, model_fun);
		
		pla_out_error_it = compute_error_on_points(pla_hypo_fun, true_labels, out_points);
		pla_out_errors(i) = pla_out_error_it;
		
		svm_out_error_it = compute_error_on_points(svm_weight, true_labels, out_points);
		svm_out_errors(i) = svm_out_error_it;
		support_vector_counts(i) = support_vector_count;
		
	endfor
	
end



%[hypo_funs, test_errors, out_errors ] = runPLAExperiment(N, 1000);
%mean(test_errors)
%mean(out_errors)

#{
points_all_on_same_side = true;
model_fun = generateRandomLine();
while (points_all_on_same_side)
	points = generatePoints(N);
	labels = getLabels(points, model_fun);
	points_all_on_same_side = (abs(sum(labels)) == N);
endwhile

[svm_weight, svm_in_error, support_vector_count] = runSVM(points, labels);
[pla_hypo_fun, it, pla_in_error] = runPLA(points, labels, 10000);

model_fun
svm_weight
pla_hypo_fun

plotPoints(points, labels)
plotY(points, pla_hypo_fun, 'g')
plotY(points, svm_weight, 'k')
hold off;
#}

%#{
[pla_test_errors, svm_test_errors, pla_out_errors, svm_out_errors, support_vector_counts ] = compare_PLA_SVM(N, COMPARE_RUNS);

%disp(fprintf('\nMean pla_test_error: %f ', mean(pla_test_errors)));
%disp(fprintf('\nMean svm_test_errors: %f ', mean(svm_test_errors)));
%disp(fprintf('\nMean pla_out_errors: %f ', mean(pla_out_errors)));\
%disp(fprintf('\nMean svm_out_errors: %f ', mean(svm_out_errors)));

fprintf('\n');

disp(mean(pla_test_errors));
disp(mean(svm_test_errors));
disp(mean(pla_out_errors));
disp(mean(svm_out_errors));
disp(mean(support_vector_counts));

fprintf('\n');

fflush(stdout) 

svm_better_pla = sum(svm_out_errors < pla_out_errors) / COMPARE_RUNS;
svm_better_pla
%#}












