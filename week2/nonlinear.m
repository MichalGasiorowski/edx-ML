clear all;
N = 1000;
global noise_factor = 0.1;
global f_weights = [-0.6; 0; 0; 0; 1; 1];
%	value = sign(sum(x(2:3) .* x(2:3)) - 0.6);
global ident_transform = [1 1 1 1 1 1];

more_transform = [1 1 1 1 1 1];
no_transform = [1 1 1 0 0 0];
EXP_NUMBER = 1000;
OUT_ERROR_N =10000;

% transform (1, x1, x2) -> (w1*1, w2*x1, w3*x2, w4*x1x2, w5*x1^2, w6*x2^2)
function non_lin = non_linear_transform(p)
	non_lin = [p p(:,2) .* p(:,3) p(:,2) .* p(:,2) p(:,3) .* p(:, 3)];
end

% features (1, x1, x2, x1x2, x1^2, x2^2)
function labels = getLabels(features, fun)
	y = fun' * features';
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

function output_with_noise = introduce_noise(output)
	global noise_factor;
	L = length(output);
	output_with_noise = output;
	perm = randperm(L);
	indices = ceil(noise_factor * L);
	output_with_noise(perm(1:indices)) *= -1;
end

function data = linearRegression(N, OUT_ERROR_N, transform_w)
	global f_weights;
	global ident_transform;
	raw_X = non_linear_transform(generatePoints(N));
	
	X = raw_X .* ident_transform;
	non_lin_X = raw_X .* transform_w;
	

	correct_labels = getLabels(X, f_weights);
	noise_labels = introduce_noise(correct_labels);
	
	hypo_f = pinv(non_lin_X'*non_lin_X)*non_lin_X'*noise_labels;
	
	ERRS = [0 0];
	if (OUT_ERROR_N > 1)
		ERRS = calculateErrors(f_weights, hypo_f, X, non_lin_X, OUT_ERROR_N, transform_w);
	else
		ERRS = [0 0];
	endif
	
	data = {X, non_lin_X, hypo_f, ERRS};
	%data = [L ;X ;w ;calculateErrors(L, w, X, OUT_ERROR_N)];
end

function ERRS = calculateErrors(correct_f, hypo_f, X, non_lin_X, OUT_ERROR_N, transform_w)
	global f_weights;
	global ident_transform;
	E_IN = 0;
	E_OUT = 0;
	
	%correct_f
	%hypo_f
	%X
	%non_lin_X
	correct_in_labels = getLabels(X, correct_f);
	in_hypo_labels = getLabels(non_lin_X, hypo_f);
	
	
	miss_in_label = (correct_in_labels != in_hypo_labels);
	E_IN = sum(miss_in_label) / size(X, 1);
	
	
	raw_out_X = non_linear_transform(generatePoints(OUT_ERROR_N));
	out_X = raw_out_X .* ident_transform;
	out_non_lin_X = raw_out_X .* transform_w;
	
	correct_out_labels = getLabels(out_X, f_weights);
	
	%noise_labels = introduce_noise(correct_labels);

	hipo_out_labels = getLabels(out_non_lin_X, hypo_f);
	
	mis_out_label = (correct_out_labels != hipo_out_labels);
	
	E_OUT = sum(mis_out_label) / OUT_ERROR_N;
	ERRS = [E_IN E_OUT];
end


function data = runNonLinearExperiment(N, transform_w, EXP_NUMBER, OUT_ERROR_N)
	hypo_f = zeros(EXP_NUMBER, 6);
	ERRS = zeros(EXP_NUMBER, 2);
	for i=1:EXP_NUMBER
		LR = linearRegression(N, OUT_ERROR_N, transform_w);
		hypo_f(i, :) = [LR{3}]';
		ERRS(i, :) = [LR{4}];
		
	endfor
	data = {hypo_f, ERRS};
end

%data = linearRegression(N, 10000, no_transform);

data = runNonLinearExperiment(N, more_transform, EXP_NUMBER, OUT_ERROR_N);

mean(data{1})
mean(data{2})
%data = linearRegression(N, 10000, ident_transform);
%y = rand(10, 1);
%introduce_noise(y)