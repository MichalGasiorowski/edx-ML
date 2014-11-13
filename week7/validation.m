
N = 35;

function [in out] = load_data()
	in = load("in.dta");
	out = load("out.dta");
end

% p in [1 x1 x2] form
function non_lin = non_linear_transform(p)
	non_lin = [p p(:,2).^2 p(:,3).^2 p(:,2).*p(:,3) abs(p(:,2) - p(:,3)) abs(p(:,2) + p(:,3))];
end

function labels = getLabels(features, fun)
	y = fun' * features';
	positive_side = (y >= 0)';
	negative_side = ((y < 0) .* -1)';
	labels = positive_side .+ negative_side;
end

function i_pos = computeMissLabel(points, labels, weight)
	wLabels = getLabels(points, weight);
	misLabel = (labels != wLabels);
	[i_pos j] = find(misLabel != 0);
end

function [train validation] = splitdata(IN, k)
	train = IN(1:k, :);
	validation = IN(k+1:end, :);
end

function [hypo_f, IN_VALID_ERROR, OUT_ERROR] = linearRegression(train_range, valid_range, non_lin_k)
	[IN OUT] = load_data();
	
	IN_train = IN(train_range, :);
	IN_validation = IN(valid_range, :);
	%[IN_train IN_validation] = splitdata(IN, train_size);
	
	N_TRAIN = size(IN_train, 1);
	N_VALIDATION = size(IN_validation, 1);
	N_OUT = size(OUT, 1);
	
	
	IN_X_TRAIN = [ones(N_TRAIN, 1) IN_train(:, 1:2)];
	IN_Y_TRAIN = IN_train(:, 3);
	
	IN_X_VALIDATION = [ones(N_VALIDATION, 1) IN_validation(:, 1:2)];
	IN_Y_VALIDATION = IN_validation(:, 3);
	
	 
	
	IN_train_x_trans = non_linear_transform(IN_X_TRAIN)(:, 1:non_lin_k+1);
	IN_valid_x_trans = non_linear_transform(IN_X_VALIDATION)(:, 1:non_lin_k+1);
	
	OUT_X = [ones(N_OUT, 1) OUT(:, 1:2)];
	OUT_Y = OUT(:, 3);
	
	OUT_TRANS = non_linear_transform(OUT_X)(:, 1:non_lin_k+1);
	
	%size(IN_TRANS)
	%lambda = 10^k;
	%IN_TRANS'*IN_TRANS
	
	hypo_f = pinv(IN_train_x_trans'*IN_train_x_trans )*IN_train_x_trans'*IN_Y_TRAIN;
	
	in_valid_pos = computeMissLabel(IN_valid_x_trans, IN_Y_VALIDATION, hypo_f);
	IN_VALID_ERROR = length(in_valid_pos) / length(IN_Y_VALIDATION);
	
	out_pos = computeMissLabel(OUT_TRANS, OUT_Y, hypo_f);
	OUT_ERROR = length(out_pos) / length(OUT_Y);
	
end

IN_VALID_ERROR_TAB = [];
for k=3:7
	%[hypo_f, IN_VALID_ERROR, OUT_ERROR] = linearRegression(26:35, 1:25, k);
	[hypo_f, IN_VALID_ERROR, OUT_ERROR] = linearRegression(1:25, 26:35, k);
	IN_VALID_ERROR_TAB = [IN_VALID_ERROR_TAB; [k IN_VALID_ERROR OUT_ERROR]];
endfor
IN_VALID_ERROR_TAB

%[hypo_f, IN_VALID_ERROR] = linearRegression(25, 5)








