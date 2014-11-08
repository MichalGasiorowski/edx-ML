
N_in = 35;


function [in out] = load_data()
	in = load("in.dta");
	out = load("out.dta");
end

% p in [1 x1 x2] form
function non_lin = non_linear_transform(p)
	non_lin = [p p(:,2).*p(:,3) p(:,2).^2 p(:,3).^2 abs(p(:,2) - p(:,3)) abs(p(:,2) + p(:,3))];
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


function [hypo_f, IN_ERROR, OUT_ERROR] = linearRegression(k)
	[IN OUT] = load_data();
	IN_N = size(IN, 1);
	IN_X = [ones(IN_N, 1) IN(:, 1:2)];
	IN_Y = IN(:, 3);
	IN_TRANS = non_linear_transform(IN_X);
	OUT_X = [ones(size(OUT, 1), 1) OUT(:, 1:2)];
	OUT_Y = OUT(:, 3);
	OUT_TRANS = non_linear_transform(OUT_X);
	
	%size(IN_TRANS)
	lambda = 10^k;
	%IN_TRANS'*IN_TRANS
	
	hypo_f = pinv(IN_TRANS'*IN_TRANS .+ lambda*eye(size(IN_TRANS, 2)) )*IN_TRANS'*IN_Y;
	in_pos = computeMissLabel(IN_TRANS, IN_Y, hypo_f);
	IN_ERROR = length(in_pos) / length(IN_Y);
	out_pos = computeMissLabel(OUT_TRANS, OUT_Y, hypo_f);
	OUT_ERROR = length(out_pos) / length(OUT_Y);
	
	
end

[hypo_f, IN_ERROR, OUT_ERROR] = linearRegression(-1);

outs = [];
for k = -10:10
	[hypo_f, IN_ERROR, OUT_ERROR] = linearRegression(k);
	outs = [outs; k OUT_ERROR];
endfor

outs







