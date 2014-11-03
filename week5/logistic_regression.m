clear all;
N = 100;
eta = 0.01;
eps = 0.01;

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
	wLabels = target_fun(points, weight);
	misLabel = (labels != wLabels);
	[i_pos j] = find(misLabel != 0);
end

function P = generatePoints(N) % generate N random points uniformly from [-1 1] X [-1 1] 
	P = [ones(N, 1) (rand(N, 2) .- 0.5) .* 2 ];
end

function plotData(points, labels, L, weights)
	%y = L(1:2) * points(:,1:2)';
	y = (-L(2)/L(3)) .* points(:, 2) .+ (-L(1) / L(3) );
	[i_pos j] = find(labels == 1);
	[i_neg j] = find(labels == -1);
	
	plot( points(i_pos, 2) , points(i_pos, 3) , 'bx' );
	axis([-1 1 -1 1]);
	xlabel ('x1');
	ylabel ('x2');
	title ('Logistic Regression');
	hold on;
	plot( points(i_neg, 2) , points(i_neg, 3) , '@11' );
	plot( points(:, 2) , y , 'g' );
	
	y2 = (-weights(2)/weights(3)) .* points(:, 2) .+ (-weights(1) / weights(3) );
	plot( points(:, 2) , y2 , 'k' ) ;

	legend ({'+1', '-1', 'target function','Logistic Regression function'}, 'location', 'northeastoutside');
	hold off;
end

function partial = partial_Ein(xn, yn, w)
	partial = (-yn*xn)/(1 + e^(yn*xn*w));
end

function cross_E_err = compute_X_Err(x, y, w)
	N = size(x, 1);
	cross_E_err = (1.0/N) .* sum(log(1 .+ e.^(-y .* x*w')));
end

function [weights, cross_E_err, iter] = runLogisticRegression(N, eta, eps,max_epoch)
	
	points = generatePoints(N);
	L = generateRandomLine();
	labels = getLabels(points, L);
	%plotData(points, labels, L);
	Ws = zeros(10, 3);
	w = [0; 0; -1];
	partial = 0;
	eps_err = 1;
	for iter=1:max_epoch
		perm = randperm(N);
		points = points(perm,:);
		labels = labels(perm);
		for i=1:N
			xn = points(i,:);
			yn = labels(i);
			partial = partial_Ein(xn, yn, w);
			w = w .- (eta * partial');
		endfor
		Ws(iter,:) = w;
		
		if(iter > 1)
			%eps_err = sum((Ws(iter-1,:) .- Ws(iter,:)).^2)
			eps_err = norm(Ws(iter-1,:) .- Ws(iter,:));
			if(eps_err < eps)
				break
			endif
		endif
	endfor
	
	%disp('eps_err'), eps_err
	weights = Ws(iter, :);
	plotData(points, labels, L, weights)
	%weights
	%L
	out_points = points = generatePoints(10*N);
	out_labels = getLabels(out_points, L);
	cross_E_err = compute_X_Err(out_points, out_labels, weights);
end

errs_arr = [];
iter_arr = [];
for i =1:10
	[weights, cross_E_err, iter] = runLogisticRegression(N, eta, eps,1000);
	iter_arr = [iter_arr iter];
	errs_arr = [errs_arr cross_E_err];
endfor
disp('mean(iter_arr)'), mean(iter_arr)
disp('mean(errs_arr)'), mean(errs_arr)






