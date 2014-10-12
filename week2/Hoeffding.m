
% 1000 virtual coins
% Flip each coin independently 10 times.
% Focus on 3 coins:
% c_1 		-> first coin flipped
% c_rand	-> coin chosen at random from 1000 coins
% c_min 	-> coin which had minimum freq. of head (pick the earlier one)
% v_1, v_rand, v_min -> fraction of heads obtained for these respective coins
% run experiment 10**5 times


N = 1000; %number of coin flips
M = 10; % each coin flipped this many times
NN = 10^5;

% 0 -> tails 1 -> heads
% fraction of heads is mean
function flips = runExperiment(coinNum, repeatFlips)
	flips = zeros(coinNum, repeatFlips);
	
	flips = ceil(rand(coinNum, repeatFlips)*2 -1);

end

% calculate v1, v_rand, v_min for specific experiment
function stats = calculateStats(flips)
	avg = mean(flips, 2);
	v_1 = avg(1);
	rand_idx = ceil(rand(1) * size(flips, 1));
	v_rand = avg(rand_idx);
	v_min = min(avg);
	stats = [v_1, v_rand, v_min];
	
end

function megaStats = runBigExperiment(howManyTimes, coinNum, repeatFlips)
	mm = zeros(howManyTimes, 3);
	for i =1:howManyTimes
		mm(i, :) = calculateStats(runExperiment(coinNum, repeatFlips));	
	endfor
	
	megaStats = mean(mm);

end

RR = runBigExperiment(NN, N, M);

