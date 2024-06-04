function [shuffled] = randomize_mat(vec);

shuffled = vec(randperm(length(vec)),:);