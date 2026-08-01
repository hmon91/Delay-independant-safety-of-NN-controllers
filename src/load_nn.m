function net = load_nn(dataDir)
%LOAD_NN  Load the trained tanh FFNN controller (1-10-10-1) used in the paper.
%   net = LOAD_NN() reads W1.csv, W2.csv, W3.csv from ../data and returns a
%   network struct with zero biases and tanh/tanh/linear activations.
%   net = LOAD_NN(dataDir) reads the weights from a custom directory.
if nargin < 1 || isempty(dataDir)
    dataDir = fullfile(fileparts(mfilename('fullpath')), '..', 'data');
end
W1 = readmatrix(fullfile(dataDir, 'W1.csv'));   % 10 x 1  (input  -> hidden 1)
W2 = readmatrix(fullfile(dataDir, 'W2.csv'));   % 10 x 10 (hidden 1 -> hidden 2)
W3 = readmatrix(fullfile(dataDir, 'W3.csv'));   % 1  x 10 (hidden 2 -> output)
net.W   = {W1, W2, W3};
net.b   = {zeros(size(W1,1),1), zeros(size(W2,1),1), zeros(size(W3,1),1)};
net.act = {'tanh', 'tanh', 'linear'};
end
