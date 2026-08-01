function u = nn_forward(y, net)
%NN_FORWARD  Forward pass of a fully-connected (feedforward) NN controller.
%   u = NN_FORWARD(y, net) evaluates u = Phi(y) for the network described by the
%   struct NET with fields:
%       net.W   - 1xL cell of weight matrices  {W1, W2, ..., WL}
%       net.b   - 1xL cell of bias vectors     {b1, b2, ..., bL}
%       net.act - 1xL cell of activation names {'tanh','relu','linear', ...}
%
%   The trained controller shipped in data/ is a 1-10-10-1 tanh network with
%   zero biases, so Phi maps the scalar output y = C*x to the scalar input u.
%
%   See also LOAD_NN, SECTOR_BOUND_FFNN.
a = y(:);
for L = 1:numel(net.W)
    z = net.W{L}*a + net.b{L};
    switch lower(net.act{L})
        case 'tanh',   a = tanh(z);
        case 'relu',   a = max(0, z);
        case 'linear', a = z;
        otherwise
            error('nn_forward:activation', 'Unknown activation "%s".', net.act{L});
    end
end
u = a;   % scalar for the shipped network
end
