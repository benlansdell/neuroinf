function [nll,g] = SimultaneousLogisticLoss(W,X,Y)
	[n,p] = size(X);
	nTasks = size(Y,2);
	W = reshape(W,p,nTasks);
	g = zeros(p,nTasks);
	nll = 0;
	for t = 1:nTasks
	    [nll_sub,g(:,t)] = PoissonLoss(W(:,t),X,Y(:,t));
	    nll = nll + nll_sub;
	end
	g = g(:);
end