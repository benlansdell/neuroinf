function h = hessian(X,y,b_hat)
	nT = size(X,1);
	h = -X'*spdiags(exp(X*b_hat), 0, nT, nT)*X;
end