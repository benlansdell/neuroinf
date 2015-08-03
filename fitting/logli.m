function [ll, d, h] = logli(X,y,b_hat)
	ll = -l(X,y,b_hat);
	d = -dl(X,y,b_hat);
	h = -hessian(X,y,b_hat);
end