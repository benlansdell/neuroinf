function d = dl(X,y,b_hat)
	d = (y*X)'-X'*exp(X*b_hat);
end