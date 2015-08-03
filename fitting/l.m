function ll = l(X,y,b_hat,r,lambda)
	if size(b_hat,1) == 1
		b_hat = b_hat';
	end
	if size(y,1) ~= 1
		y = y';
	end
	%size(b_hat)
	%size(y)
	%size(X)
	ll = y*X*b_hat-sum(exp(X*b_hat)+log(y+0.00001)');
end