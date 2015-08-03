function [ll, d, h] = PoissonLoss(b_hat,X,y)
    ll = -l(X,y,b_hat);
    d = -dl(X,y,b_hat);
    h = -hessian(X,y,b_hat);
end

function ll = l(X,y,b_hat)
    ll = y*X*b_hat-sum(exp(X*b_hat)+log(y+0.00001)');
end

function d = dl(X,y,b_hat)
    d = (y*X)'-X'*exp(X*b_hat);
end

function h = hessian(X,y,b_hat)
    nT = size(X,1);
    h = -X'*spdiags(exp(X*b_hat), 0, nT, nT)*X;
end