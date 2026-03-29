function c = kkb2v3(x, f, w, n)
%KKB2 KKB met veeltermen van maximale graad $n$
%   Stelt discrete kleinste-kwadratenveeltermbenadering op
%   van graad n voor de punten $(x_i,f_i)$ met gewichten $w_i$,
%   $i = 1,\ldots,N$.
    x = x(:); f = f(:); w = w(:);
    X = bsxfun(@power, x, 0:2*n);
    
    A = w'*X;
    A = hankel(A(1:n+1), A(n+1:end)');
    b = (w.*f)'*X(:,1:n+1);
    c = A \ b';
end
