function [ Xupdate ] = miniL1L2( X, Y, W, idx, rho )
%miniL1L2 This function is used to calculate the X^(k+1) at new iterate when 
%           g(X^k) = 0 using a closed form solution.
%Input    X:    X^k at last iterate
%         Y:    Input matrix, including observations
%         W:    Real weight matrix
%         idx:  Index Hankel matrix
%         rho:  Penalty parameter
%Output   Xupdate:  X^(k+1)

[L, K] = size(X);
N = L + K - 1;

xk = [X(1:L, 1); transpose(X(L, 2:end))];
A = Y - X;
obj = A .* W;
sumH = accumarray(idx(:), obj(:));
sumW = accumarray(idx(:), W(:));
x = zeros(N, 1);
for i = 1 : N
    
    if (real(sumH(i, 1)) == 0) && (imag(sumH(i, 1)) ==0)
        
        x(i) = xk(i);
        
    elseif real(sumH(i, 1)) == 0

        if imag(sumH(i, 1)) > rho * i
            x(i) = xk(i) + 1j * (imag(sumH(i, 1))/sumW(i, 1) - rho * i/sumW(i, 1));
        elseif imag(sumH(i, 1)) < -rho * i
            x(i) = xk(i) + 1j * (imag(sumH(i, 1))/sumW(i, 1) + rho * i/sumW(i, 1));
        else
            x(i) = xk(i);
        end
        
    elseif imag(sumH(i, 1)) ==0

        if real(sumH(i, 1)) > rho * i
            x(i) = xk(i) + real(sumH(i, 1))/sumW(i, 1) - rho * i/sumW(i, 1);
        elseif real(sumH(i, 1)) < -rho * i
            x(i) = xk(i) + real(sumH(i, 1))/sumW(i, 1) + rho * i/sumW(i, 1);
        else
            x(i) = xk(i);
        end
        
    else
        if real(sumH(i,1)) > rho*i/sqrt(1 + (imag(sumH(i,1)/real(sumH(i, 1))))^2)
            x(i) = xk(i) + real(sumH(i,1))/sumW(i,1) - rho*i/(sumW(i,1)* sqrt(1 + (imag(sumH(i,1)/real(sumH(i, 1))))^2));
        elseif real(sumH(i,1)) < -rho*i/sqrt(1 + (imag(sumH(i,1)/real(sumH(i, 1))))^2)
            x(i) = xk(i) + real(sumH(i,1))/sumW(i,1) + rho*i/(sumW(i,1)* sqrt(1 + (imag(sumH(i,1)/real(sumH(i, 1))))^2));
        else
            x(i) = xk(i);
        end
        
        if imag(sumH(i,1)) > rho*i/sqrt(1 + (real(sumH(i,1)/imag(sumH(i, 1))))^2)
            x(i) = x(i) + 1j * (imag(sumH(i,1))/sumW(i,1) - rho*i/(sumW(i, 1) * sqrt(1 + (real(sumH(i,1)/imag(sumH(i, 1))))^2)));
        elseif imag(sumH(i,1)) <- rho*i/sqrt(1 + (real(sumH(i,1)/imag(sumH(i, 1))))^2)
            x(i) = x(i) + 1j * (imag(sumH(i,1))/sumW(i,1) + rho*i/(sumW(i, 1) * sqrt(1 + (real(sumH(i,1)/imag(sumH(i, 1))))^2)));
        else
            x(i) = x(i);
        end
    end

end

Xupdate = hankel(x(1:L), x(L:end));
end

