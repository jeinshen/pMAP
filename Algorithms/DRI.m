function [x, it, Xerr, W] = DRI(y, pars )
        r = pars.r;
        N = pars.N;
        L = pars.L - 1;
        ftol = pars.tol;
		maxit = pars.maxit;
		mu=0.1;				%parameter. Must be in ]0,2[
        %mu = 0;
		gamma=0.51*mu;		%parameter. Must be in ]mu/2,1[
        %gamma = 0;
		%note: setting mu=0 and gamma=0 yields the Douglas-Rachford method
		Tnoisy=toeplitz(y(L+1:N),y(L+1:-1:1)); 
        W = pars.W;
        w = [W(1:L+1, 1); W(L+1,2:end)'];
        W = toeplitz(w(L+1:N),w(L+1:-1:1));
        %W = (fliplr(W));
        %W=ones(N-L,L+1)*(L+1);	%matrix of weights for the weighted Frobenius norm
		%for indexcol=2:L+1
		%	for indexrow=1:L+2-indexcol
		%		W(indexrow,indexcol-1+indexrow)=L+2-indexcol;
		%		W(N-L-indexrow+1,L+3-indexcol-indexrow)=L+2-indexcol;
		%	end
		%end
		Tdenoised=Tnoisy;		%the noisy matrix is the initial estimate
		mats=Tdenoised;			%auxiliary matrix
        it = 1;
        Xerror = 1;
       
        while (it < maxit) && (Xerror > ftol)
			[U, S, V]=svd(mats+gamma*(Tdenoised-mats)+mu*(Tnoisy-Tdenoised).*W,0);
			Tdenoisedr=U(:,1:r)*S(1:r,1:r)*(V(:,1:r))';	%SVD truncation -> Tdenoised has rank K
            %[U, S, V]=svds(mats+gamma*(Tdenoised-mats)+mu*(Tnoisy-Tdenoised).*W,r);
			%Tdenoisedr=U*S*V';	%SVD truncation -> Tdenoised has rank K
			mats=mats-Tdenoisedr+Toeplitzation(2*Tdenoisedr-mats); 
            it = it + 1;
            Xerror = norm(Tdenoised - Tdenoisedr, 'fro');
            Tdenoised = Tdenoisedr;
            Xerr(it, 1) = Xerror;
        end	
        
		%at this point, Tdenoised has rank K but is not exactly Toeplitz
		Tdenoised=Toeplitzation(Tdenoised);
        x = [transpose(Tdenoised(1, (L+1):-1:2));Tdenoised(:, 1)];
end

