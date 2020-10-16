function coef = coefABC(obj,coefEz_y,coefHy_z,order)
%#codegen

% tmp = Courant/sqrt(mu_r*eps_r)
tmp1 = sqrt(mean(obj.(coefEz_y)(1:2,1))*obj.(coefHy_z)(1,1));
tmp2 = sqrt(mean(obj.(coefEz_y)(end-(0:1),1))*obj.(coefHy_z)(end,1));

if order==1
	coef = [(tmp1-1)./(tmp1+1),(tmp2-1)./(tmp2+1)];
else
	coef = [[-(1/tmp1-2+tmp1);-2*(tmp1-1/tmp1);4*(tmp1+1/tmp1)]/(1/tmp1+2+tmp1), ...
			[-(1/tmp2-2+tmp2);-2*(tmp2-1/tmp2);4*(tmp2+1/tmp2)]/(1/tmp2+2+tmp2)];
end

end
