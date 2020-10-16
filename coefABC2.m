function coef = coefABC2(S,order)
    if order==1
        coef(1,2) = (S-1)/(S+1);
        coef(1,1) = coef(2,1); % (start:end)
    else
        coef(1:3,2) = [-(1./S-2+S);
                         -2*(S-1./S);
                          4*(S+1./S)]./(1./S+2+S);
        coef(1:3,1) = coef(1:3,2); % (start:end),(0,+/-1,+/-2)
    end
end