function [Pn, p0] = MGmm_SD_Queue(k,Kmax,lambda,s,Rc)
%%Inputs:
% k: Number of Users
% Kmax: Max Number of Users
% lambda: Users Arriaval Rate
% s: File size for single user
% Rc: Bit Rate for different user states

if all(Rc)
    fn = Rc./Rc(1);
    CumProduct = cumprod(fn);
    
    % p0 calculations
    DD = zeros(1,Kmax);
    for n = 1:1:Kmax
        DD(n) = ((lambda*s/Rc(1))^n)/(factorial(n)*CumProduct(n));
    end
    p0 = 1/(1 + sum(DD));
    
    Pn = (((lambda*(s/Rc(1)))^k)/(factorial(k)* CumProduct(k)))*p0;
    
else
    error('Error: Rc contains zero vaule');
end

end
