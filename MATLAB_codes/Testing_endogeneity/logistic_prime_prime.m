function fr_prime_prime = logistic_prime_prime(r)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
fr_prime_prime = -logistic_prime(r)*(1+2*exp(-r)*logistic(r));

end

