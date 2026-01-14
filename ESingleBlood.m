function E=ESingleBlood(x,k_1,k_2,eta,A,beta,a_1,a_0)

gamma=x-(a_1+a_0)/2;

term1=k_1./(1+k_1*gamma);


n_star=eta+k_2*eta/k_1*(1+k_1*gamma)./(1-exp(-k_2*(A-x)));

E=term1.*(n_star.*(1-beta*gamma)-eta*(1+beta/k_1))-eta*k_2;