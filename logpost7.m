function logposterior = logpost7(variables,timedata,ind,child,attempt)
% timedata = log of rection times
theta = variables(8:end);
a = length(theta)/2;
sigma = variables(1);
mu_0 = variables(2);
tau_0 = variables(3);
phi_0 = variables(4);
theta_0 = theta(1:a);
mu_1 = variables(5);
tau_1 = variables(6);
phi_1 = variables(7);
theta_1 = theta(a+1:end);
theta_0=theta_0(:);
theta_1=theta_1(:);
child=child(:);
attempt=attempt(:);
logprior_theta_0 = sum(-log(tau_0)-0.5*((theta_0-(mu_0+phi_0.*child))./tau_0).^2);
logprior_theta_1 = sum(-log(tau_1)-0.5*((theta_1-(mu_1+phi_1.*child))./tau_1).^2);
logprior_sigma = log(double(sigma>0));
logprior_mu_0 = 0;
logprior_tau_0 = log(double(tau_0>0));
logprior_phi_0 = 0;
logprior_mu_1 = 0;
logprior_tau_1 = log(double(tau_1>0));
logprior_phi_1 = 0;
mu=theta_0(ind)+attempt.*theta_1(ind);
loglikelihoodl=-length(timedata)*log(sigma)-0.5*sum((timedata(:)-mu).^2)/sigma^2;
logprior_0=logprior_theta_0+logprior_tau_0+logprior_mu_0+logprior_phi_0;
logprior_1=logprior_theta_1+logprior_tau_1+logprior_mu_1+logprior_phi_1;
if (sigma>0)&&(tau_0>0)&&(tau_1>0)
    logposterior=loglikelihoodl+logprior_0+logprior_1+logprior_sigma;
    else
    logposterior=-inf;
end   
end
