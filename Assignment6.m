clc
clear
Y = [607, 583, 521, 494, 369, 782, 570, 678, 467, 620, 425, 395,...
    346, 361, 310, 300, 382, 294, 315, 323, 421, 339, 398, 328, 335,...
    291, 329, 310, 294, 321, 286, 349, 279, 268, 293, 310, 259, 241,...
    243, 272, 247, 275, 220, 245, 268, 357, 273, 301, 322, 276, 401,...
    368, 149, 507, 411, 362, 358, 355, 362, 324, 332, 268, 259, 274,...
    248, 254, 242, 286, 276, 237, 259, 251, 239, 247, 260, 237, 206,...
    242, 361, 267, 245, 331, 357, 284, 263, 244, 317, 225, 254, 253,...
    251, 314, 239, 248, 250, 200, 256, 233, 427, 391, 331, 395, 337,...
    392, 352, 381, 330, 368, 381, 316, 335, 316, 302, 375, 361, 330,...
    351, 186, 221, 278, 244, 218, 126, 269, 238, 194, 384, 154, 555,...
    387, 317, 365, 357, 390, 320, 316, 297, 354, 266, 279, 327, 285,...
    258, 267, 226, 237, 264, 510, 490, 458, 425, 522, 927, 555, 550,...
    516, 548, 560, 545, 633, 496, 498, 223, 222, 309, 244, 207, 258,...
    255, 281, 258, 226, 257, 263, 266, 238, 249, 340, 247, 216, 241,...
    239, 226, 273, 235, 251, 290, 473, 416, 451, 475, 406, 349, 401,...
    334, 446, 401, 252, 266, 210, 228, 250, 265, 236, 289, 244, 327,...
    274, 223, 327, 307, 338, 345, 381, 369, 445, 296, 303, 326, 321,...
    309, 307, 319, 288, 299, 284, 278, 310, 282, 275, 372, 295, 306,...
    303, 285, 316, 294, 284, 324, 264, 278, 369, 254, 306, 237, 439,...
    287, 285, 261, 299, 311, 265, 292, 282, 271, 268, 270, 259, 269,...
    249, 261, 425, 291, 291, 441, 222, 347, 244, 232, 272, 264, 190,...
    219, 317, 232, 256, 185, 210, 213, 202, 226, 250, 238, 252, 233,...
    221, 220, 287, 267, 264, 273, 304, 294, 236, 200, 219, 276, 287,...
    365, 438, 420, 396, 359, 405, 397, 383, 360, 387, 429, 358, 459,...
    371, 368, 452, 358, 371];
ind = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 5, 5,...
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6,...
    7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10,...
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11,...
    11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12,...
    12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,...
    12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14,...
    14, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 15, 16, 16, 16, 16,...
    16, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 20,...
    20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,...
    20, 20, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 22, 22, 22, 22, 22,...
    22, 22, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,...
    24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,...
    24, 24, 24, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 26,...
    26, 26, 26, 26, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28, 28, 28, 28,...
    28, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30,...
    31, 31, 31, 31, 31, 32, 32, 32, 32, 32, 33, 34, 34, 34, 34, 34, 34,...
    34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34];
child = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0,...
    0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0];
zy = log(Y); %log of time data for numerical stability
n = length(unique(ind)); %calculating number of individuals
ind_mean = zeros(1,n);
for i = 1:n
    ind_mean(i) = mean(zy(ind==i)); %calculating individuals' means
end
phi = 1;
x0_a6 = [mean(zy),std(zy),std(zy),phi,ind_mean]; %[mean,sigma,tau,ind means]
logpdf_a6 = @(x)logpost_a6(x,zy,ind,child); %function for posterior calculation
N = 1000000; % no of samples
%using built-in slicesampler for MCMC sample generation 
samples_a6 = slicesample(x0_a6,N,'logpdf',logpdf_a6,'burnin',1000); 
x0_a5 = [mean(zy),std(zy),std(zy),ind_mean]; %[mean,sigma,tau,ind means]
logpdf_a5 = @(x)logpost_a5(x,zy,ind); %function for posterior calculation
%using built-in slicesampler for MCMC sample generation 
samples_a5 = slicesample(x0_a5,N,'logpdf',logpdf_a5,'burnin',1000); 

% for assignment 6, i have used a6 with parameters
mu_a6 = samples_a6(:,1);
sigma_a6 = samples_a6(:,2);
tau_a6 = samples_a6(:,3);
phi_a6 = samples_a6(:,4);

% for assignment 5, i have used a5 with parameters
mu_a5 = samples_a5(:,1);
sigma_a5 = samples_a5(:,2);
tau_a5 = samples_a5(:,3);

% Task 1- modelling the effect of different reaction times
plothistogram(phi_a6); title('Phi Original Scale (logy)')

% Task 2- Posterior of tau for Assignment 5 and 6
figure; 
subplot(1,2,1); plothistogram(tau_a5); title('Tau (Assignment 5)')
subplot(1,2,2); plothistogram(tau_a6); title('Tau (Assignment 6)')

% Task 3- Plot two prior distributions, one for kids and one for adults
prior_child = mean(mu_a6 + phi_a6) + mean(tau_a6).*randn(1,N);
prior_adult = mean(mu_a6) + mean(tau_a6).*randn(1,N);
% Compare against a single prior of theta in Assignment 5.
assign_5 = mean(mu_a5) + mean(tau_a5).*randn(1,N);
figure; 
subplot(1,3,1); plothistogram(prior_child); title('Prior Distribution for Child')
subplot(1,3,2); plothistogram(prior_adult); title('Prior Distribution for Adult')
subplot(1,3,3); plothistogram(assign_5); title('Prior Distribution for Assignment 5')

% % Task 4- posterior predictive distribution
posterior_child = mean(mu_a6 + phi_a6) + mean(tau_a6).*randn(1,N) + mean(sigma_a6).*randn(1,N);
posterior_adult = mean(mu_a6) + mean(tau_a6).*randn(1,N) + mean(sigma_a6).*randn(1,N);
% for a bernoulli process, we have beta distribution for posterior predictive distribution
% posterior is given as beta(z+1, N-z+1). z is number of kids in total individuals N.
unknown = betarnd(sum(child)+1,length(child)-sum(child)+1);
posterior_unknown = mean(mu_a6 + phi_a6.*unknown) + mean(tau_a6).*randn(1,N) + mean(sigma_a6).*randn(1,N);
figure; 
subplot(1,3,1); plothistogram(posterior_child); title('Posterior Distribution for Child')
subplot(1,3,2); plothistogram(posterior_adult); title('Posterior Distribution for Adult')
subplot(1,3,3); plothistogram(posterior_unknown); title('Posterior Distribution for Unknown')

%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%

function logposterior = logpost_a6(variables,timedata,ind,child)
% timedata = log of rection times
mu = variables(1);
sigma = variables(2);
tau = variables(3);
phi = variables(4);
theta = variables(5:end);
n = length(theta);
logprior_mu = 0;
logprior_sigma = log(double(sigma>0));
logprior_tau = log(double(tau>0));
logprior_phi = 0;
logprior_theta = 0;
for i = 1:n
    logprior_theta = logprior_theta - (theta(i)-(mu+(phi*child(i))))^2/2/tau^2 - log(tau);
end
loglikelihood = 0;
for i = 1:n
    z = timedata(ind==i);
    n2 = length(z);
    for j = 1:n2
        loglikelihood =  loglikelihood - (z(j)-theta(i))^2/2/sigma^2 - log(sigma);
    end
end
logposterior = logprior_mu + logprior_tau + logprior_sigma + logprior_theta + logprior_phi + loglikelihood;
end

function logposterior = logpost_a5(variables,timedata,ind)
mu = variables(1);
sigma = variables(2);
tau = variables(3);
theta = variables(4:end);
n = length(theta);
logprior_mu = 0;
logprior_sigma = log(double(sigma>0));
logprior_tau = log(double(tau>0));
logprior_theta = 0;
for i = 1:n
    logprior_theta = logprior_theta - (theta(i)-mu)^2/2/tau^2 - log(tau);
end
loglikelihood = 0;
for i = 1:n
    z = timedata(ind==i);
    n2 = length(z);
    for j = 1:n2
        loglikelihood =  loglikelihood - (z(j)-theta(i))^2/2/sigma^2 - log(sigma);
    end
end
logposterior = logprior_mu + logprior_tau + logprior_sigma + logprior_theta + loglikelihood;
end

function HDIofMCMC = HDI(sampleVec,credMass)
sortedPts = sort(sampleVec);
ciIdxInc = ceil(credMass*length(sortedPts));
nCIs = length(sortedPts) - ciIdxInc; %number of intervals
ciWidth = zeros(1, nCIs); %width of interval
for i = 1:nCIs
    ciWidth(i) = sortedPts(i + ciIdxInc) - sortedPts(i);
end
[~,id] = min(ciWidth);
HDImin = sortedPts(id);
HDImax = sortedPts(id + ciIdxInc);
HDIlim = [HDImin, HDImax];
HDIofMCMC = HDIlim;
end

function plot = plothistogram(Parameter)
parameter_mean = mean(Parameter);
[M,edges] = histcounts(Parameter);
[~, id] = max(M);
parameter_mode = edges(id);
d = sort(Parameter);
parameter_median = d(end/2);
phi_CI = HDI(Parameter,0.95);
histogram(Parameter,'normalization','pdf','EdgeColor','none'); 
hold on; xl = xlim(); yl = ylim(); xWidth = xl(2)-xl(1); yHeight = yl(2)-yl(1);
x = (xl(1) + 0.8 * xWidth); y = (yl(1) + 0.8 * yHeight);
data = sprintf('Mode = %g \n%s %g',parameter_mode,'Median = ',parameter_median);
text(x,y,data); 
lower = phi_CI(1,1); upper = phi_CI(1,2);
yL = get(gca,'YLim');
line([parameter_mean parameter_mean],yL,'LineWidth',2); 
line([lower lower],yL,'LineWidth',2); 
line([upper upper],yL,'LineWidth',2); 
text(parameter_mean,y*1.2,sprintf('Mean = %g',parameter_mean),'HorizontalAlignment','center')
text(lower,y*1.1,sprintf('HDI low = %g',lower),'HorizontalAlignment','right')
text(upper,y*1.1,sprintf('HDI up = %g',upper),'HorizontalAlignment','left')
set(findall(gcf,'-property','FontSize'),'FontSize',16); hold off;
end
