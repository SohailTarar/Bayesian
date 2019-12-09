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
zy = log(Y); %log of time data for numerical stability
n = length(unique(ind)); %calculating number of individuals
ind_mean = zeros(1,n);
for i = 1:n
    ind_mean(i) = mean(zy(ind==i)); %calculating individuals' means
end
x0 = [mean(zy),std(zy),std(zy),ind_mean]; %[mean,sigma,tau,ind means]
logpdf = @(x)logpost(x,zy,ind); %function for posterior calculation
N = 1000000; % no of samples
w  = 0.1; %jump
%using built-in slicesampler for MCMC sample generation 
samples = slicesample(x0,N,'logpdf',logpdf,'width',w); 

% Task A-1: Dude reaction time
mu = samples(:,1);
sigma = samples(:,2);
tau = samples(:,3);
individuals = samples(:,4:end);
individuals_unLog = exp(individuals + 0.5 .* (sigma).^2);
dude = individuals_unLog(:,4);
dude_mean = mean(dude);
[M,edges] = histcounts(dude);
[~, id] = max(M);
dude_mode = edges(id);
d = sort(dude);
dude_median = d(end/2);
dude_CI = HDI(dude,0.95);
figure; h = histogram(dude,'normalization','pdf'); hold on;
xlabel('Reaction Time'); ylabel('pdf'); title('Dude Reaction Time');
xl = xlim(); yl = ylim(); xWidth = xl(2)-xl(1); yHeight = yl(2)-yl(1);
x = (xl(1) + 0.7 * xWidth); y = (yl(1) + 0.8 * yHeight);
data = sprintf('Dude Mean = %g \n%s %g \n%s %g \n%s [%g    %g]',...
    dude_mean,'Dude Mode = ',dude_mode,'Dude Median = ',dude_median,...
    'Dude 95% HDI = ',dude_CI);
text(x,y,data); 
set(findall(gcf,'-property','FontSize'),'FontSize',16)

% Task A-2: Group reaction time using random guy
% i. expected reaction time
group = exp(mu + 0.5 .* sigma.^2 + 0.5 .* tau.^2);
group_mean = mean(group);
[M,edges] = histcounts(group);
[~, id] = max(M);
group_mode = edges(id);
d = sort(group);
group_median = d(end/2);
group_CI = HDI(group,0.95);
figure; histogram(group,'normalization','pdf'); hold on;
xlabel('Reaction Time'); ylabel('pdf'); 
title('Expected Reaction Time for Random Guy');
xl = xlim(); yl = ylim(); xWidth = xl(2)-xl(1); yHeight = yl(2)-yl(1);
x = (xl(1) + 0.7 * xWidth); y = (yl(1) + 0.8 * yHeight);
data = sprintf('Expected Mean = %g \n%s %g \n%s %g \n%s [%g    %g]',...
    group_mean,'Expected Mode = ',group_mode,'Expected Median = ',...
    group_median,'Expected 95% HDI = ',group_CI);
text(x,y,data); 
set(findall(gcf,'-property','FontSize'),'FontSize',16)

% ii. predicted reaction time
% N = 1000000;
Random_Guy = zeros(1,N);
for j = 1:N
    mu_j = samples(j,1);
    sigma_j = samples(j,2);
    tau_j = samples(j,3);
    mean_j = exp(mu_j + 0.5 .* sigma_j.^2 + 0.5 .* tau_j.^2);
    variance_j = exp(2*mu_j + (sigma_j).^2)*(exp(sigma_j.^2)-1);
    std_deviation_j = sqrt(variance_j);
    Random_Guy(j) = std_deviation_j.*randn(1) + mean_j; %generating a 
    %random measurement using mean and standard deviation
end
random_mean = mean(Random_Guy);
[M,edges] = histcounts(Random_Guy);
[~, id] = max(M);
random_mode = edges(id);
r = sort(Random_Guy);
random_median = r(end/2);
random_CI = HDI(Random_Guy,0.95);   
figure; histogram(Random_Guy,'normalization','pdf');
xlabel('Reaction Time'); ylabel('pdf'); 
title('Predicted Reaction Time for Random Guy');
xl = xlim(); yl = ylim(); xWidth = xl(2)-xl(1); yHeight = yl(2)-yl(1);
x = (xl(1) + 0.7 * xWidth); y = (yl(1) + 0.8 * yHeight);
data = sprintf('Predicted Mean = %g \n%s %g \n%s %g \n%s [%g    %g]',...
    random_mean,'Predicted Mode = ',random_mode,'Predicted Median = ',...
    random_median,'Predicted 95% HDI = ',random_CI);
text(x,y,data); 
set(findall(gcf,'-property','FontSize'),'FontSize',16)

% Task A-3: 
n = length(unique(ind)); %number of individuals
ind_mean_without_hierarchical = zeros(1,n);
for k = 1:n
    ind_mean_without_hierarchical(k) = mean(Y(ind==k));
end
ind_mean_with_hierarchical = mean(individuals_unLog); 
figure;
scatter((1:1:n),ind_mean_without_hierarchical,400,'.','b'); hold on;
scatter((1:1:n),ind_mean_with_hierarchical,300,'.','r');
set( gca, 'XGrid', 'on' );  xticks([0:1:35]);
xlabel('Individual ID'); ylabel('Reaction Time'); 
title('comparison of sample means with & without hierarchical model');
hline = refline(0,group_mean);
hline.Color = 'r'; 
legend('Without Hierarchical','With Hierarchical','Group Mean');
set(findall(gcf,'-property','FontSize'),'FontSize',16)

function logposterior = logpost(variables,timedata,ind)
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
% Computes highest density interval from a sample of representative values,
% estimated as shortest credible interval.
% sampleVec is a vector of representative values from a probability distribution
% credMass is a scalar between 0 and 1, indicating the mass within the credible
% interval that is to be estimated.
% HDIlim is a vector containing the limits of the HDI 
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

