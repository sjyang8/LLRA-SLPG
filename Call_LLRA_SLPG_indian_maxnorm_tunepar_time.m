function []=Call_LLRA_SLPG_indian_maxnorm_tunepar_time()
%the main call funtion for running the whole algorithm
%written by Shujun Yang (yangsj3@sustech.edu.cn; sjyang8-c@my.cityu.edu.hk)

%% parameter setting
per_ratio=0.05; % The training percentage per class equals 5%
img_name='Indian_pines'; % data set name
num_Pixel=64;  % number of superpixels
par.lambda=0.1;% parameter of \lambda
par.beta=50;   % parameter of \beta
par.k=1;       % parameter for radius r
random_iters=10; %% random 10 times for SVM classifier

[ave_OA_SVM2 ave_AA_SVM2 ave_Kappa_SVM2 ave_TPR_SVM2 mean_time]=demo_LLRA_SLPG_maxnorm_tunepar_parfor_time(img_name,par,num_Pixel,per_ratio,random_iters);
 disp(['average performance for ' num2str(random_iters) ' times: OA(mean)=' num2str(ave_OA_SVM2)  ',AA(mean)=' num2str(ave_AA_SVM2) ...
',Kappa(mean)=' num2str(ave_Kappa_SVM2)]);
disp(['average time for ' num2str(random_iters) ' times (seconds):' num2str(mean_time)]);
end