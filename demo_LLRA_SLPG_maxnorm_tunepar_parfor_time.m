function[ave_OA_SVM2 ave_AA_SVM2 ave_Kappa_SVM2 ave_TPR_SVM2 mean_time]=demo_LLRA_SLPG_maxnorm_tunepar_parfor_time(img_name,par,num_Pixel,per_ratio,random_iters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the main call funtion for 3-step algorithm framework:
% (1) superpixel segmentation and (2)  low-rank representation (3) SVM classification

%written by Shujun Yang (yangsj3@sustech.edu.cn; sjyang8-c@my.cityu.edu.hk)

generate_train_index_ratio(per_ratio,img_name); %% generate data

time_start=tic;
addpath(genpath('.\classification_code\'));
addpath('./common/');
addpath(genpath(cd));
addpath('./Entropy Rate Superpixel Segmentation/');
%% prepare for data 
%% read file
R=importdata([img_name '_corrected.mat']);
gt=importdata([img_name '_gt.mat']);
R=double(R);
[m,n,d]=size(R);
indexes_v=[1:m*n];
gt_v=reshape(gt,[m*n 1]);
zeros_index=find(gt_v==0);
gt_v(zeros_index)=[];
indexes_v(zeros_index)=[];
GroundT(1:length(indexes_v),1)=indexes_v;
GroundT(1:length(indexes_v),2)=gt_v;

C=max(gt(:));


flag_select=1; %%set as 1: select the build graph with superpixel
lambda=par.lambda;
beta = par.beta;
par.sigma=1.0;
par.D1=m; %% the height of the map
par.D2=n;

%% step1-- over segmentation based on ERS
%% this part of code could reference the superPCA
data3D=R; 
%% max normalize before segment
data3D = data3D./max(data3D(:));

%% super-pixels segmentation
labels = cubseg(data3D,num_Pixel);
labels=labels+1;

%% cut into multiple sub images according to segment labels
[sub_I,size_subimage,position_2D]=Partition_into_subimages_superpixel(data3D,labels);%% input is maxnormlized feature
%% parepare for data matrix
X=[];
for cur=1:num_Pixel
    X=[X;sub_I{cur}];
    sub_cluster_n(cur)=size(sub_I{cur},1);
end
%% save path setting
if flag_select==1
	save_path=[img_name 'SP' num2str(num_Pixel) '_results_feature_LLRA_SLPG_maxnorm_parfor_time/' ];
	save_path2=[img_name 'SP' num2str(num_Pixel) 'SVM_results_LLRA_SLPG_maxnorm_parfor_time/' 'per_C' num2str(per_ratio) '/' 'lambda' num2str(lambda) 'beta' num2str(beta) 'k' num2str(par.k) '/' ];
end
if ~exist(save_path,'dir')
    mkdir(save_path);
end

if ~exist(save_path2,'dir')
    mkdir(save_path2);
end
time_end1=toc(time_start);

%% step2--Graph regularized Low rank model for feature extraction
 %% generate the clean part for image
%% parameter setting
maxiter=10^4;
rho=1.1;
maxtau = 10^12;
para.maxiter=maxiter;
para.lambda=lambda;%  sparse error term coffecient
para.beta=beta; %
para.rho=rho;
para.maxtau=maxtau;
para.DEBUG=1;
L=[];E=[];
if flag_select==1
	res_file_name1= [img_name 'SP' num2str(num_Pixel) 'lambda' num2str(lambda) 'beta' num2str(beta) 'k' num2str(par.k) '_maxnorm_parfor_time.mat'];
	save_feature_file=[img_name 'SP' num2str(num_Pixel) 'sigma' num2str(par.sigma) 'k' num2str(par.k) '_maxnorm_parfor_time.mat'];
end


if ~exist([save_path save_feature_file])
	if flag_select==1
		[G1,D1,W_sparse1,usedtime]=construct_graph_parfor(X',sub_cluster_n,position_2D,par,num_Pixel);
		save([save_path save_feature_file], 'G1','D1','W_sparse1','usedtime');
	end
else
	load([save_path save_feature_file]);
end
time_end2=toc(time_start);

if ~exist([save_path res_file_name1])
	if flag_select==1
		[L,E,Li,Ei,Xi,conv_iter,obj2]=LLRA_SLPG(X',sub_cluster_n,G1,para,num_Pixel);
	end
    
	%% remap the low rank L to the original postition, then we get the low ranked 2D data map L_spe_m for the original data

    L_trans=L';
    L_spe_m=zeros(size(L_trans));
    for cur=1:num_Pixel
        L_spe_m(position_2D{cur},:)=L_trans((sum(sub_cluster_n(1:cur-1))+1:sum(sub_cluster_n(1:cur))),:);
    end
	%% save feature file
	if flag_select==1
		save([save_path res_file_name1],'Li','Ei','conv_iter','obj2','L','E','L_spe_m','par','para','G1','D1','W_sparse1','time_end1','time_end2');
	end
else
    load([save_path res_file_name1]);

end
time_end3=toc(time_start);
L_spe_m=reshape(L_spe_m,[m n d]);

%% step3--we run SVM for per_ratio data

[res,accracy_SVM1,TPR_SVM1,Kappa_SVM1,accracy_SVM2,TPR_SVM2,Kappa_SVM2,...
Predict_SVM1,Predict_SVM2,time_result] = my_Classification_V2_CK_ratio_multiple_iters_time(R,L_spe_m,gt,random_iters,img_name,per_ratio);
if flag_select==1
	res_file_name2= [img_name 'lambda' num2str(lambda) 'beta' num2str(beta) 'k' num2str(par.k) 'per_C' num2str(per_ratio) '_maxnorm_parfor_time.mat'];	
end	
time_end4=toc(time_start);
																													  
save([save_path2 res_file_name2],'res','accracy_SVM1','TPR_SVM1','Kappa_SVM1','accracy_SVM2','TPR_SVM2','Kappa_SVM2',...
'Predict_SVM1','Predict_SVM2','time_end1','time_end2','time_end3','time_result','time_end4');
ave_OA_SVM2=res.ave_OA_SVM2;
ave_AA_SVM2=res.ave_AA_SVM2;
ave_Kappa_SVM2=res.ave_Kappa_SVM2;
ave_TPR_SVM2=res.ave_TPR_SVM2;

mean_time=time_end3+time_result.t1+mean(time_result.t2)+time_result.t3;
end