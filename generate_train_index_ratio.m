function [] = generate_train_index_ratio(per_ratio,img_name)

time_start=tic;
addpath(genpath('.\classification_code\'));
R=load(['.\data\' img_name '_corrected.mat']);
GT_map=importdata(['.\data\' img_name '_gt.mat']);
C=max(unique(GT_map));

random_times=10;

C_total=[];
for class=1:C
    no_class(class)=length(find(GT_map==class));
    C_total=[C_total  no_class(class)];
end
CTrain=ceil(C_total*per_ratio);
random_iter=10;
save_path=['./train_indexes2/' img_name  '/'];
if ~exist(save_path,'dir')
    mkdir(save_path);
end

for i=1:random_iter
    [loc_train, loc_test, CTest] = Generating_training_testing(GT_map,CTrain);
     save([save_path '/' img_name '_train_test' num2str(i) '_' num2str(per_ratio) '.mat'], 'CTrain','loc_train', 'loc_test', 'CTest');
end
end

