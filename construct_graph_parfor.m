function [G,D,W_sparse,usedtime]=construct_graph_parfor(X,sub_cluster_n,position_2D,par,S)
% ---------------------------------------------
%% Input:
% X=[X1,X2,...,XS] -d*N
% sub_cluster_n=[n1,n2,n3,...nC]; N=\sum_{i=1}^{C}n_i
%% par.sigma (default 1.0--same with svm classifier)
%% par.k     (default 4, 4 nearest neighbors, according to pixel position) i-par.k:i+par.k~j-par.k:j+par.k
%% par.D1    ([D1,D2]=size(original map));
%% output:
% G 	- 	N*N (sparse graph)
%---------------------------------------------
%written by Shujun Yang (yangsj3@sustech.edu.cn; sjyang8-c@my.cityu.edu.hk)
%---------------------------------------------
[d,N]=size(X);
X=X./max(X(:));

W_sparse=sparse(N,N);
D=sparse(N,N);
G=sparse(N,N);
tic
position_2D_all=[];
tic
pos_superpixel_id=[];
for i=1:S
	[d1,d2]=size(position_2D{i});
	num_sub_pixels=max(d1,d2);
    pos_superpixel_id=[pos_superpixel_id; i*ones(num_sub_pixels,1)];
	position_2D_all=[position_2D_all; position_2D{i}];
	
end
[num_pixels,d2]=size(position_2D_all);
W_sparse_cell=cell(N,1);

parfor i=1:num_pixels%% the actual id of 1st pixel in laplacian graph G
	pos1=position_2D_all(i);
	ith=pos_superpixel_id(i);
	[x,y]=pos1d_to_2d(pos1,par.D1);
	%% k nearest neighbors
	x_min=max(1,x-par.k);
	x_max=min(x+par.k,par.D1);
	y_min=max(1,y-par.k);
	y_max=min(y+par.k,par.D2);
    W_sparse_tmp=sparse(1,N);
    [W_sparse(i,:)]=inner_loop(X,i,x_min,x_max,y_min,y_max,par,position_2D,ith,sub_cluster_n,W_sparse_tmp);
end
W_sparse=(W_sparse+W_sparse')/2;

parfor i=1:N;
	Ds_diag(i)=sum(W_sparse(i,:));
end;
D=diag(Ds_diag);
G=D-W_sparse;
usedtime=toc
fprintf(1,'the laplacian graph computes using %d seconds\n',usedtime);
delete(gcp('nocreate'));
end
function [W_sparse_tmp]=inner_loop(X,actual_id1_G,x_min,x_max,y_min,y_max,par,position_2D,ith,sub_cluster_n,W_sparse_tmp)
	for k1=x_min:x_max
			for k2=y_min:y_max
				[pos_1d]=pos2d_to_1d(k1,k2,par.D1);
				[IDx]=find(position_2D{ith}==pos_1d);
				[d1,d2]=size(IDx);
				if min(d1,d2)>0
					actual_id2_G=sum(sub_cluster_n(1:ith-1))+IDx;
					tmp_vector=(X(:,actual_id1_G)-X(:,actual_id2_G));
					W_sparse_tmp(1,actual_id2_G)=exp(-par.sigma*(tmp_vector'*tmp_vector));
				end
			end
	end
end