function [L,E,Li,Ei,Xi,conv_iter,obj2]=LLRA_SLPG(X,sub_cluster_n,G,para,S)
% Solve the Local Low-Rank Approximation with Superpixel-guided Locality Preserving Graph for ...
% Hyperspectral Image Classification
%% min_{L,E} \sum_{i=1}^{S} ||Li||_{*}+\lambda ||Ei||_1+\beta tr(L^{T}GL)
%% s.t. X=L+E, J=L; L=[L1,L2,...,LS];
% ---------------------------------------------
%% Input:
% X=[X1,X2,...,XS] -d*N
% sub_cluster_n=[n1,n2,n3,...nC]; N=\sum_{i=1}^{C}n_i
% G 	- 	N*N 
% para  -	parameter
%			para.maxiter	-   maximum number of iterations
%			para.lambda		- 	parameter for sparse error term (l1 norm)
%			para.beta		-	parameter for local graph term
%			para.rho 		-	rho>=1, ratio used to increase tau
%			para.tau		- 	stepsize for dual variable updating in ADMM
%			para.maxtau		-	maximum stepsize
%			para.DEBUG		-	0 or 1
%% output:
%	L=[L1,L2,...LC];  - (d*N) matrix (d*N: dimension*num_samples)
%	E=[E1,E2,...EC];  - (d*N) matrix (d*N: dimension*num_samples)
%	Li				  -  cell for store the Li{ith}(stands for L_i)
%	Ei				  -  cell for store the Ei{ith}(stands for E_i)
%	Xi				  -  cell for store the Xi{ith}(stands for X_i)
% 	conv_iter		  -  the running iterations until the algorithm converges
%	isfinite_flag	  -  check whether exists abnormal num
%	obj2 			  -  the objective function value
%---------------------------------------------
%written by Shujun Yang (yangsj3@sustech.edu.cn; sjyang8-c@my.cityu.edu.hk)
%---------------------------------------------
%% initialization
maxiter = para.maxiter;
lambda  = para.lambda;
beta    = para.beta;
rho     = para.rho;
maxtau  = para.maxtau;
DEBUG   = para.DEBUG;
tau     = 10^(-4);
tol     = 1e-6;

Y1=zeros(size(X));
Y2=zeros(size(X));

L=[];E=[];J=[];
for ith=1:S
    ith
    Xi{ith}=X(:,(sum(sub_cluster_n(1:ith-1))+1:sum(sub_cluster_n(1:ith))));
    Li{ith}=zeros(size(Xi{ith}));
    Ji{ith}=zeros(size(Xi{ith}));
    Ei{ith}=zeros(size(Xi{ith}));
    Y1i{ith}=zeros(size(Xi{ith}));
    Y2i{ith}=zeros(size(Xi{ith}));
    L=[L Li{ith}];
    E=[E Ei{ith}];
    J=[J Ji{ith}];
end
obj=[];obj2=[];

%% main loop
for iter=1:maxiter
	tic
	newL=[]; newE=[];
	for ith=1:S
		Y1i{ith}=Y1(:,(sum(sub_cluster_n(1:ith-1))+1:sum(sub_cluster_n(1:ith))));
        Y2i{ith}=Y2(:,(sum(sub_cluster_n(1:ith-1))+1:sum(sub_cluster_n(1:ith))));
        % updata Li
		A=Xi{ith}-Ei{ith}+Y1i{ith}/tau;
		B=Ji{ith}+Y2i{ith}/tau;
		Li{ith}=Do(1/(2*tau),(A+B)/2);
		% updata Ei
		Ei{ith}=So(lambda/tau,Xi{ith}-Li{ith}+Y1i{ith}/tau);
		newL=[newL Li{ith}];
        newE=[newE Ei{ith}];
	end
	L=newL;
    E=newE;
	%% update J
	J=(tau*L-Y2)/(beta*(G'+G)+tau*speye(size(G)));
	for ith=1:S
        Ji{ith}=J(:,(sum(sub_cluster_n(1:ith-1))+1:sum(sub_cluster_n(1:ith))));
    end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% - update the auxiliary lagrange Variables
    %% updata Y
	Y1=Y1+tau*(X-L-E);
	Y2=Y2+tau*(J-L);
	%% updata tau
    tau=min(rho*tau,maxtau);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    derr1=X-L-E;
    derr2=L-J;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if DEBUG
		obj(iter)=beta*trace(J*G*J');
		for ith=1:S
			obj(iter)=obj(iter)+norm_nuclear(Li{ith})+lambda*norm(Ei{ith},1);
		end
		obj2(iter)=obj(iter)+tau/2*(norm(X-L-E+Y1/tau,'fro'))^2+tau/2*(norm(J-L+Y2/tau,'fro'))^2;
		if (iter == 1) || (mod(iter, 10) == 0) || ((max(max(abs(derr1))) < tol) && (max(max(abs(derr2)))<tol))
            fprintf(1, 'iter: %d \t err: %f \t  err2: %f \t rank(L): %f \t max(E): %f \t obj: %f \n', ...
                iter,max(max(abs(derr1))),max(max(abs(derr2))),rank(L),max(E(:)),obj2(iter));
        end
	else
		if (iter == 1) || (mod(iter, 10) == 0) || ((max(max(abs(derr1))) < tol) && (max(max(abs(derr2)))<tol))
            fprintf(1, 'iter: %d \t err: %f \t  err2: %f \t rank(L): %f \t max(E): %f \n', ...
                iter,max(max(abs(derr1))),max(max(abs(derr2))),rank(L),max(E(:)));
        end
	end
	if (iter >10 && max(max(abs(derr1)))<tol ) && (max(max(abs(derr2)))<tol )
        current_iter=iter;
        break;
    end
	toc
end
if current_iter<para.maxiter
    conv_iter=current_iter;
else
    conv_iter=para.maxiter;
end
end
function r = So(tau, X)
% shrinkage operator
r = sign(X) .* max(abs(X) - tau, 0);
end

function r = Do(tau, X)
% shrinkage operator for singular values
[U, S, V] = svd(X, 'econ');
r = U*So(tau, S)*V';
end


