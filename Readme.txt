------------------------------------------------------------------------------------------
	                   Readme for the LLRA-SLPG Package
	 		       version Oct 5, 2022
------------------------------------------------------------------------------------------
The package includes the MATLAB code of the LLRA-SLPG algorithm in paper "Local Low-Rank Approximation With Superpixel-Guided Locality Preserving Graph for Hyperspectral Image Classification" [1].

[1] S.-J. Yang, Y. Zhang, Y.-H. Jia, W.-J. Zhang. "Local Low-Rank Approximation With Superpixel-Guided Locality Preserving Graph for Hyperspectral Image Classification." IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing (2022).

1)  Get Started
For LLRA-SLPG,  you can call the "Call_LLRA_SLPG_indian_maxnorm_tunepar_time" function to run the algorithm for Indian_pines dataset. After above procedure, you can see the on-screen instructions for showing the result, such as 
“average performance for 10 times: OA(mean)=*,AA(mean)=*,Kappa(mean)=*". 
“average time for 10 times (seconds): *”.
2). Details

For LLRA-SLPG,  the "Call_LLRA_SLPG_indian_maxnorm_tunepar_time" will automatically create a new folder for saving the detailed results,  i.e. "Indian_pinesSP64SVM_results_LLRA_SLPG_maxnorm_parfor_time/per_C0.05/lambda0.1beta50k1/". In such folder, 
for example, the result file "Indian_pineslambda0.1beta50k1per_C0.05_maxnorm_parfor_time.mat" saves the results.  In the following, we also show the detailed description for the 
variables that stored in the result file.
accracy_SVM1 	-------- 	OA on original feature (10 times)
accracy_SVM2  	-------- 	OA on denoising feature by LLRA_SLPG (10 times)
Kappa_SVM1     	--------	Kappa on original feature (10 times)
Kappa_SVM2    	--------	Kappa on denoising feature by LLRA_SLPG (10 times)
TPR_SVM1         	--------	Accuracy of each class on original feature (10 times)
TPR_SVM2  	--------    Accuracy of each class on denoising feature by LLRA_SLPG (10 times)
Predict_SVM1	--------	The prediction on original feature (10 times)
Predict_SVM2	--------	The prediction on denoising feature by LLRA_SLPG (10 times)
res   -------  the average results for 10 times
res.ave_OA_SVM1 -------- mean OA on original feature
res.ave_OA_SVM2 -------- mean OA on denoising feature by LLRA_SLPG 
res.ave_AA_SVM1 -------- mean AA on original feature
res.ave_AA_SVM2 -------- mean AA on denoising feature by LLRA_SLPG 
res.ave_Kappa_SVM1 -------- mean Kappa on original feature
res.ave_Kappa_SVM2 -------- mean Kappa on denoising feature by LLRA_SLPG 
res.ave_TPR_SVM1 -------- mean Accuracy of each class on original feature
res.ave_TPR_SVM2 -------- mean Accuracy of each class on denoising feature by LLRA_SLPG 
Dependancies:
1)Entropy Rate Superpixel Segmentation
2)MATLAB toolboxes on my PC that you may need:
-----------Deep Learning Toolbox
-----------Image Processing Toolbox
-----------Mapping Toolbox
-----------Optimization Toolbox
-----------Parallel Computing Toolbox
-----------Statistics and Machine Learning Toolbox
-----------Symbolic Math Toolbox

Acknowledgement
1) Thanks to the paper "SuperPCA: A Superpixelwise PCA Approach for Unsupervised Feature Extraction of Hyperspectral Imagery", we refer some codes of it that saved in the folder "common".
2) Thanks to the paper "Simultaneous spatial and spectral low-rank representation of hyperspectral images for classification", we refer the classification codes of it that saved in the folder "classification_code".

ATTN: 
- This package is free for academic usage. You can run it at your own risk. 

- This package was developed by Ms. Shu-Jun Yang (yangsj3@sustech.edu.cn; sjyang8-c@my.cityu.edu.hk). For any problem concerning the code, please feel free to contact Ms. Yang.

------------------------------------------------------------------------------------------
