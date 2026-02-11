% A Column Generation Approach to Exact Experimental Design
% Authors: Selin Ahipasaoglu, Stefano Cipolla, Jacek Gondzio
% arXiv ID: 2507.03210 (July 2025)

clear all;
close all;
clc;
%The path on which all the data are:
problems_path = '/scratch/sc9c23/Dataset/Gen1';%'/scratch/sc9c23/large_D_Optimal';   %'./Dataset/Gen1'; 
% Kurtosis Parameters
kur_params = logspace(0, 2, 5);  
%Finds all the problems and stores their names in a struct
d = dir(fullfile(problems_path,'Generator_10M_*'));
% N
NN = [0];
expN = length(d)*length(kur_params)*length(NN);
%File to Sotre Results
resultsIBOSS = table( strings(expN,1), NaN(expN,1), NaN(expN,1), NaN(expN,1),NaN(expN,1),NaN(expN,1),...
                                       NaN(expN,1),...
                                     'VariableNames', {'Problem', 'n', 'm','Time','Kurtosis','N','Obj'});

exp_num =0;
tol                                 = 1e-5;
init_type                       = 'KY';
for k = 1:length(d) 
    load(fullfile(problems_path,d(k).name));
    
    for jj =1:length(kur_params)
         kur_param = kur_params(jj);
         X                            = genmatrix;
         d(k).name
         [n, m ]                      = size(X);
         X                            = sinh((1 / kur_param) * asinh(X));
         computed_kur                 = kurtosis_summary(X);  
         for ll = 1:length(NN)
            


             ll
             N = 2*n+ NN(ll);
             exp_num                  = exp_num +1;

             time_IBOSS               = 0;
             tic
             fit                      = iboss_od_mex(X.', ones(m,1), N);
             time_IBOSS               = time_IBOSS + toc;
             X_S                      = X(:,fit.index);
             obj_IBOSS                = log(det(X_S*X_S.'));
             
             % Record Performance ColGen for Exact
             resultsIBOSS.Problem(exp_num)                 = d(k).name;
             resultsIBOSS.n(exp_num )                      = n;
             resultsIBOSS.m(exp_num )                      = m;
             resultsIBOSS.Time(exp_num )                   = time_IBOSS;
             resultsIBOSS.N(exp_num )                      = N ;
             resultsIBOSS.Obj(exp_num )                    = obj_IBOSS;
             resultsIBOSS.Kurtosis(exp_num )               = computed_kur ;
         end   
          
    end
end

writetable(resultsIBOSS, '../Results/IBOSS_10M_mex.csv');