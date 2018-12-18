%DIVERS: Biological (non-technical) and technical variance and covariance
%   contribution estimates from absolute bacterial abundance data
%
%   *User required to specify input and saving directories
%
%  INPUTS: 1) Two absolute abundance tables (data_X, data_Y) of
%             equal size
%               a) data_X and data_Y are technical replicates of the same
%               biological samples (measured at every time point of a
%               longtitudinal microbiome study)
%    
%           *Assumes taxon (OTU) identifiers are provided in the first
%           column and full taxonomies are provided in the last column
%
%   OUTPUTS: 1) Variance decomposition of each taxon (DIVERS_dual_variances.txt)
%               a) Average abundances of each taxon
%               b) Total abundances variances of each taxon
%               c) Biological and technical variances of each taxon
%
%            2) Covariance decomposition for all pairs of taxa 
%               a) Total correlation matrix between all pairs of taxa
%               (DIVERS_dual_cormat_total.txt)
%               b) Biological correlation matrix between all pairs of taxa
%               (DIVERS_dual_cormat_B.txt)
%               c) Technical correlation matrix between all pairs of taxa
%               (DIVERS_dual_cormat_N.txt)
%
%           *Covariance decomposition output reflects abundant OTUs (log10
%           mean absolute abundance > -4). This value was informed by the
%           variance decomposition results. 
%
%           *For large data sets, filtering of abundant OTUs may be
%           required before covariance decomposition analysis
%
%  Brian Ji, 2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Directory containing absolute abundance tables
file_dir = ['./absolute_abundances/'];

%Directory to save DIVERS output to
save_dir = ['./matlab/DIVERS_scripts/DIVERS_output/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%data_X (Samples that represent one of the two technical replicates from each time point)
T_X = readtable([file_dir 'data_X.txt'],'Delimiter','\t');
otu_ids = table2array(T_X(:,1));
tax = table2array(T_X(:,end));
data_X = table2array(T_X(:,2:end-1));
[Mx,Nx] = size(data_X);

%data_Y (Samples that represent the second the two technical replicates from each time point)
T_Y = readtable([file_dir 'data_Y.txt'],'Delimiter','\t');
data_Y = table2array(T_Y(:,2:end-1));
[My,Ny] = size(data_Y);




%% Calculate marginal means and variances of each taxon 

%Perform 500 different re-sampling iterations
marg_means = [];
marg_vars = [];

for i = 1:5e2
    
    data_X_perm = [];
    
    %Randomly draw a sample from data_X, data_Y or data_Z
    for j = 1:Nx
        flip = rand;
        if flip > 1/2
            data_X_perm(:,j) = data_X(:,j);
        else
            data_X_perm(:,j) = data_Y(:,j);
        end
    end
    
    %Estimate mean and variance from this iteration
    marg_means(:,i) = mean(data_X_perm,2);
    marg_vars(:,i) = var(data_X_perm')';
end

%Average over all re-sampling iterations
means = mean(marg_means,2);
vars_total = mean(marg_vars,2);


%% Variance decomposition of taxa abundances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Performing variance decomposition...'])


covs_XY = [];
vars_XmY = [];

%Perfor 500 re-sampling iterations
for i = 1:5e2

    data_X_perm = [];
    data_Y_perm = [];
    
    %Randomly permute X and Y, Z stays the same
    for j = 1:Nx
        flip = rand;
        if flip > .5
            data_X_perm(:,j) = data_X(:,j);
            data_Y_perm(:,j) = data_Y(:,j);
        else
            data_X_perm(:,j) = data_Y(:,j);
            data_Y_perm(:,j) = data_X(:,j);
        end
    end
    
    %Temporal + spatial variance
    cov_XY = [];
    for j = 1:Mx
        covmat_xy = cov(data_X_perm(j,:),data_Y_perm(j,:));
        cov_xy = covmat_xy(1,2);
        cov_XY(j) = cov_xy;
    end
    covs_XY(:,i) = cov_XY';
     
    %Technical
    var_XmY = .5 * var((data_X_perm - data_Y_perm)')';
    vars_XmY(:,i) = var_XmY;
    
end

disp(['Variance decomposition complete!'])


%% Average over re-sampling iterations

%Biological variance
vars_B = mean(covs_XY,2);
vars_B = max(vars_B,0); %Ignore negative estimates of variability

%Technical variance
vars_N = mean(vars_XmY,2);

%Fraction of each component to total variance
vf_N = vars_N ./ (vars_N + vars_B);
vf_B = vars_B ./ (vars_N + vars_B);


%% Covariance decomposition of taxa abundances %%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['Performing covariance decomposition...'])


%% Calculate total covariances
covs_total = [];

%Perform 500 re-sampling iterations
for i = 1:5e2

    data_X_perm = [];

    %Randomly choose a sample from data_X or data_Y
    for j = 1:Nx
        flip = rand;
        if flip > 1/2
            data_X_perm(:,j) = data_X(:,j);
        else
            data_X_perm(:,j) = data_Y(:,j);
        end
    end
    
    %Mean center the data 
    data_X_mc = data_X_perm - repmat(mean(data_X_perm,2),1,Nx);
    
    %Calculate covariance matrix
    covmat_X_perm = 1/(Nx-1) * data_X_mc * data_X_mc';
    
    %Store covariance matrix for this iteration
    covs_total(:,:,i) = covmat_X_perm;  
    
end
covmat_total = mean(covs_total,3);
clear covs_total;


%% Biological covariance decomposition
crosscovs_B = [];

%Perform 500 re-sampling iterations
for i = 1:5e2

    data_X_perm = [];
    data_Y_perm = [];
    
    %One realization of the data
    for j = 1:Nx
        flip = rand;
        if flip > .5
            data_X_perm(:,j) = data_X(:,j);
            data_Y_perm(:,j) = data_Y(:,j);
        else
            data_X_perm(:,j) = data_Y(:,j);
            data_Y_perm(:,j) = data_X(:,j);
        end
    end
    
    %Mean center matrices;
    data_X_mc = data_X_perm - repmat(mean(data_X_perm,2),1,Nx);
    data_Y_mc = data_Y_perm - repmat(mean(data_Y_perm,2),1,Ny);
    
    %Calculate all pairwise covariances for this realization (2
    %permutations)
    covmat_XY = 1/(Nx-1) * data_X_mc*data_Y_mc'; %(X_i,Y_j)
    covmat_YX = 1/(Nx-1) * data_Y_mc*data_X_mc'; %(Y_i,X_j)
    covmat = 1/2 * (covmat_XY + covmat_YX);
    
    %Store covmat for each realization
    crosscovs_B(:,:,i) = covmat;    

end

%Average over all iterations
covmat_B = mean(crosscovs_B,3);
clear crosscovs_B;



 %% Technical covariance decomposition
crosscovs_N = [];

%Peform 500 re-sampling iterations
for i = 1:5e2

    data_X_perm = [];
    data_Y_perm = [];
    
    %One realization of the data
    for j = 1:Nx
    flip = rand;
        if flip > .5
            data_X_perm(:,j) = data_X(:,j);
            data_Y_perm(:,j) = data_Y(:,j);
        else
            data_X_perm(:,j) = data_Y(:,j);
            data_Y_perm(:,j) = data_X(:,j);
        end
    end
   
    data_XmY = data_X_perm - data_Y_perm;
    
    %Mean center matrices;
    data_XmY_mc = data_XmY - repmat(mean(data_XmY,2),1,Nx);
    
    %Calculate all pairwise covariances for this realization (1
    %permutation)
    covmat = 1/(Nx-1) * data_XmY_mc*data_XmY_mc';
    
    %Store covmat for each realization
    crosscovs_N(:,:,i) = .5 * covmat;   

end

%Average over all iterations
covmat_N = mean(crosscovs_N,3);  
clear crosscovs_N;


%% Re-scale covariances to obtain correlations
L = length(vars_total);

%Calculate product of marginal standard deviations for each pair of OTUs
sigxsigy = [];
for i = 1:L
   for j = 1:i
       if vars_total(i) > 0 && vars_total(j) > 0
            sigxsigy(i,j) = sqrt(vars_total(i))*sqrt(vars_total(j));
            sigxsigy(j,i) = sqrt(vars_total(i))*sqrt(vars_total(j));
       else
            sigxsigy(i,j) = 0;
            sigxsigy(j,i) = 0;
       end
   end
end

cormat_total = covmat_total ./ sigxsigy;
cormat_B = covmat_B ./ sigxsigy;
cormat_N = covmat_N ./ sigxsigy;



    
clear covs_XY vars_XmY sigxsigy data_X_perm data_Y_perm  data_X_mc data_Y_mc data_XmY_mc data_XmY 
clear covmat_X_perm covmat_XY covmat_YX


%% Saving

save([save_dir 'matData/DIVERS_dual.mat']);


%% Write to output

%Variances
vars_table = table(otu_ids, means, [vars_B + vars_N], vars_B, vars_N, tax);
vars_table.Properties.VariableNames = {'OTU_ID','Average_abundance','Total_variances','Biological_variances','Technical_variances','Taxonomy'};
writetable(vars_table,[save_dir 'DIVERS_dual_variances.txt'],'Delimiter','\t');

%Covariances (for higher abundance taxon)
high_inds = find(log10(means) > -4); %This cutoff should be informed by technical variance decomposition results
high_otus = table(otu_ids(high_inds),tax(high_inds));
high_otus.Properties.VariableNames = {'OTU_ID','Taxonomy'};
writetable(high_otus,[save_dir 'DIVERS_dual_cormat_otus.txt'],'Delimiter','\t');


dlmwrite([save_dir 'DIVERS_dual_cormat_total.txt'],cormat_total(high_inds,high_inds),'Delimiter','\t');
dlmwrite([save_dir 'DIVERS_dual_cormat_B.txt'],cormat_B(high_inds,high_inds),'Delimiter','\t');
dlmwrite([save_dir 'DIVERS_dual_cormat_N.txt'],cormat_N(high_inds,high_inds),'Delimiter','\t');


 
 

