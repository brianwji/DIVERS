% DIVERS: Temporal, spatial, and technical variance and covariance
%   contribution estimates from absolute bacterial abundance data
%
%   *User required to specify input and saving directories
%
%  INPUTS: 1) Three absolute abundance tables (data_X, data_Y, data_Z) of
%             equal size
%               a) data_X and data_Y are technical replicates of the same
%               biological samples (measured at every time point of a
%               longtitudinal microbiome study)
%               b) data_Z is a second replicate (from a second spatial
%               location at every time point of a longitudinal microbiome
%               study)
%    
%           *Assumes taxon (OTU) identifiers are provided in the first
%           column and full taxonomies are provided in the last column
%
%   OUTPUTS: 1) Variance decomposition of each taxon (DIVERS_variances.txt)
%               a) Average abundances of each taxon
%               b) Total abundances variances of each taxon
%               c) Temporal, spatial, technical variances of each taxon
%
%            2) Covariance decomposition for all pairs of taxa 
%               a) Total correlation matrix between all pairs of taxa
%               (DIVERS_cormat_total.txt)
%               b) Temporal correlation matrix between all pairs of taxa
%               (DIVERS_cormat_T.txt)
%               c) Spatial correlation matrix between all pairs of taxa
%               (DIVERS_cormat_S.txt)
%               d) Technical correlation matrix between all pairs of taxa
%               (DIVERS_cormat_N.txt)
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
file_dir = ['/Path/To/.../absolute_abundances/'];

%Directory to save DIVERS output to
save_dir = ['/Path/To/.../DIVERS_output/'];

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

%data_Z (Samples that represent the second spatial replicate from each time point)
T_Z = readtable([file_dir 'data_Z.txt'],'Delimiter','\t');
data_Z = table2array(T_Z(:,2:end-1));
[Mz,Nz] = size(data_Z);




%% Calculate marginal means and variances of each OTU 

%Perform 500 different re-sampling iterations
marg_means = [];
marg_vars = [];

for i = 1:5e2
    
    data_X_perm = [];
    
    %Randomly draw a sample from data_X, data_Y or data_Z
    for j = 1:Nx
        flip = rand;
        if flip > (2/3)
            data_X_perm(:,j) = data_X(:,j);
        elseif flip > (1/3)
            data_X_perm(:,j) = data_Y(:,j);
        else
            data_X_perm(:,j) = data_Z(:,j);
        end
    end
    
    %Estimate mean and variance from this iteration
    marg_means(:,i) = mean(data_X_perm,2);
    marg_vars(:,i) = var(data_X_perm')';
end

%Average over all re-sampling iterations
means = mean(marg_means,2);
vars_total = mean(marg_vars,2);


%% Variance decomposition of OTU abundances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Performing variance decomposition...'])


covs_XZ = [];
covs_XmZY = [];
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
    
    %Temporal variance
    cov_XZ = [];
    for j = 1:Mx
        covmat_xz = cov(data_X_perm(j,:),data_Z(j,:));
        cov_xz = covmat_xz(1,2);
        cov_XZ(j) = cov_xz;
    end
    covs_XZ(:,i) = cov_XZ';
    
    %Spatial sampling variance
    cov_XmZY = [];
    for j = 1:Mx
        covmat_xmzy = cov(data_X_perm(j,:)-data_Z(j,:),data_Y_perm(j,:));
        cov_xmzy = covmat_xmzy(1,2);
        cov_XmZY(j) = cov_xmzy;
    end
    covs_XmZY(:,i) = cov_XmZY';
    
    %Technical
    var_XmY = .5 * var((data_X_perm - data_Y_perm)')';
    vars_XmY(:,i) = var_XmY;
    
end
disp(['Variance decomposition complete!'])


%% Average over re-sampling iterations

%Temporal variance
vars_T = mean(covs_XZ,2);
vars_T = max(vars_T,0); %Ignore negative estimates of variability

%Spatial sampling variance
vars_S = mean(covs_XmZY,2);
vars_S = max(vars_S,0); %Ignore negative estimates of variability

%Technical variance
vars_N = mean(vars_XmY,2);

%Fraction of each component to total variance
vf_N = vars_N ./ (vars_N + vars_S + vars_T);
vf_S = vars_S ./ (vars_N + vars_S + vars_T);
vf_T = vars_T ./ (vars_N + vars_S + vars_T);



%% Covariance decomposition of OTU abundances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Performing covariance decomposition...'])


%% Calculate total covariances
covs_total = [];

%Perform 500 re-sampling iterations
for i = 1:5e2

    data_X_perm = [];

    %Randomly choose a sample from of data_X, data_Y, data_Z
    for j = 1:Nx
        flip = rand;
        if flip > (2/3)
            data_X_perm(:,j) = data_X(:,j);
        elseif flip > (1/3)
            data_X_perm(:,j) = data_Y(:,j);
        else
            data_X_perm(:,j) = data_Z(:,j);
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


%% Temporal covariance decomposition
crosscovs_T = [];

%Perform 500 re-sampling iterations
for i = 1:5e2

    data_X_perm = [];
    data_Y_perm = [];
    data_Z_perm = data_Z;
    
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
    data_Z_mc = data_Z_perm - repmat(mean(data_Z_perm,2),1,Nz);
    
    %Calculate all pairwise covariances for this realization (2
    %permutations)
    covmat_XZ = 1/(Nx-1) * data_X_mc*data_Z_mc'; %(X_i,Z_j)
    covmat_ZX = 1/(Nx-1) * data_Z_mc*data_X_mc'; %(Z_i,X_j)
    covmat = 1/2 * (covmat_XZ + covmat_ZX);
    
    %Store covmat for each realization
    crosscovs_T(:,:,i) = covmat;    
    
end

%Average over all iterations
covmat_T = mean(crosscovs_T,3);
clear crosscovs_T;



%% Spatial sampling covariance decomposition
crosscovs_S = [];

%Peform 500 re-sampling iterations
for i = 1:5e2

    data_X_perm = [];
    data_Y_perm = [];
    data_Z_perm = data_Z;
    
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
    
    data_XmZ = data_X_perm - data_Z_perm;
    
    %Mean center matrices;
    data_XmZ_mc = data_XmZ - repmat(mean(data_XmZ,2),1,Nx);
    data_Y_mc = data_Y_perm - repmat(mean(data_Y_perm,2),1,Ny);
    
    %Calculate all pairwise covariances for this realization (2
    %permutations)
    covmat_XZY = 1/(Nx-1) * data_XmZ_mc*data_Y_mc';
    covmat_YXZ = 1/(Nx-1) * data_Y_mc*data_XmZ_mc';
    covmat = 1/2 * (covmat_XZY + covmat_YXZ);
    
    %Store covmat for each realization
    crosscovs_S(:,:,i) = covmat;  
    
end

%Average over all iterations
covmat_S = mean(crosscovs_S,3);
clear crosscovs_S;



 %% Technical covariance decomposition
crosscovs_N = [];

%Peform 500 re-sampling iterations
for i = 1:5e2

    data_X_perm = [];
    data_Y_perm = [];
    data_Z_perm = data_Z;
    
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

disp(['Covariance decomposition complete!'])


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
cormat_T = covmat_T ./ sigxsigy;
cormat_S = covmat_S ./ sigxsigy;
cormat_N = covmat_N ./ sigxsigy;

    
clear covs_XZ covs_XmZY vars_XmY sigxsigy data_X_perm data_Y_perm data_Z_perm data_X_mc data_Y_mc data_Z_mc data_XmY_mc data_XmZ_mc data_XmY data_XmZ



%% Saving

save([save_dir 'DIVERS.mat']);

%% Write to output

%Variances
vars_table = table(otu_ids, means, [vars_T + vars_S + vars_N], vars_T, vars_S, vars_N, tax);
vars_table.Properties.VariableNames = {'OTU_ID','Average_abundance','Total_variances','Temporal_variances','Spatial_variances','Technical_variances','Taxonomy'};
writetable(vars_table,[save_dir 'DIVERS_variances.txt'],'Delimiter','\t');

%Covariances (for higher abundance taxon)
high_inds = find(log10(means) > -4); %This cutoff should be informed by technical variance decomposition results
high_otus = table(otu_ids(high_inds),tax(high_inds));
high_otus.Properties.VariableNames = {'OTU_ID','Taxonomy'};
writetable(high_otus,[save_dir 'DIVERS_cormat_otus.txt'],'Delimiter','\t');

dlmwrite([save_dir 'DIVERS_cormat_total.txt'],cormat_total(high_inds,high_inds),'Delimiter','\t');
dlmwrite([save_dir 'DIVERS_cormat_T.txt'],cormat_T(high_inds,high_inds),'Delimiter','\t');
dlmwrite([save_dir 'DIVERS_cormat_S.txt'],cormat_S(high_inds,high_inds),'Delimiter','\t');
dlmwrite([save_dir 'DIVERS_cormat_N.txt'],cormat_N(high_inds,high_inds),'Delimiter','\t');


 
 

