%%  This script performs variance and covariance decomposition in absolute 
%%  abundances using original OTU table as input %%%%%%%%%%%%%%%%%%%%%%%%%%




%% Do some preprocessing of OTU table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read in the metadata

M = readtable('/Users/brianji/Documents/dv_lab/Manuscripts/noise/all_original_data/meta_data/brian_metadata.xlsx');
%rows 1-88 are samples we are interested in (22x3 DIVERS samples, 12 additional
%spatial replicates, 10 additional technical replicates = 88 samples)

M_weights = table2array(M(1:88,10)); % sample weights in mg
M_names = table2array(M(1:88,16)); %sample names
M_days = table2array(M(1:88,5)); %day of time series study


%% Set samples in preferred chronological/preferred order

%First 20 days of DIVERS samples, 12 additional spatial replicates, 10 additional
%technical replicates, two additional days of DIVERS sampes
target_names = cat(1,M_names(1:60),M_names(67:88));
target_names = cat(1,target_names,M_names(61:66));
target_weightz = cat(1,M_weights(1:60),M_weights(67:88));
target_weightz = cat(1,target_weightz,M_weights(61:66));
target_days = cat(1,M_days(1:60),M_days(67:88));
target_days = cat(1,target_days,M_days(61:66));

%Convert weights to numeric
target_weights = [];
for i = 1:length(target_weightz)
    target_weights(i) = str2double(target_weightz(i));
end

%These will be final indices - OTU data will be arranged in this order
time_inds = 1:60;
space_inds = [40 41 61:72];
tech_inds = [41 42 73:82];


%% Read in OTU table data

T = readtable('/Users/brianji/Documents/dv_lab/Manuscripts/noise/all_original_data/raw_data/brian/brian_uclust_55K.txt','Delimiter','\t');
samples = T.Properties.VariableNames; samples = samples(2:end-1); %get sample names
otu_ids = table2array(T(:,1)); %get otu ids
tax = table2array(T(:,end)); %get taxonomic assignments
data = table2array(T(:,2:end-1)); %get OTU count data
[M,N] = size(data);


%% Get OTU data in chronological/preferred order

%Find mapping between target names and sample names
target_inds = [];
for i = 1:length(target_names)
   target_inds = [target_inds;find(strcmp(target_names(i),samples))]; 
end

%Rearrange OTU table in chronological/preferred order
data = data(:,target_inds);

%Rename stuff
samples = target_names;
weights = target_weights;
days = target_days;

%% Split time series data into replicates

Z_inds = time_inds(1:3:60);
X_inds = time_inds(2:3:60);
Y_inds = time_inds(3:3:60);


%% Calculate total bacterial densities in each sample

%Relative abundance of spike-in in each sample
spike_otu_abunds = data(1,:) ./ sum(data,1);

%Calculate scaled absolute density per sample based on spike-in abundance and
%sample weights
abs_abunds = (1-spike_otu_abunds) ./ (spike_otu_abunds .* target_weights);

%Renormalize abundances to mean 1
norm_abs_abunds = abs_abunds ./ mean(abs_abunds);
abund_mat = repmat(norm_abs_abunds,M,1);


%% Convert OTU relative abundances into absolute abundances

%Renormalize relative abundances of OTUs ignoring spike-in OTU
data_rel = data(2:end,:) ./ repmat(sum(data(2:end,:),1),M-1,1); 

%Multiply relative abundances by absolute bacterial densities in each sample 
data_abs = data_rel .* abund_mat(2:end,:);















%% Perform variance decomposition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Split time series samples
data_X = data_abs(:,X_inds);
data_Y = data_abs(:,Y_inds);
data_Z = data_abs(:,Z_inds);
[Mx,Nx] = size(data_X);
[My,Ny] = size(data_Y);
[Mz,Nz] = size(data_Z);


%% Calculate marginal means and variances of each OTU 

%Perform 1,000 different re-sampling iterations
marg_means = [];
marg_vars = [];

for i = 1:1e3
    
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

covs_XZ = [];
covs_XmZY = [];
vars_XmY = [];

%Loop over 1,000 re-sampling iterations
for i = 1:1e3

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
    
    i
end


%% Average over re-sampling iterations

%Temporal variance
vars_T = mean(covs_XZ,2);
vars_T = max(vars_T,0); %Ignore negative estimates of variability

%Spatial sampling variance
vars_S = mean(covs_XmZY,2);
vars_S = max(vars_S,0); %Ignore negative estimates of variability

%Technical variance
vars_N = mean(vars_XmY,2);

% Fraction of each component to total variance
vf_N = vars_N ./ (vars_N + vars_S + vars_T);
vf_S = vars_S ./ (vars_N + vars_S + vars_T);
vf_T = vars_T ./ (vars_N + vars_S + vars_T);



%% Variance decomposition of total bacterial densities %%%%%%%%%%%%%%%%%%%%

% Total absolute abundance densities
abunds = abund_mat(1,:);

%Include days 27 and 48 in this calculation
abunds_X = abunds([X_inds 84 87]);
abunds_Y = abunds([Y_inds 85 88]);
abunds_Z = abunds([Z_inds 83 86]);
abunds_all = [abunds_X; abunds_Y; abunds_Z];

%% First calculate marginal variance of total bacterial densities

%Perform 1,000 re-sampling iterations
for i = 1:1e3
    
    abunds_all_perm = [];
    
    %Randomly shuffle X, Y and Z
    for j = 1:length(abunds_X)
        ord = randperm(3);
        abunds_all_perm(:,j) = abunds_all(ord,j);
    end
    
    var_perm = mean(var(abunds_all_perm')');
    vars_perm(i) = var_perm;
    
end

%Average over 1,000 iterations
abunds_vars_total = mean(vars_perm);


%% Now perform variance decomposition

abunds_covs_XZ = [];
abunds_covs_XmZY = [];
abunds_vars_XmY = [];

%Loop over 1,000 re-sampling iterations
for i = 1:1e3
    
    abunds_X_perm = [];
    abunds_Y_perm = [];
    
    for j = 1:length(abunds_X)
        flip = rand;
        if flip > .5
            abunds_X_perm(j) = abunds_X(j);
            abunds_Y_perm(j) = abunds_Y(j);
        else
            abunds_X_perm(j) = abunds_Y(j);
            abunds_Y_perm(j) = abunds_X(j);
        end
    end 
    
    %Temporal variance
    abunds_covmat_XZ = cov(abunds_X_perm,abunds_Z);
    abunds_cov_XZ = abunds_covmat_XZ(1,2);
    abunds_covs_XZ(i) = abunds_cov_XZ;
    
    %Spatial sampling variance
    abunds_covmat_XmZY = cov(abunds_X_perm-abunds_Z,abunds_Y_perm);
    abunds_cov_XmZY = abunds_covmat_XmZY(1,2);
    abunds_covs_XmZY(i) = abunds_cov_XmZY;
    
    %Technical variance
    abunds_var_XmY = .5 * var([abunds_X_perm - abunds_Y_perm])';
    abunds_vars_XmY(i) = abunds_var_XmY;
    
    i
end

%% Average over re-sampling iterations

%Temporal variance
abunds_vars_T = mean(abunds_covs_XZ,2);
abunds_vars_T = max(abunds_vars_T,0);  %Ignore negative estimates of variability

%Spatial sampling variance
abunds_vars_S = mean(abunds_covs_XmZY,2);
abunds_vars_S = max(abunds_vars_S,0);  %Ignore negative estimates of variability

%Technical variance
abunds_vars_N = mean(abunds_vars_XmY,2);

% Fraction of each component to total variance
abunds_vf_T = abunds_vars_T / (abunds_vars_T + abunds_vars_S + abunds_vars_N);
abunds_vf_S = abunds_vars_S / (abunds_vars_T + abunds_vars_S + abunds_vars_N);
abunds_vf_N = abunds_vars_N / (abunds_vars_T + abunds_vars_S + abunds_vars_N);

%Error in the estimate
abunds_vf_T_std = std(abunds_covs_XZ / (abunds_vars_T + abunds_vars_S + abunds_vars_N));
abunds_vf_S_std = std(abunds_covs_XmZY / (abunds_vars_T + abunds_vars_S + abunds_vars_N));
abunds_vf_N_std = std(abunds_vars_XmY / (abunds_vars_T + abunds_vars_S + abunds_vars_N));










%% Covariance decomposition of OTU abundances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Calculate total covariances
covs_total = [];

%Loop over 1,000 re-sampling iterations
for i = 1:1e3

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
    
    i
end
covmat_total = mean(covs_total,3);


%% Temporal covariance decomposition
crosscovs_T = [];

%Loop over 1,000 re-sampling iterations
for i = 1:1e3

    data_X_perm = [];
    data_Y_perm = [];
    data_Z_perm = data_Z;
    
    %One realization of the data
    for j = 1:Nx
        flip = rand;
        if flip < .5
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
    i
end

%Average over all iterations
covmat_T = mean(crosscovs_T,3);



%% Spatial sampling covariance decomposition
crosscovs_S = [];

%Loop over 1,000 re-sampling iterations
for i = 1:1e3

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
    
    i
end

%Average over all iterations
covmat_S = mean(crosscovs_S,3);



 %% Technical covariance decomposition
crosscovs_N = [];

%Loop over 1,000 re-sampling iterations
for i = 1:1e3

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
    i
end

%Average over all iterations
covmat_N = mean(crosscovs_N,3);   


%% Scale covariances to correlations
L = length(vars_total);

%Calculate product of marginal standard deviations for each pair of OTUs
sigxsigy = [];
for i = 1:L
   for j = 1:i
       if vars_total(i) && vars_total(j) > 0
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
