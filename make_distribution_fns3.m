% --- 1 denotes lognormal distribution, 2 denotes beta distribution

data.inc          = [[373 520 691],  [1]];

data.pHIV_coinf   = [[0.55 0.59 0.65],  [2]];
data.HIV_prev     = [[6.9 7.5 8.0]*1e6/56e6,  [2]];
data.nART_inits   = [1.3*[0.8 1 1.2],  [1]];
data.pART_vs      = [[51 55 68]/100,  [2]];

vec = [0.8 1 1.2];
data.pCD4_ARTn    = [[0.76 0.12 0.12]'*vec,  [2 2 2]'];
data.pCD4_ART_out = [[0.41 0.25 0.30]'*vec,  [2 2 2]'];
data.pCD4_ART_inp = [[0.34 0.22 0.35]'*vec,  [2 2 2]']; 

vec = [0.7 1 1.3];
data.hosp         = [5*vec,  [1]];

data.mortality    = [[35 37 39; 51 73 99],  [1 1]'];

vec = [0.9 1 1.1];
data.notifs       = [227999/56e6*1e5*vec,  [1]];


% --- Concatenate the data, taking care that all future comparisons are made in the same order as this concatenation
alldat = [data.inc; data.pHIV_coinf; data.HIV_prev; data.nART_inits; data.pART_vs; data.pCD4_ARTn; data.pCD4_ART_out; data.pCD4_ART_inp; data.hosp; data.mortality; data.notifs];

% --- Get the distribution parameters
dp = zeros(size(alldat,1),2);
inds1 = find(alldat(:,4)==1);
for ii = 1:length(inds1)
    [~, dp(inds1(ii),:)] = get_distribution_fns(alldat(inds1(ii),1:3), 'lognorm');
end
inds2 = find(alldat(:,4)==2);
for ii = 1:length(inds2)
    [~, dp(inds2(ii),:)] = get_distribution_fns(alldat(inds2(ii),1:3), 'beta');
end
dp = dp';

% --- Weights on those data considered most important
weights1 = ones(1,18); weights1([1,4,6:14,16,17]) = 10; 
weights1 = weights1/max(weights1);

% --- Compose the likelihood function - see notes 'Prob distributions for likelihoods'
lhd.fn  = @(y) sum((-log(y(inds1).*dp(2,inds1)) - (log(y(inds1)) - dp(1,inds1)).^2./(2*dp(2,inds1)).^2).*weights1(inds1)) + ...
               sum((dp(1,inds2)-1).*log(y(inds2)) + (dp(2,inds2)-1).*log(1-y(inds2)) - betaln(dp(1,inds2), dp(2,inds2)).*weights1(inds2));
lhd.sgn = -Inf;
     
% --- And compose the least-squares function - see notes 'Prob distributions for likelihoods'
weights2 = ones(1,18); weights2([1,4,6:14,16,17]) = 10; weights2 = weights2/max(weights2);
lsq.fn  = @(y) sum((1-y./alldat(:,2)').^2.*weights2)*100;
lsq.sgn = +Inf;

% --- Finally, set up the boundaries and vector addresses
bounds                        = zeros(2,122);
bounds(2,:)                   = Inf;
bounds(2,xi.p_Dx)             = 1;
bounds(2,xi.r_careseeking)    = 4;

bounds(:,xi.r_CD4prog)             = [0 0.9;0 0.9]'; 
bounds(:,xi.r_cs2)                 = [2 50]';
bounds(:,xi.p_HIV_rltve_cs)        = [1 5]';
bounds(:,xi.p_HIV_rltve_pdx)       = [1 5]';
bounds(:,xi.p_ETB_rltve_pdx)       = [1 5]';
bounds(:,xi.prop_sputum)           = [0.8 1; 0.4 0.8]';
bounds(:,xi.Alere_sn)              = [0 0.05; 0.046 0.237; 0.152 0.389; 0.422 0.696; 0.317 0.518; 0.046 0.237; 0.152 0.389; 0.422 0.696]';
bounds(:,xi.Fuji_sn)               = [0.20 0.50; 0.297 0.585; 0.444 0.725; 0.714 0.914; 0.530 0.831; 0.297 0.585; 0.444 0.725; 0.714 0.914]';
bounds(:,xi.prop_symp)             = [0.66 0.98; 0.76 1.00]';
bounds(:,xi.prop_offered_Xpert)    = ([0.8 1]'*p.prop_offered_Xpert);
bounds(:,xi.prop_urine)            = [0.95 1]';
bounds(:,xi.Xpert_sn)              = [0.76 0.92; 0.70 0.86; 0.70 0.86; 0.70 0.86; 0.76 0.92; 0.70 0.86; 0.70 0.86; 0.70 0.86]';
bounds(:,xi.Tx_init)               = [0.70 0.92; 0.90 1.00; 0.77 0.90]'; 
bounds(:,xi.clinical)              = ([0.8 1.2]'*p.clinical);
bounds(:,xi.mort_inp)              = [0.5 2]';
bounds(:,xi.rel)                   = [0.60 0.80]';
bounds(:,xi.Fast)                  = [0.115 0.16; 0.33 0.51; 0.67 0.82; 0.90 1.00; 0.115 0.16; 0.33 0.51; 0.67 0.82; 0.90 1.00]';
bounds(:,xi.reactivation)          = [0.0003 0.0024; 0.003 0.006; 0.085 0.15; 0.1 0.3; 0.0003 0.0024; 0.003 0.006; 0.085 0.15; 0.1 0.3]';
bounds(:,xi.relapse_selfcure)      = [0.10 0.20]';
bounds(:,xi.relapse_treatcure)     = [0.001 0.003]';
bounds(:,xi.stabil)                = [0.4 0.6]';
bounds(:,xi.cure)                  = [0.77 0.87; 0.75 0.85; 0.75 0.85; 0.75 0.85; 0.75 0.85; 0.75 0.85; 0.75 0.85; 0.75 0.85]';
bounds(:,xi.imm)                   = [0.25 0.85]';
bounds(:,xi.ptb)                   = [0.80 0.90; 0.50 0.70; 0.50 0.70; 0.50 0.70; 0.50 0.70; 0.50 0.70; 0.50 0.70; 0.50 0.70]';
bounds(:,xi.mort_HIV_excess_low)   = [0.02 0.06]';
bounds(:,xi.mort_HIV_excess_mid)   = [0.13 0.21]';
bounds(:,xi.mort_HIV_excess_high)  = [0.42 0.83]';
bounds(:,xi.vs)                    = [0.80 0.95]';
bounds(:,xi.inp_leave)             = [36.5 73; 33.2 60.8; 28.1 52.1]';
bounds(:,xi.reactivation_inp)      = [0.0003 0.0024; 0.003 0.006; 0.085 0.15; 0.1 0.3]';
bounds(:,xi.Fast_inp)              = [0.115 0.16; 0.33 0.51; 0.67 0.82; 0.90 1.00]';
bounds(:,xi.ptb_inp)               = [0.50 0.70; 0.50 0.70; 0.50 0.70; 0.50 0.70]';
bounds(:,xi.cure_inp)              = [0.75 0.85]';
bounds(:,xi.self_cure)             = [0.15 0.25; 0 0.05; 0 0.05; 0 0.05;0.15 0.25; 0 0.05; 0 0.05; 0 0.05]';
prm.bounds = bounds;




