clear all;

%States:
   %h0    = HIV negative
   %hlow  = HIV CD4 >350 (no ART)
   %hmid  = HIV CD4 200-350 (no ART)
   %hhigh = HIV CD4 <200 (no ART)
   %com   = community 
   %inp   = inpatient 
   %outp  = outpatient 
   %ptb   = pulmonary TB 
   %etb   = extra-pulmonary TB 

states0 = {'U_com','L_com'};                                                                 hivs_com = {'h0','hlow','hmid','hhigh','ART_p','ART_n_low','ART_n_mid','ART_n_high'};
states1 = {'I_com','E_com','SelfCure_com','TreatCure_com','Tx_com','Dx_com'};                hivs_inp = {'ART_p_inp','ART_n_low_inp','ART_n_mid_inp','ART_n_high_inp'};
states2 = {'U_inp','L_inp'};                                                                 tbtypes = {'ptb','etb'};
states3 = {'I_inp','E_inp','SelfCure_inp','TreatCure_inp','Tx_inp','AdmissionDx_inp'};
states4 = {'Dx_outp'};




[i, s, d, lim] = get_addresses({states0, hivs_com}, [], [], [], 0);
[i, s, d, lim] = get_addresses({states1, hivs_com, tbtypes}, i, s, d, lim);
mk = lim;
[i, s, d, lim] = get_addresses({states2, hivs_inp}, i, s, d, lim);
[i, s, d, lim] = get_addresses({states3, hivs_inp, tbtypes}, i, s, d, lim);
[i, s, d, lim] = get_addresses({states4, hivs_com(5:end), tbtypes}, i, s, d, lim);

d = char(d);

%%% Include the auxiliaries 
names = {'inc_total','mortTB','hospitalisation','HIVcom_high', ...
    'HIVcom_mid','HIVcom_low','HIVinp_high','HIVinp_mid','HIVinp_low','Tx_inits', 'morts_hhigh', 'morts_inp'};
lgths = [          3,       2,                1,             1,...
               1,           1,             1,          1,           1,         1,             1,           1];

for ni = 1:length(names)
    inds = lim + [1:lgths(ni)];
    i.aux.(names{ni}) = inds;
    lim = inds(end);
end
i.nx = lim;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Extra states %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct the list of inpatient states
tmp = [states2, states3]; tmp2 = [];
for ti = 1:length(tmp)
    tmp2 = [tmp2, s.(tmp{ti})];
end
s.inpatients = unique(tmp2);

%%% HIV states
s.HIV = unique([s.hlow, s.hmid, s.hhigh, s.ART_p, s.ART_n_low, s.ART_n_mid, s.ART_n_high, s.ART_p_inp, s.ART_n_low_inp, s.ART_n_mid_inp, s.ART_n_high_inp]);
s.ART = unique([s.ART_p, s.ART_n_low, s.ART_n_mid, s.ART_n_high, s.ART_p_inp, s.ART_n_low_inp, s.ART_n_mid_inp, s.ART_n_high_inp]);
s.ART_vs = [s.ART_p, s.ART_p_inp];
s.HIV_nvs = unique([s.hlow, s.hmid, s.hhigh, s.ART_n_low, s.ART_n_mid, s.ART_n_high, s.ART_n_low_inp, s.ART_n_mid_inp, s.ART_n_high_inp]);

%%% TB disease states
s.disease_overall  = unique([s.I_inp, s.E_inp, s.AdmissionDx_inp, s.I_com, s.Dx_com, s.E_com, s.Dx_outp]);
s.disease_overall_HIVn = intersect(s.disease_overall, s.h0);
s.disease_inp      = unique([s.I_inp, s.E_inp, s.AdmissionDx_inp]);
s.disease_com      = unique([s.I_com, s.Dx_com, s.E_com, s.Dx_outp]);
s.infectious       = intersect(s.disease_overall, s.ptb);
s.disease_inp_low  = intersect(s.disease_inp, s.ART_n_low_inp);
s.disease_inp_mid  = intersect(s.disease_inp, s.ART_n_mid_inp);
s.disease_inp_high = intersect(s.disease_inp, s.ART_n_high_inp);
s.disease_inp_ART  = intersect(s.disease_inp, s.ART_p_inp);

%%% Prevalence
s.prevalent_overall             = [s.I_com, s.I_inp, s.E_com, s.E_inp, s.Tx_inp, s.Tx_com, s.Dx_com, s.AdmissionDx_inp, s.Dx_outp];
s.prevalent_inpatient           = [s.I_inp, s.E_inp, s.Tx_inp, s.AdmissionDx_inp];
s.prevalent_outpatient          = [s.I_com, s.E_com, s.Tx_com, s.Dx_com, s.Dx_outp];


%%% Other
s.I_total                       = [s.I_com, s.I_inp];
s.Tx                            = [s.Tx_com, s.Tx_inp];
s.Dx_com_h0                     = intersect(s.Dx_com, s.h0);
s.Dx_com_hp                     = intersect(s.Dx_com, s.HIV);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Selectors and aggregators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Incidence (total) %%%%%
tmp = zeros(3,i.nstates);
tmp(1,s.I_total) = 1;
tmp(2,intersect(s.I_total,s.h0)) = 1;
tmp(3,intersect(s.I_total,s.HIV)) = 1;
agg.inc_total = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.I_total,:)  = 1;
tmp(:, s.I_total) = 0;
sel.inc_total = tmp - diag(diag(tmp));


%%%%% Number of hospitalisations %%%%%
tmp = zeros(i.nstates);
tmp(s.inpatients,setdiff(1:i.nstates,s.inpatients)) = 1;
sel.hospitalisation = tmp - diag(diag(tmp));


%%%%% HIV ART initiation: outpatient %%%%%
inds = setdiff(s.ART, s.inpatients);

tmp = zeros(1,i.nstates);
tmp(1,inds) = 1;
agg.HIVcom = sparse(tmp);

tmp = zeros(1,i.nstates);
tmp(inds, s.hhigh) = 1;                                         % Outpatient ART initiation from hhigh
tmp(mk+1:end,1:mk) = 0; tmp(1:mk,mk+1:end) = 0;
sel.HIVcom_high = tmp - diag(diag(tmp));

tmp = zeros(1,i.nstates);
tmp(inds, s.hmid) = 1;                                          % Outpatient ART initiation from hmid
tmp(mk+1:end,1:mk) = 0; tmp(1:mk,mk+1:end) = 0;
sel.HIVcom_mid = tmp - diag(diag(tmp));

tmp = zeros(1,i.nstates);
tmp(inds, s.hlow) = 1;                                          % Outpatient ART initiation from hlow
tmp(mk+1:end,1:mk) = 0; tmp(1:mk,mk+1:end) = 0;
sel.HIVcom_low = tmp - diag(diag(tmp));


%%%%% HIV ART initiation: inpatients %%%%%
inds = intersect(s.ART, s.inpatients);

tmp = zeros(1,i.nstates);
tmp(1,inds) = 1;
agg.HIVinp = sparse(tmp);


tmp = zeros(i.nstates);
tmp(inds,[s.hlow]) = 1;                               % Inpatient ART initiation from hlow (ignoring those from ART_n_low as notification data won't recount those)
sel.HIVinp_low = tmp - diag(diag(tmp));

tmp = zeros(i.nstates);
tmp(inds,[s.hmid]) = 1;                               % Inpatient ART initiation from hmid
sel.HIVinp_mid = tmp - diag(diag(tmp));

tmp = zeros(i.nstates);
tmp(inds,[s.hhigh]) = 1;                              % Inpatient ART initiation from hhigh
sel.HIVinp_high = tmp - diag(diag(tmp));


%%%%% TB treatment initiations per year %%%%%
tmp = zeros(1,i.nstates);
tmp(1,s.Tx) = 1;
agg.Tx_inits = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.Tx,:) = 1;
tmp(:,s.Tx) = 0;
sel.Tx_inits = sparse(tmp - diag(diag(tmp)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% HIV SPECIFIC %%%%%
r.mort_HIV_excess_low  = 0.04;                                                                                     %Excess mortality in hlow/not virally suppressed 
r.mort_HIV_excess_mid  = 0.17;                                                                                     %Excess mortality in hmid/not virally suppressed 
r.mort_HIV_excess_high = 0.60;                                                                                     %Excess mortality in hhigh/not virally suppressed 
r.mort_HIV_excess_ART  = 0;                                                                                        %Excess mortality in ART virally suppressed 
p.vs                   = 0.88;                                                                                     %Among PLHIV on ART, proportion that are virally suppressed

%%%%% COMMUNITY %%%%% 
%               -->     'h0',   'hlow',  'hmid',  'hhigh',   'ART_p',  'ART_n_low',  'ART_n_mid', 'ART_n_high'
r.reactivation      = [0.001     0.005     0.1      0.2       0.001       0.005          0.1           0.2];      %Rate of reactivation of latent TB infection
p.Fast              = [0.14      0.40      0.80     0.99      0.14        0.40           0.80          0.90];     %Proportion of TB infections undergoing rapid progression
r.self_cure         = [1/6       0         0        0         1/6         0              0             0];        %Rate of spontaneous cure 
p.ptb               = [0.89      0.6       0.6      0.6       0.6         0.6            0.6           0.6];      %Proportion of TB cases with pulmonary TB     
p.cure              = [0.82      0.8       0.8      0.8       0.8         0.8            0.8           0.8];      %Proportion cured after first-line treatment 
r.careseeking2      = 12*[1,     2,        2,       2,        1,          2,             2,            2];        %Rate of repeat careseeking for TB symptoms following missed diagnosis

p.imm               = 0.79;                                                                                       %Reduction in susceptibility due to previous infection
r.relapse_selfcure  = 0.14;                                                                                       %Rate of relapse following self-cure                                                                                            
r.relapse_treatcure = 0.002;                                                                                      %Rate of relapse following treatment cure          
r.stabil            = 0.5;                                                                                        %Stabilisation of relapse risk following treatment           
r.mort              = 1/64;                                                                                       %Rate of background mortality
r.Dx_rou            = 52;                                                                                         %TB treatment initiation delay (takes a week to begin treatment)                           
r.Tx                = 2;                                                                                          %Rate of first-line treatment completion
r.routine           = 3;                                                                                          %Rate of TB testing during ART follow up (every 4 months) 
p.rel               = 0.66;                                                                                       %TB infectiousness in HIV+ (not virally suppressed) vs. HIV-                  

 
%%%%% INPATIENT %%%%% 
%               -->   'ART_p_inp',  'ART_n_low_inp',  'ART_n_mid_inp',  'ART_n_high_inp'
r.reactivation_inp  = [0.001             0.005             0.1              0.2];                                 %Rate of reactivation of latent TB infection, amongst inpatients 
p.Fast_inp          = [0.14              0.40              0.80             0.99];                                %Proportion of TB infections undergoing rapid progression, amongst inpatients
p.ptb_inp           = [0.6               0.6               0.6              0.6];                                 %Proportion of TB cases with pulmonary TB , amongst inpatients

%               -->   'low',  'mid',  'high'
r.inp_leave         = [45.6    42.4    36.5];                                                                     %Average duration of hospitalisation (rate of leaving the hospital)                                                               
p.Tx_init           = [0.86    1.0     0.88];                                                                     %Proportion of diagnosed TB cases successfully initiating treatment

r.self_cure_inp     = 0;                                                                                          %Rate of spontaneous cure, amongst inpatients
p.cure_inp          = 0.8;                                                                                        %Proportion cured after first-line treatment, amongst inpatients 
r.Dx                = 365;                                                                                        %TB treatment initiation delay, amongst inpatients (Takes 1 day to begin treatment)                                                                                    
                                                                                          

%%%%% Test performances %%%%% 
%               --> 'Outpatient', 'inpatient'
p.prop_symp          = [0.82          0.95];                                                                       %Proportion having TB symptoms 
p.prop_offered_Xpert = [0.80          1.00];                                                                       %Proportion of symptomatics receiving an Xpert test

%               --> 'Inpatient', 'outpatient', 'Routine'
p.xpert_extra        = [0             0            0];                                                             %Indicator matrix to switch expanded Xpert coverage on or off

%               -->   'HIV-',  'HIV+'
p.prop_sputum       = [0.90     0.5];                                                                             %Proportion of patients that can produce a sputum sample
p.clinical          = [0.20     0.3];                                                                             %Proportion of negative Xpert results that are clinically diagnosed

p.prop_urine        = 0.99;
r.mort_inp          = 1; 

%               -->    HIV-   hlow   hmid    hhigh    ARTp    ARTlow    ARTmid    ARThigh
p.Xpert_sn          = [0.86   0.79   0.79    0.79     0.86    0.79      0.79      0.79];                          %Xpert sensitivity
p.Xpert_sp          = [0.99   0.98   0.98    0.98     0.99    0.98      0.98      0.98];                          %Xpert specificity 
p.Fuji_sn           = [0.30   0.440  0.606   0.842    0.704   0.440     0.606     0.842];                         %Future LAM test sensitivity
p.Fuji_sp           = [0.974  0.970  0.896   0.850    0.908   0.970     0.896     0.850];                         %Future LAM test specificity 
p.Alere_sn          = [0      0.122  0.264   0.573    0.423   0.122     0.264     0.573];                         %Current LAM test sensitivity 
p.Alere_sp          = [1      0.972  0.928   0.941    0.950   0.972     0.928     0.941];                         %Current LAM test specificity

%               -->   HIV/CD4 status,  Inp/Outp/Routine,  Alere/Fuji
LAM_mat             = zeros(8,                3,              2);                                                 %Indicator matrix to switch on and off different scenarios 

              
             

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Set up the input vector
names1 = {'r_beta','p_HIV_suscep','r_careseeking','r_HIV','r_ART_outp','r_hosp','r_CD4prog','p_Dx','r_mort_TB','r_mort_HIV','r_cs2'};
lgths1 = [       1,             1,              1,      1,           3,       3,          2,     2,          1,           1,      1];

names2 = {'p_HIV_rltve_cs','p_HIV_rltve_pdx','p_ETB_rltve_pdx','prop_sputum','prop_symp','prop_offered_Xpert','prop_urine','Tx_init','clinical','mort_inp','rel'};
lgths2 = [               1,                1,                1,            2,          2,                   2,            1,        3,         2,        1,    1];%,       8,        8,         8];
 
names3 = {'relapse_selfcure','relapse_treatcure','stabil','imm','mort_HIV_excess_low','mort_HIV_excess_mid','mort_HIV_excess_high','Fast','reactivation'};   
lgths3 = [                 1,                  1,       1,    1,                    1,                    1,                     1,     8,             8];   

names4 = {'cure','ptb','vs', 'reactivation_inp','Fast_inp','ptb_inp','cure_inp','Fuji_sn','Xpert_sn','Alere_sn','self_cure', 'inp_leave'};
lgths4 = [     8,    8,   1,                  4,         4,        4,         1,        8,         8,         8,          8,           3];


names = [names1, names2, names3, names4];
lgths = [lgths1, lgths2, lgths3, lgths4];


xi = []; lim = 0;
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    xi.(names{ii}) = inds;
    lim = inds(end);
end
xi.calibrating = xi.r_mort_HIV;                                            % Last index of params being calibrated
xi.nx = lim;

make_distribution_fns3;
Get_HIV_data2;


prm.p = p; prm.r = r; prm.bounds = bounds;
ref.i = i; ref.s = s; ref.d = d; ref.xi = xi;
gps.tbtypes = tbtypes; gps.hivs_com = hivs_com; gps.hivs_inp = hivs_inp;


obj = @(x) get_objective2(x, prm, ref, sel, agg, gps, lsq, HIVdat, LAM_mat);

% Starting point for calibration

x0 = [11.8884    1.8587    1.9877    0.0188    0.0203    0.0969    0.1938    0.0372    0.1272    0.2399    0.2427    0.9000    0.7628    0.3973    0.3702    0.1383   12.9369    1.0227    4.9496    4.4986];

x0(xi.prop_sputum)           = p.prop_sputum; 
x0(xi.prop_symp)             = p.prop_symp;   
x0(xi.prop_offered_Xpert)    = p.prop_offered_Xpert; 
x0(xi.prop_urine)            = p.prop_urine;  
x0(xi.Tx_init)               = p.Tx_init; 
x0(xi.clinical)              = p.clinical; 
x0(xi.mort_inp)              = r.mort_inp; 
x0(xi.rel)                   = p.rel; 
x0(xi.relapse_selfcure)      = r.relapse_selfcure; 
x0(xi.relapse_treatcure)     = r.relapse_treatcure;
x0(xi.stabil)                = r.stabil; 
x0(xi.imm)                   = p.imm; 
x0(xi.mort_HIV_excess_low)   = r.mort_HIV_excess_low; 
x0(xi.mort_HIV_excess_mid)   = r.mort_HIV_excess_mid; 
x0(xi.mort_HIV_excess_high)  = r.mort_HIV_excess_high; 
x0(xi.Fast)                  = p.Fast; 
x0(xi.reactivation)          = r.reactivation; 
x0(xi.cure)                  = p.cure;
x0(xi.ptb)                   = p.ptb; 
x0(xi.vs)                    = p.vs; 
x0(xi.reactivation_inp)      = r.reactivation_inp; 
x0(xi.Fast_inp)              = p.Fast_inp;
x0(xi.ptb_inp)               = p.ptb_inp; 
x0(xi.cure_inp)              = p.cure_inp; 
x0(xi.Fuji_sn)               = p.Fuji_sn;     
x0(xi.Xpert_sn)              = p.Xpert_sn;    
x0(xi.Alere_sn)              = p.Alere_sn;    
x0(xi.self_cure)             = r.self_cure; 
x0(xi.inp_leave)             = r.inp_leave; 


vec = linspace(0.8, 1.2, 11);



opt = 1;
if opt == 0
    % No optimisation - for testing purposes
    x1 = x0;
elseif opt == 1
    % Least-squares optimisation
    obj = @(x) get_objective2(x, prm, ref, sel, agg, gps, lsq, HIVdat, LAM_mat);
    x1 = fminsearchbnd(obj, x0, bounds(1,:), bounds(2,:), optimset('PlotFcns',@optimplotfval));
else
    nsam = 1e6;

    
    obj(x0);
    
    
    obj = @(x) get_objective2(x, prm, ref, sel, agg, gps, lhd, HIVdat, LAM_mat);
    
    load cov0b; 
    
Mat = 0.0001*eye(122);     
Mat(1:20,1:20) = cov0;
cov01 = Mat;

    
    [xsto, outsto, history, accept_rate] = MCMC_adaptive(obj, x0, nsam, 1, [], [], true, cov01);
        
    
    fprintf('\n');
    save MCMC_res_MAIN2c;
end

save('estims.mat');

