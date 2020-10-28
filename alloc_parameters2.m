function [p, r] = alloc_parameters2(x, p, r, xi)


     
p_HIV_rltve_cs  = 2;
p_HIV_rltve_pdx = 2;
p_ETB_rltve_pdx = 2;

% --- Non-calibration parameters
if length(x) > xi.calibrating
    p_HIV_rltve_cs  = x(xi.p_HIV_rltve_cs);
    p_HIV_rltve_pdx = x(xi.p_HIV_rltve_pdx);
    p_ETB_rltve_pdx = x(xi.p_ETB_rltve_pdx);
    r.careseeking2  = x(xi.r_cs2)*repmat([1, p_HIV_rltve_cs*ones(1,3)],1,2);
    
    % Yields
    p.prop_sputum = x(xi.prop_sputum);
    p.prop_urine  = x(xi.prop_urine);
     
    % Test sensitivities
    p.Xpert_sn = x(xi.Xpert_sn);
    p.Alere_sn = x(xi.Alere_sn);
    p.Fuji_sn  = x(xi.Fuji_sn);
     
    % other diagnostic values 
    p.prop_symp          = x(xi.prop_symp);
    p.prop_offered_Xpert = x(xi.prop_offered_Xpert);   
    p.Tx_init            = x(xi.Tx_init);
    p.clinical           = x(xi.clinical);
    r.mort_inp           = x(xi.mort_inp); 
    
    % natural history
    p.rel                   = x(xi.rel);
    p.Fast                  = x(xi.Fast);
    r.reactivation          = x(xi.reactivation);
    r.relapse_selfcure      = x(xi.relapse_selfcure);
    r.relapse_treatcure     = x(xi.relapse_treatcure);
    r.stabil                = x(xi.stabil);
    p.cure                  = x(xi.cure);
    p.imm                   = x(xi.imm);
    p.ptb                   = x(xi.ptb);
    r.mort_HIV_excess_low   = x(xi.mort_HIV_excess_low);
    r.mort_HIV_excess_mid   = x(xi.mort_HIV_excess_mid);
    r.mort_HIV_excess_high  = x(xi.mort_HIV_excess_high);
    p.vs                    = x(xi.vs);
    r.inp_leave             = x(xi.inp_leave); 
    r.reactivation_inp      = x(xi.reactivation_inp);
    p.Fast_inp              = x(xi.Fast_inp);
    p.ptb_inp               = x(xi.ptb_inp);
    p.cure_inp              = x(xi.cure_inp); 
    r.self_cure             = x(xi.self_cure);
        
end

% HIV vectors have order: 'h0', 'hlow', 'hmid', 'hhigh', 'ART_p', 'ART_n_low', 'ART_n_mid', 'ART_n_high'
r.beta         = x(xi.r_beta);
p.HIV_suscep   = x(xi.p_HIV_suscep);
r.careseeking  = x(xi.r_careseeking)*repmat([1, p_HIV_rltve_cs*ones(1,3)],1,2);
r.hiv          = x(xi.r_HIV);
r.ART_outp     = x(xi.r_ART_outp);                                          
r.hosp         = x(xi.r_hosp);                                              
r.CD4prog      = x(xi.r_CD4prog);                                           

mat = ones(2,8); 
mat(1,1)     = x(xi.p_Dx(1));
mat(1,2:end) = x(xi.p_Dx(1))/p_HIV_rltve_pdx;
mat(2,:)     = x(xi.p_Dx(1))/p_ETB_rltve_pdx;

p.Dx_TB_com  = mat;

r.mort_TB = x(xi.r_mort_TB); 
r.mort_HIV = x(xi.r_mort_HIV);


