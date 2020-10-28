function [out, aux] = get_objective2(x, prm, ref, sel, agg, gps, calfn, HIVdat, LAM_mat)

r = prm.r; p = prm.p;
i = ref.i; s = ref.s; xi = ref.xi;


mat = diff([prm.bounds(1,1:length(x)); x; prm.bounds(2,1:length(x))],1);

if min(mat(:)) < 0
    out = calfn.sgn; aux = nan;
else
    [p, r] = alloc_parameters2(x, p, r, xi);
    
    %----------Setting up the necessary models ------------------------------
    
    % p.cure scale up (p.cure = 60% in 1989-1999, linear increase to 80% from 2000-2016)
    p2 = p; r2 = r;
    M2 = make_model_wEPTB(p2, r2, i, s, gps, LAM_mat);
    
    % With HIV, but poorer TB care
    p1 = p; r1 = r;
    p1.cure      = p.cure*0.75;
    p1.Dx_TB_com = p.Dx_TB_com/1.2;
    r1.ART_outp  = 0*r1.ART_outp;
    M1 = make_model_wEPTB(p1, r1, i, s, gps, LAM_mat);
    
    % Equilibrium: Without HIV or ART
    p0 = p1; r0 = r1;
    r0.hiv  = 0;
    r0.hosp = [0 0 0];
    M0 = make_model_wEPTB(p0, r0, i, s, gps, LAM_mat);
    
    % Solving equilibrium model
    init = zeros(1,i.nx); seed = 1e-6;
    init(i.U_com.h0) = (1-seed);
    init(i.I_com.h0.ptb) = seed;
    
    % To equilibrium in absence of HIV
    geq0 = @(t,in)goveqs_basis2(t, in, M0, i, s, r0, p0, sel, agg);
    [t0, soln0] = ode15s(geq0, [0:3e3], init, odeset('NonNegative',[1:i.nstates]));
    
    % Now introduce HIV (but no ART)
    init = soln0(end,:);
    geq1 = @(t,in)goveqs_basis_HIV(t, in, M1, i, s, r1, p1, sel, agg, HIVdat);
    [t1, soln1] = ode15s(geq1, [1989:2001], init, odeset('NonNegative',[1:i.nstates]));
    
    % Introducing ART and improvements in TB programme
    init = soln1(end,:);
    [t2, soln2] = ode15s(@(t,in) goveqs_scaleup(t, in, M1, M2, [2000 2016], i, s, r2, p2, sel, agg, HIVdat), [2000:2016], init, odeset('NonNegative',[1:i.nstates]));
    
    sfin   = soln2(end,:);
    getdif = @(cols) diff(soln2(:,cols),[],1)*1e5;
    
    
    % --- Calibration outputs
    
    % TB Incidence
    inc = getdif(i.aux.inc_total);
    
    % Proportion TB/HIV coinfection
    pHIV_coinf = inc(end,3)/inc(end,1);
    
    % HIV prevalence
    HIV_prev = sum(sfin(s.HIV))/sum(sfin(1:i.nstates));
    
    % Number of ART initiations 
    tmpai = getdif([i.aux.HIVcom_high, i.aux.HIVcom_mid, i.aux.HIVcom_low, i.aux.HIVinp_high, i.aux.HIVinp_mid, i.aux.HIVinp_low]);
    nART_inits  = sum(tmpai(end,:))/1e3;
        
    % Proportion HIV that are virally suppressed
    pARTvs = sum(sfin(s.ART_vs))/sum(sfin(s.HIV));
    
    % CD4 distribution amongst those not on ART
    pCD4_ARTn = [sum(sfin(s.hlow)),sum(sfin(s.hmid)),sum(sfin(s.hhigh))]/sum(sfin(setdiff(s.HIV,s.ART)));
    
    % TB mortality
    rTBmort  = getdif(i.aux.mortTB);
    
    % Proportion of high/mid/low from ART initiations in outpatients
    mat = getdif([i.aux.HIVcom_high, i.aux.HIVcom_mid, i.aux.HIVcom_low]);
    vec = mat(end,:);
    pCD4_outp = vec/sum(vec);
    
    % Proportion of high/mid/low from ART initiations in hospital
    mat = getdif([i.aux.HIVinp_high, i.aux.HIVinp_mid, i.aux.HIVinp_low]);
    vec = mat(end,:);
    pCD4_hosp = vec/sum(vec);
    
    % Numbers being hospitalised
    tmp = diff(soln2(:,i.aux.hospitalisation),1);
    nhosp = tmp(end)/sum(sfin(s.HIV))*100;
    
    % TB treatment initiations
    tmp = getdif(i.aux.Tx_inits);
    Tx_inits = tmp(end);
    
    y = [inc(end,1), pHIV_coinf, HIV_prev, nART_inits, pARTvs, pCD4_ARTn, pCD4_outp, pCD4_hosp, nhosp, rTBmort(end,:), Tx_inits];
    out = calfn.fn(y);
    
    % --- Other outputs
    % TB prevalence
    prev = sum(sfin(s.prevalent_overall))*1e5;
    
    % --- Store them as outputs
    aux.inc          = inc(end,1);
    aux.pHIV_coinf   = pHIV_coinf;
    aux.HIV_prev     = HIV_prev;
    aux.nART_inits   = nART_inits;
    aux.pART_vs      = pARTvs;
    aux.pCD4_ARTn    = pCD4_ARTn;
    aux.pCD4_ART_out = pCD4_outp;
    aux.pCD4_ART_inp = pCD4_hosp;
    aux.hosp         = nhosp;
    aux.mortality    = rTBmort(end,:);
    aux.notifs       = Tx_inits;
    aux.ART_inits    = [sum(tmpai(end,1:3)), sum(tmpai(end,4:6))];
    aux.soln         = [soln2, t2];
    aux.M2           = M1;
    aux.prev         = prev;
    
    aux.tmp          = getdif(i.aux.Tx_inits);
end

