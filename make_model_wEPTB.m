function M = make_model_wEPTB(p, r, i, s, gps, LAM_mat)


% --- Get the linear rates ------------------------------------------------

% --- COMMUNITY MODEL -----------------------------------------------------
m = zeros(i.nstates);
for ih = 1:length(gps.hivs_com)
    hivs_com = gps.hivs_com{ih};
    
    f1    = @(st) i.(st).(hivs_com);
    L_com = f1('L_com');
    
    % --- Reactivation
    source = L_com; destins = [i.I_com.(hivs_com).ptb, i.I_com.(hivs_com).etb]; rates = r.reactivation(ih)*[p.ptb(ih), 1-p.ptb(ih)]';
    m(destins, source) = m(destins, source) + rates;
    
    for it = 1:length(gps.tbtypes)
        tbtypes = gps.tbtypes{it};
        
        f2             = @(st) i.(st).(hivs_com).(tbtypes);
        I_com          = f2('I_com');
        Dx_com         = f2('Dx_com');
        Tx_com         = f2('Tx_com');
        E_com          = f2('E_com');
        SelfCure_com   = f2('SelfCure_com');
        TreatCure_com  = f2('TreatCure_com');
        
        % --- Primary careseeking
        source = I_com; destin = Dx_com; rate = r.careseeking(ih);
        m(destin, source) = m(destin, source) + rate;
        
       
        % --- Diagnosis
        p_Dx_0           = p.Dx_TB_com(it, ih);
        p_Dx_LAM         = p.prop_urine*(p.Fuji_sn(ih)*LAM_mat(ih,3,2) + p.Alere_sn(ih)*LAM_mat(ih,3,1));
        p_clinical_r     = [p.clinical(1), p.clinical(2)*ones(1,7); p.clinical(1), p.clinical(2)*ones(1,7)];
        p_Dx_clinical    = (1-(p_Dx_0+p_Dx_LAM)).*p_clinical_r(it,ih);
        p_prop_sputum_r = [p.prop_sputum(1), p.prop_sputum(2)*ones(1,7)];         
        p_Dx_Xpert_extra = p_prop_sputum_r(ih).*p.Xpert_sn(ih)*(2-it)*p.xpert_extra(3); 
        p_Dx = 1 - ((1-p_Dx_0)*(1-p_Dx_clinical)*(1-p_Dx_LAM)*(1-p_Dx_Xpert_extra)); 
        p_Tx_link = p_Dx*p.Tx_init(3);
        
       
        source = Dx_com; destin = Tx_com;
        rate = r.Dx_rou*p_Tx_link;
        m(destin, source) = m(destin, source) + rate;
        
        source = Dx_com; destin = E_com;
        rate = r.Dx_rou*(1-p_Tx_link);
        m(destin, source) = m(destin, source) + rate;
        

        
        % --- Treatment
        % Successful cure
        source = Tx_com; destin = TreatCure_com;  rate = r.Tx*p.cure(ih);
        m(destin, source) = m(destin, source) + rate;
        
        % Failed treatment
        source = Tx_com; destin = E_com;   rate = r.Tx*(1-p.cure(ih));
        m(destin, source) = m(destin, source) + rate;
        
        
        % --- Secondary careseeking
        source = E_com; destin = Dx_com; rate = r.careseeking2(ih);
        m(destin, source) = m(destin, source) + rate;
        
        
        % --- Relapse
        sources = SelfCure_com;  destin = I_com;  rates = r.relapse_selfcure;
        m(destin, sources) = m(destin, sources) + rates;
        
        sources = TreatCure_com;  destin = I_com;  rates = r.relapse_treatcure;
        m(destin, sources) = m(destin, sources) + rates;
        
        
        % --- Stabilisation of relapse
        source = SelfCure_com;  destin = TreatCure_com; rate = r.stabil;
        m(destin, source) = m(destin, source) + rate;
        
        
        % --- Self cure
        rate = r.self_cure(ih);
        for source = intersect(intersect(s.disease_com,s.(tbtypes)),s.(hivs_com))
            destin = SelfCure_com;
            m(destin, source) = m(destin, source) + rate;
        end
    end
end


% --- ART initiation (community) 

uninfected_com = [s.U_com, s.L_com, s.SelfCure_com, s.TreatCure_com, s.Tx_com];
uninfected_inp = [s.U_inp, s.L_inp, s.SelfCure_inp, s.TreatCure_inp, s.Tx_inp];
community = [s.U_com, s.L_com, s.I_com, s.E_com, s.Tx_com, s.SelfCure_com, s.TreatCure_com];
inpatient = [s.U_inp, s.L_inp, s.I_inp, s.E_inp, s.Tx_inp, s.SelfCure_inp, s.TreatCure_inp];


mART = zeros(i.nstates);
hsets = {'low','mid','high'};
for ih = 1:length(hsets)
    hset = hsets{ih};
    s_h   = s.(['h',hset]);
    ART_n = s.(['ART_n_',hset]);
    
    % >> First time staring ART
    % Amongst TB negatives
    inds = sub2ind([i.nstates, i.nstates], intersect(uninfected_com,s.ART_p), intersect(uninfected_com,s_h));
    mART(inds) = mART(inds) + r.ART_outp(ih)*p.vs;
    inds = sub2ind([i.nstates, i.nstates], intersect(uninfected_com, ART_n), intersect(uninfected_com,s_h));
    mART(inds) = mART(inds) + r.ART_outp(ih)*(1-p.vs);
    
    % Amongst pre-careseeking TB cases (I_com_sym/I_com_asym)
    inds = sub2ind([i.nstates, i.nstates], intersect(s.Dx_outp, s.ART_p), intersect(s.I_com,s_h));
    mART(inds) = mART(inds) + r.ART_outp(ih)*p.vs;
    inds = sub2ind([i.nstates, i.nstates], intersect(s.Dx_outp, ART_n),   intersect(s.I_com,s_h));
    mART(inds) = mART(inds) + r.ART_outp(ih)*(1-p.vs);
    
    % Amongst between-careseeking TB cases (E_com)
    inds = sub2ind([i.nstates, i.nstates], intersect(s.Dx_outp, s.ART_p), intersect(s.E_com,s_h));
    mART(inds) = mART(inds) + r.ART_outp(ih)*p.vs;
    inds = sub2ind([i.nstates, i.nstates], intersect(s.Dx_outp, ART_n),   intersect(s.E_com,s_h));
    mART(inds) = mART(inds) + r.ART_outp(ih)*(1-p.vs);
    
    % Amongst TB positives awaiting diagnosis (Dx_com)
    inds = sub2ind([i.nstates, i.nstates], intersect(s.Dx_outp, s.ART_p), intersect(s.Dx_com,s_h));
    mART(inds) = mART(inds) + r.ART_outp(ih)*p.vs;
    inds = sub2ind([i.nstates, i.nstates], intersect(s.Dx_outp, ART_n),   intersect(s.Dx_com,s_h));
    mART(inds) = mART(inds) + r.ART_outp(ih)*(1-p.vs);
    
    % >> Re-initiating ART
    % Amongst TB negatives
    inds = sub2ind([i.nstates, i.nstates], intersect(uninfected_com,s.ART_p), intersect(uninfected_com,ART_n));
    mART(inds) = mART(inds) + r.ART_outp(ih)*p.vs;
    inds = sub2ind([i.nstates, i.nstates], intersect(uninfected_com, ART_n), intersect(uninfected_com,ART_n));
    mART(inds) = mART(inds) + r.ART_outp(ih)*(1-p.vs);
    
    % Amongst pre-careseeking TB cases (I_com_sym/I_com_asym)
    inds = sub2ind([i.nstates, i.nstates], intersect(s.Dx_outp, s.ART_p), intersect(s.I_com,ART_n));
    mART(inds) = mART(inds) + r.ART_outp(ih)*p.vs;
    inds = sub2ind([i.nstates, i.nstates], intersect(s.Dx_outp, ART_n),   intersect(s.I_com,ART_n));
    mART(inds) = mART(inds) + r.ART_outp(ih)*(1-p.vs);
    
    % Amongst between-careseeking TB cases (E_com)
    inds = sub2ind([i.nstates, i.nstates], intersect(s.Dx_outp, s.ART_p), intersect(s.E_com,ART_n));
    mART(inds) = mART(inds) + r.ART_outp(ih)*p.vs;
    inds = sub2ind([i.nstates, i.nstates], intersect(s.Dx_outp, ART_n),   intersect(s.E_com,ART_n));
    mART(inds) = mART(inds) + r.ART_outp(ih)*(1-p.vs);
    
    % Amongst TB positives awaiting diagnosis (Dx_com)
    inds = sub2ind([i.nstates, i.nstates], intersect(s.Dx_outp, s.ART_p), intersect(s.Dx_com,ART_n));
    mART(inds) = mART(inds) + r.ART_outp(ih)*p.vs;
    inds = sub2ind([i.nstates, i.nstates], intersect(s.Dx_outp, ART_n),   intersect(s.Dx_com,ART_n));
    mART(inds) = mART(inds) + r.ART_outp(ih)*(1-p.vs);
end
M.ART = sparse(mART - diag(sum(mART,1)));


% --- TB screening as part of ART (community) 

hsets = {'p','n_low','n_mid','n_high'};

for ih = 1:length(hsets)
    
    ih_test = ih+4;
    
    for it = 1:length(gps.tbtypes)
        tbtypes = gps.tbtypes{it};

        % --- Testing at ART initiation 
        p_dx_Xpert       = p.prop_sputum(2)*p.prop_offered_Xpert(1)*(p.Xpert_sn(ih_test)*(2-it));
        p_dx_LAM         = p.prop_urine*(p.Fuji_sn(ih_test)*LAM_mat(ih_test,2,2) + p.prop_urine*p.Alere_sn(ih_test)*LAM_mat(ih_test,2,1));
        p_dx_clinical    = (1-(p_dx_Xpert+p_dx_LAM))*p.clinical(2);
        p_dx_Xpert_extra = p.prop_sputum(2)*p.Xpert_sn(ih_test)*(2-it)*p.xpert_extra(2); 
        p_dx = 1 - ((1-p_dx_Xpert)*(1-p_dx_clinical)*(1-p_dx_LAM)*(1-p_dx_Xpert_extra)); 
        p.Tx_link_outp = p.prop_symp(1)*p.Tx_init(1)*p_dx;
                     
         
        hset = hsets{ih};
        ART = s.(['ART_',hset]);
        
        f4 = @(st) intersect(intersect(st, ART),s.(tbtypes));
        
        Tx_com       = f4(s.Tx_com);
        Dx_outp      = f4(s.Dx_outp);
        E_com        = f4(s.E_com);
        I_com        = f4(s.I_com);
        Dx_com       = f4(s.Dx_com);
        
        inds = sub2ind([i.nstates, i.nstates], Tx_com, Dx_outp);
        m(inds) = m(inds) + r.Dx*p.Tx_link_outp;
        
        inds = sub2ind([i.nstates, i.nstates], E_com, Dx_outp);
        m(inds) = m(inds) + r.Dx*(1-p.Tx_link_outp);
       
        
        % --- Testing during ART (zeroed out - assume no routine TB testing)
        inds = sub2ind([i.nstates, i.nstates], Tx_com, I_com);
        m(inds) = m(inds) + r.Dx*r.routine*p.Tx_link_outp*0;
        
        inds = sub2ind([i.nstates, i.nstates], Tx_com, E_com);
        m(inds) = m(inds) + r.Dx*r.routine*p.Tx_link_outp*0;
        
        inds = sub2ind([i.nstates, i.nstates], Tx_com, Dx_com);
        m(inds) = m(inds) + r.Dx*r.routine*p.Tx_link_outp*0;
    end
end



% -------------------------------------------------------------------------

% --- INPATIENT MODEL -----------------------------------------------------

% --- INPATIENT cascade

for ih = 1:length(gps.hivs_inp)
    hivs_inp = gps.hivs_inp{ih};
    
    f3    = @(st) i.(st).(hivs_inp);
    L_inp = f3('L_inp');
    
    % --- Reactivation
    source = L_inp; destins = [i.I_inp.(hivs_inp).ptb, i.I_inp.(hivs_inp).etb]; rates = r.reactivation_inp(ih)*[p.ptb_inp(ih), 1-p.ptb_inp(ih)]';
    m(destins, source) = m(destins, source) + rates;
    
    for it = 1:length(gps.tbtypes)
        tbtypes = gps.tbtypes{it};
        
        f4                   = @(st) i.(st).(hivs_inp).(tbtypes);
        I_inp                = f4('I_inp');
        Tx_inp               = f4('Tx_inp');
        E_inp                = f4('E_inp');
        SelfCure_inp         = f4('SelfCure_inp');
        TreatCure_inp        = f4('TreatCure_inp');
        AdmissionDx_inp      = f4('AdmissionDx_inp');
        
        % --- Treatment outcomes - success and failure
        source  = Tx_inp;
        destins =      [TreatCure_inp, E_inp];
        rates   = r.Tx*[p.cure_inp,    1-p.cure_inp]';
        m(destins, source) = m(destins, source) + rates;
        
        % --- Relapse
        sources = [SelfCure_inp, TreatCure_inp];
        destin  = I_inp;
        rates   = [r.relapse_selfcure, r.relapse_treatcure];
        m(destin, sources) = m(destin, sources) + rates;
        
        % --- Stabilisation of relapse
        source = SelfCure_inp;  destin = TreatCure_inp; rate = r.stabil;
        m(destin, source) = m(destin, source) + rate;
        
        % --- Self cure
        rate = r.self_cure_inp;
        for source = intersect(intersect(s.disease_inp,s.(tbtypes)),s.(hivs_inp))
            destin = SelfCure_inp;
            m(destin, source) = m(destin, source) + rate;
        end
        
        ih_test = ih + 4;
        
        % --- Diagnosis   
        
        p_dx_Xpert       = p.prop_sputum(2)*p.prop_offered_Xpert(2)*(p.Xpert_sn(ih_test)*(2-it));
        p_dx_LAM         = p.prop_urine*(p.Fuji_sn(ih_test)*LAM_mat(ih_test,1,2) +  p.prop_urine*p.Alere_sn(ih_test)*LAM_mat(ih_test,1,1));
        p_dx_clinical    = (1-(p_dx_Xpert+p_dx_LAM))*p.clinical(2);
        p_dx_Xpert_extra = p.prop_sputum(2)*p.Xpert_sn(ih_test)*(2-it)*p.xpert_extra(1);       
        p_dx = 1 - ((1-p_dx_Xpert)*(1-p_dx_clinical)*(1-p_dx_LAM)*(1-p_dx_Xpert_extra));
        p_Tx_link_inp = p.prop_symp(2)*p.Tx_init(2)*p_dx;
        
              
        inds    = sub2ind([i.nstates, i.nstates], Tx_inp, AdmissionDx_inp);
        m(inds) = m(inds) + r.Dx*p_Tx_link_inp;
        
        inds    = sub2ind([i.nstates, i.nstates], E_inp, AdmissionDx_inp);
        m(inds) = m(inds) + r.Dx*(1-p_Tx_link_inp);
        
    end
end


% --- ART initiation upon hospitalisation 

hsets = {'low','mid','high'};
for ih = 1:length(hsets)
    hset         = hsets{ih};
    source_stubs = {s.(['h',hset]), s.(['ART_n_',hset])};
    
    for istub = 1:length(source_stubs)
        s_h   = source_stubs{istub};
        ART_n = s.(['ART_n_',hset,'_inp']);
        ART_n_com = s.(['ART_n_',hset]);
        
        % > First time initiating ART
        % Amongst TB negatives
        inds = sub2ind([i.nstates, i.nstates], intersect(uninfected_inp, s.ART_p_inp), intersect(uninfected_com, s_h));
        m(inds) = m(inds) + r.hosp(ih)*p.vs;
        inds = sub2ind([i.nstates, i.nstates], intersect(uninfected_inp, ART_n), intersect(uninfected_com, s_h));
        m(inds) = m(inds) + r.hosp(ih)*(1-p.vs);
        % Amongst pre-careseeking TB cases (I_com_sym and I_com_asym)
        inds = sub2ind([i.nstates, i.nstates], intersect(s.AdmissionDx_inp, s.ART_p_inp), intersect(s.I_com, s_h));
        m(inds) = m(inds) + r.hosp(ih)*p.vs;
        inds = sub2ind([i.nstates, i.nstates], intersect(s.AdmissionDx_inp, ART_n), intersect(s.I_com, s_h));
        m(inds) = m(inds) + r.hosp(ih)*(1-p.vs);
        % Amongst between-careseeking TB cases (E_com)
        inds = sub2ind([i.nstates, i.nstates], intersect(s.AdmissionDx_inp, s.ART_p_inp), intersect(s.E_com, s_h));
        m(inds) = m(inds) + r.hosp(ih)*p.vs;
        inds = sub2ind([i.nstates, i.nstates], intersect(s.AdmissionDx_inp, ART_n), intersect(s.E_com, s_h));
        m(inds) = m(inds) + r.hosp(ih)*(1-p.vs);
        % Amongst TB positives awaiting diagnosis (Dx_com)
        inds = sub2ind([i.nstates, i.nstates], intersect(s.AdmissionDx_inp, s.ART_p_inp), intersect(s.Dx_com, s_h));
        m(inds) = m(inds) + r.hosp(ih)*p.vs;
        inds = sub2ind([i.nstates, i.nstates], intersect(s.AdmissionDx_inp, ART_n), intersect(s.Dx_com, s_h));
        m(inds) = m(inds) + r.hosp(ih)*(1-p.vs);
        
        % > Re-initiating ART
        % Amongst TB negatives
        inds = sub2ind([i.nstates, i.nstates], intersect(uninfected_inp, s.ART_p_inp), intersect(uninfected_com, ART_n_com));
        m(inds) = m(inds) + r.hosp(ih)*p.vs;
        inds = sub2ind([i.nstates, i.nstates], intersect(uninfected_inp, ART_n), intersect(uninfected_com, ART_n_com));
        m(inds) = m(inds) + r.hosp(ih)*(1-p.vs);
        % Amongst pre-careseeking TB cases (I_com_sym and I_com_asym)
        inds = sub2ind([i.nstates, i.nstates], intersect(s.AdmissionDx_inp, s.ART_p_inp), intersect(s.I_com, ART_n_com));
        m(inds) = m(inds) + r.hosp(ih)*p.vs;
        inds = sub2ind([i.nstates, i.nstates], intersect(s.AdmissionDx_inp, ART_n), intersect(s.I_com, ART_n_com));
        m(inds) = m(inds) + r.hosp(ih)*(1-p.vs);
        % Amongst between-careseeking TB cases (E_com)
        inds = sub2ind([i.nstates, i.nstates], intersect(s.AdmissionDx_inp, s.ART_p_inp), intersect(s.E_com, ART_n_com));
        m(inds) = m(inds) + r.hosp(ih)*p.vs;
        inds = sub2ind([i.nstates, i.nstates], intersect(s.AdmissionDx_inp, ART_n), intersect(s.E_com, ART_n_com));
        m(inds) = m(inds) + r.hosp(ih)*(1-p.vs);
        % Amongst TB positives awaiting diagnosis (Dx_com)
        inds = sub2ind([i.nstates, i.nstates], intersect(s.AdmissionDx_inp, s.ART_p_inp), intersect(s.Dx_com, ART_n_com));
        m(inds) = m(inds) + r.hosp(ih)*p.vs;
        inds = sub2ind([i.nstates, i.nstates], intersect(s.AdmissionDx_inp, ART_n), intersect(s.Dx_com, ART_n_com));
        m(inds) = m(inds) + r.hosp(ih)*(1-p.vs);
        
    end
end


% --- Hospital discharge  

inds = sub2ind([i.nstates, i.nstates], intersect(community, s.ART_p), intersect(inpatient,s.ART_p_inp));
m(inds) = m(inds) + r.inp_leave(1);
inds = sub2ind([i.nstates, i.nstates], intersect(community, s.ART_n_low), intersect(inpatient,s.ART_n_low_inp));
m(inds) = m(inds) + r.inp_leave(1);
inds = sub2ind([i.nstates, i.nstates], intersect(community, s.ART_n_mid), intersect(inpatient,s.ART_n_mid_inp));
m(inds) = m(inds) + r.inp_leave(2);
inds = sub2ind([i.nstates, i.nstates], intersect(community, s.ART_n_high), intersect(inpatient,s.ART_n_high_inp));
m(inds) = m(inds) + r.inp_leave(3);



% -------------------------------------------------------------------------
% --- DEVELOPING HIV AND CD4 PROGRESSION ----------------------------------

% People developing HIV (h0 to hlow)
mHIVinc = zeros(i.nstates);
inds = sub2ind([i.nstates, i.nstates], s.hlow, s.h0);
mHIVinc(inds) = mHIVinc(inds) + r.hiv;
M.HIV_inc = sparse(mHIVinc - diag(sum(mHIVinc,1)));

% >>> CD4 progression amongst those not on ART
inds = sub2ind([i.nstates, i.nstates], s.hmid, s.hlow);
m(inds) = m(inds) + r.CD4prog(1);

inds = sub2ind([i.nstates, i.nstates], s.hhigh, s.hmid);
m(inds) = m(inds) + r.CD4prog(2);

% >>> CD4 progression amongst those on ART, not suppressed in community
inds = sub2ind([i.nstates, i.nstates], s.ART_n_mid, s.ART_n_low);
m(inds) = m(inds) + r.CD4prog(1);

inds = sub2ind([i.nstates, i.nstates], s.ART_n_high, s.ART_n_mid);
m(inds) = m(inds) + r.CD4prog(2);

% >>> CD4 progression amongst inpatients
inds = sub2ind([i.nstates, i.nstates], s.ART_n_mid_inp, s.ART_n_low_inp);
m(inds) = m(inds) + r.CD4prog(1);

inds = sub2ind([i.nstates, i.nstates], s.ART_n_high_inp, s.ART_n_mid_inp);
m(inds) = m(inds) + r.CD4prog(2);



% Collate all the rates
M.lin = sparse(m - diag(sum(m,1)));








% -------------------------------------------------------------------------
% --- Nonlinear rates -----------------------------------------------------

%%%%% COMMUNITY MODEL 

% --- Allocating transitions
m = zeros(i.nstates);

for ih = 1:length(gps.hivs_com)
    hivs_com = gps.hivs_com{ih};
    
    U_com = i.U_com.(hivs_com);                            L_com = i.L_com.(hivs_com);
    Ip_com = i.I_com.(hivs_com).ptb;                       Ie_com = i.I_com.(hivs_com).etb;
    SelfCure_p_com = i.SelfCure_com.(hivs_com).ptb;        SelfCure_e_com = i.SelfCure_com.(hivs_com).etb;
    TreatCure_p_com = i.TreatCure_com.(hivs_com).ptb;      TreatCure_e_com = i.TreatCure_com.(hivs_com).etb;
    
    m(Ip_com, [U_com, L_com, SelfCure_p_com, SelfCure_e_com, TreatCure_p_com, TreatCure_e_com]) = p.Fast(ih)*p.ptb(ih);
    m(Ie_com, [U_com, L_com, SelfCure_p_com, SelfCure_e_com, TreatCure_p_com, TreatCure_e_com]) = p.Fast(ih)*(1-p.ptb(ih));
    m(L_com,  [U_com, L_com, SelfCure_p_com, SelfCure_e_com, TreatCure_p_com, TreatCure_e_com]) = 1-p.Fast(ih);
    
    m(:,[s.L_com, s.SelfCure_com, s.TreatCure_com]) = m(:,[s.L_com, s.SelfCure_com, s.TreatCure_com])*p.imm;  % Immune protection
    
end

%%%%% INPATIENT MODEL 

% --- Allocating transitions
for ip = 1:length(gps.hivs_inp)
    hivs_inp = gps.hivs_inp{ip};
    
    U_inp = i.U_inp.(hivs_inp);                            L_inp = i.L_inp.(hivs_inp);
    Ip_inp = i.I_inp.(hivs_inp).ptb;                       Ie_inp = i.I_inp.(hivs_inp).etb;
    SelfCure_p_inp = i.SelfCure_inp.(hivs_inp).ptb;        SelfCure_e_inp = i.SelfCure_inp.(hivs_inp).etb;
    TreatCure_p_inp = i.TreatCure_inp.(hivs_inp).ptb;      TreatCure_e_inp = i.TreatCure_inp.(hivs_inp).etb;
    
    m(Ip_inp, [U_inp, L_inp, SelfCure_p_inp, SelfCure_e_inp, TreatCure_p_inp, TreatCure_e_inp]) = p.Fast_inp(ip)*p.ptb_inp(ip);
    m(Ie_inp, [U_inp, L_inp, SelfCure_p_inp, SelfCure_e_inp, TreatCure_p_inp, TreatCure_e_inp]) = p.Fast_inp(ip)*(1-p.ptb_inp(ip));
    m(L_inp,  [U_inp, L_inp, SelfCure_p_inp, SelfCure_e_inp, TreatCure_p_inp, TreatCure_e_inp]) = 1-p.Fast_inp(ip);
    
    m(:,[s.L_inp, s.SelfCure_inp, s.TreatCure_inp]) = m(:,[s.L_inp, s.SelfCure_inp, s.TreatCure_inp])*p.imm;  % Immune protection
    
end

m(:,s.HIV) = m(:,s.HIV)*p.HIV_suscep;

M.nlin = sparse(m - diag(sum(m,1)));


% --- Getting force-of-infection
m = zeros(1,i.nstates);
m(s.infectious) = r.beta;
m(intersect(s.infectious, s.HIV_nvs)) = m(intersect(s.infectious, s.HIV_nvs))*p.rel;
M.lambda = sparse(m);


% --- Get the mortality rates
vec = r.mort*ones(i.nstates,1);

% TB-associated mortality risk
vec(s.disease_overall_HIVn) = r.mort_TB;

% Excess mortality risk associated with HIV
inds      = unique([s.hlow, s.ART_n_low, s.ART_n_low_inp]);
vec(inds) = vec(inds) + r.mort_HIV + r.mort_HIV_excess_low;

inds      = unique([s.hmid, s.ART_n_mid, s.ART_n_mid_inp]);
vec(inds) = vec(inds) + r.mort_HIV + r.mort_HIV_excess_mid;

inds      = unique([s.hhigh, s.ART_n_high, s.ART_n_high_inp]);
vec(inds) = vec(inds) + r.mort_HIV + r.mort_HIV_excess_high;

inds      = unique([s.ART_p, s.ART_p_inp]);
vec(inds) = vec(inds) + r.mort_HIV + r.mort_HIV_excess_ART;

vec(s.disease_inp_ART)  = vec(s.disease_inp_ART)  + r.mort_inp;
vec(s.disease_inp_low)  = vec(s.disease_inp_low)  + r.mort_inp;
vec(s.disease_inp_mid)  = vec(s.disease_inp_mid)  + r.mort_inp;
vec(s.disease_inp_high) = vec(s.disease_inp_high) + r.mort_inp;


% Now organise into total / HIV-TB / HIV+TB
m = zeros(i.nstates,3);
m(:,1) = vec;

inds = s.disease_overall_HIVn;
m(inds,2) = vec(inds);

inds = intersect(s.disease_overall,s.HIV);
m(inds,3) = vec(inds);

M.mortvec = m;