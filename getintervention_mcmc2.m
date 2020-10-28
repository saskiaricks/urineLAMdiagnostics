clear all; load('estims.mat');
load('MCMC_res_MAIN2c.mat');

popn = 57e6;

%remove burn-in and thin
xs = xsto(10000:100:end,:);

for xt = 1:size(xs,1)    
    
    x = xs(xt,:);
    
    [p, r] = alloc_parameters2(x, p, r, xi);
  
  p1 = p; r1 = r; prm1 = prm; 
  
  
    % --- Get initial conditions
    [out, aux] = get_objective2(x, prm, ref, sel, agg, gps, lsq, HIVdat, LAM_mat);  

      
    init = aux.soln(end,1:end-1);
    
    % Current state (baseline)
    M0 = make_model_wEPTB(p1, r1, i, s, gps, LAM_mat);
    
    %------------- SCENARIO 1 --------------------------------------------
     % Scenario 1a1 - Inpatients (Alere)
     LAM_mat1a1 = zeros(8,3,2);
     LAM_mat1a1([3:4,7:8],1,1) = 1;
     M1a1 = make_model_wEPTB(p1, r1, i, s, gps, LAM_mat1a1);
    
    
    % Scenario 1a2 - Inpatients (Fuji)
    LAM_mat1a2 = zeros(8,3,2);
    LAM_mat1a2([3:4,7:8],1,2) = 1;
    M1a2 = make_model_wEPTB(p1, r1, i, s, gps, LAM_mat1a2);
    
    
    % Scenario 1b1 - Outpatients (Alere)
    LAM_mat1b1 = zeros(8,3,2);
    LAM_mat1b1([3:4,7:8],1,1) = 1;
    LAM_mat1b1([4,8],2,1) = 1;
    M1b1 = make_model_wEPTB(p1, r1, i, s, gps, LAM_mat1b1);
    
    
    % Scenario 1b2 - Outpatients  (Fuji)
    LAM_mat1b2 = zeros(8,3,2);
    LAM_mat1b2([3:4,7:8],1,2) = 1;
    LAM_mat1b2([4,8],2,2) = 1;
    M1b2 = make_model_wEPTB(p1, r1, i, s, gps, LAM_mat1b2);
    
    
    %Scenario 1c1 - Routine TB care (Fuji) - everyone
    LAM_mat1c1 = zeros(8,3,2);
    LAM_mat1c1([3:4,7:8],1,2) = 1;
    LAM_mat1c1([4,8],2,2) = 1;
    LAM_mat1c1([1:8],3,2) = 1;
    M1c1 = make_model_wEPTB(p1, r1, i, s, gps, LAM_mat1c1);
    
    
    %------------- SCENARIO 2 --------------------------------------------
     
     % SCENARIO 2 baseline (expanded Xpert Inpatients/Outpatients/Routine)
     p2a = p1;
     p2a.xpert_extra = [1 1 1];
     M2a = make_model_wEPTB(p2a, r1, i, s, gps, LAM_mat);
 
  
     % SCENARIO 2a1 Inpatients (Alere)  
     p2a1 = p1;
     LAM_mat2a1 = zeros(8,3,2); 
     LAM_mat2a1([3:4,7:8],1,1) = 1;
     p2a1.xpert_extra = [1 1 1];
     M2a1 = make_model_wEPTB(p2a1, r1, i, s, gps, LAM_mat2a1);        
     
    
     % SCENARIO 2a2 Inpatients (Fuji)
     p2a2 = p1;   
     LAM_mat2a2 = zeros(8,3,2); 
     LAM_mat2a2([3:4,7:8],1,2) = 1;
     p2a2.xpert_extra = [1 1 1]; 
     M2a2 = make_model_wEPTB(p2a2, r1, i, s, gps, LAM_mat2a2);    
   
     
     % SCENARIO 2b1 Outpatients (Alere)
     p2b1 = p1;   
     LAM_mat2b1 = zeros(8,3,2); 
     LAM_mat2b1([3:4,7:8],1,1) = 1;
     LAM_mat2b1([4,8],2,1) = 1;
     p2b1.xpert_extra = [1 1 1]; 
     M2b1 = make_model_wEPTB(p2b1, r1, i, s, gps, LAM_mat2b1); 
    
     
     % SCENARIO 2b2 Outpatients (Fuji)  
     p2b2 = p1; 
     LAM_mat2b2 = zeros(8,3,2); 
     LAM_mat2b2([3:4,7:8],1,2) = 1;
     LAM_mat2b2([4,8],2,2) = 1;
     p2b2.xpert_extra = [1 1 1]; 
     M2b2 = make_model_wEPTB(p2b2, r1, i, s, gps, LAM_mat2b2);     
 
     
     % SCENARIO 2c1 Routine TB care (Fuji)   
     p2c1 = p1; 
     LAM_mat2c1 = zeros(8,3,2); 
     LAM_mat2c1([3:4,7:8],1,2) = 1;
     LAM_mat2c1([4,8],2,2) = 1;
     LAM_mat2c1([1:8],3,2) = 1;
     p2c1.xpert_extra = [1 1 1]; 
     M2c1 = make_model_wEPTB(p2c1, r1, i, s, gps, LAM_mat2c1);     
    
     
     
     
    getsol = @(Mfinal) ode15s(@(t,in) goveqs_scaleup(t, in, M0, Mfinal, [2019 2022], i, s, r, p, sel, agg, HIVdat), [2019:2035], init(end,:), odeset('NonNegative',[1:i.nstates],'AbsTol',1e-10,'RelTol',1e-10));
    models = {M0, M1a1, M1a2, M1b1, M1b2, M1c1, M2a, M2a1, M2a2, M2b1, M2b2, M2c1};
    allsol = [];
    
    for mi = 1:length(models)
        fprintf('%0.5g ', mi);
        [~,allsol(:,:,mi)] = getsol(models{mi});
    end
    fprintf('\n');
    dallsol = diff(allsol,[],1);
    
    % --- Incidence
    incs(:,:,xt)  = squeeze(dallsol(:,i.aux.inc_total(1),:))*1e5;
    incs_total(:,:,xt)  = squeeze(dallsol(:,i.aux.inc_total(1),:))*popn;

    
    % --- Mortality
    morts(:,:,xt) = squeeze(sum(dallsol(:,i.aux.mortTB,:),2))*popn;
    mortsinp(:,:,xt) = squeeze(dallsol(:,i.aux.morts_inp,:))*popn;
   
   
end

save('intv.mat','incs','incs_total','morts','mortsinp'); 





