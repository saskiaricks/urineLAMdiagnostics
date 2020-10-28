plotting = 0;

[num, txt, raw] = xlsread('Data/HIV New Infections_National_new.xlsx','02-HIV incidence per 1000 pop');

% --- Get the HIV incidence data

[a,b] = find(strcmp(raw,'All ages estimate'));  HIV_inc = cell2mat(raw(a+1:end,b))/1e3;
[a,b] = find(strcmp(raw,'Year'));              year_inc = cell2mat(raw(a+1:end,b));

year_inc = [year_inc(1)-1; year_inc];
HIV_inc  = [0; HIV_inc];

% Normalise so level in 2018 is 1
ind = find(year_inc==2018);
HIV_inc = HIV_inc/HIV_inc(ind);

if plotting
    figure;
    subplot(1,2,1); plot(year_inc, HIV_inc/HIV_inc(end));
    title('HIV incidence');
end
HIVdat.inc = [year_inc, HIV_inc/HIV_inc(end)]';

% --- ART coverage data - actual rates will be estimated later, here just getting the shape of the curve
[num, txt, raw] = xlsread('Data/HIV Treatment cascade_National_new.xlsx','02-Coverage of people receivi'); 
[a,b] = find(strcmp(raw,'All ages estimate'));       pART = cell2mat(raw(a+1:end,b));
[a,b] = find(strcmp(raw,'Year'));               year_pART = cell2mat(raw(a+1:end,b));

% Normalise so level in 2018 is 1
ind = find(year_pART==2018);
pART = pART/pART(ind);
param = sigm_fit(year_pART, pART, [0 NaN NaN NaN],[],0);

if plotting
    subplot(1,2,2); hold on; plot(year_pART, pART);
    fsigm = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)));
    plot(year_pART, fsigm(param, year_pART));
    title('Proportion on ART');
    legend('Data','Fit');
end
HIVdat.pART_prms = param;