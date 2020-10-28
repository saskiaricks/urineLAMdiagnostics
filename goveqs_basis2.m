function [out, lam] = goveqs_basis2(t, in, M, i, s, r, p, sel, agg)

invec = in(1:i.nstates);

% Normalise by populations
lam = M.lambda*invec;
allmat = M.lin + lam*M.nlin;
out = allmat*invec;

% Implement deaths
morts = M.mortvec(:,1).*invec;
out = out - morts;

% Implement births
births = sum(morts);
out(i.U_com.h0) = out(i.U_com.h0)+births;

% Get the auxiliaries
out(i.aux.inc_total)       = agg.inc_total*(sel.inc_total.*allmat)*invec;
out(i.aux.mortTB)          = sum(M.mortvec(:,[2,3]).*[invec, invec],1); 

out(i.aux.HIVcom_high)     = agg.HIVcom*(sel.HIVcom_high.*allmat)*invec;
out(i.aux.HIVcom_mid)      = agg.HIVcom*(sel.HIVcom_mid.*allmat)*invec;
out(i.aux.HIVcom_low)      = agg.HIVcom*(sel.HIVcom_low.*allmat)*invec;
out(i.aux.HIVinp_high)     = agg.HIVinp*(sel.HIVinp_high.*allmat)*invec;
out(i.aux.HIVinp_mid)      = agg.HIVinp*(sel.HIVinp_mid.*allmat)*invec;
out(i.aux.HIVinp_low)      = agg.HIVinp*(sel.HIVinp_low.*allmat)*invec;

out(i.aux.hospitalisation) = sum((sel.hospitalisation.*allmat)*invec);
out(i.aux.Tx_inits)        = agg.Tx_inits*(sel.Tx_inits.*allmat)*invec;
out(i.aux.morts_hhigh)     = sum(morts(s.disease_inp_high)); 
out(i.aux.morts_inp)       = sum(morts(s.disease_inp));






