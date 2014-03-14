function pmc=pmc_beta_v10(beta)
% pmc=pmc_beta_v10(beta)
% pmc>1 is mutual-catalysis excess
% 16/06/2011 GARD10, by Omer Markovitch

if ~exist('beta', 'var') || isempty(beta); error('Pmc: problem with beta'); end;
[ng, bla]=size(beta);
if (ng~=bla); error('Pmc: beta is not NG*NG matrix'); end;
pmc=0;
for i=1:ng; pmc=pmc+beta(i,i); end;
pmc=sum(beta(:))/(pmc*ng);

return;
