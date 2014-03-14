function [ans, asympt]=guess_tau_v10(ctnorm, method);
% [ans, asympt]=guess_tau(ctnorm, method);
% Will guess when the decay of c(t) stops.
% It does so by simply looking at the first time c(t)<=c(t+1)
% Returning: ans=-1 for bad input or ans=0 if could not estimate, else ans>0
% Methods:
% (1) Guess only time by first rise. This is default.
% (2) Guess time by the first time c(t) is below the asymptotic value, by middle 50%.
% (3) Same as (2), instead of time returns the slope.
% 16/06/2011 GARD10, by Omer Markovitch

%Make sure the input is OK
if ~exist('ctnorm','var') | isempty(ctnorm) | length(ctnorm)<2; ans=-1; return; end;
if ~exist('method','var') | isempty(method); method=1; end;

ans=0;
asympt=0;
if method==1; ans=taubyrise(ctnorm); return; end;
if method==2; [ans, asympt]=timebyasympt(ctnorm); return; end;
if method==3; [ans, asympt]=slopebyasympt(ctnorm); return; end;

return;

function tau=taubyrise(ct);

tau=0;

for i=2:length(ct);
	if ct(i)>=ct(i-1);
		tau=i-1;
		return;
	end;
end

return;

function [time, avg]=timebyasympt(ct);

time=-1;
clength=length(ct);
from=round(0.25*clength);
to=round(0.75*clength);
avg=mean(ct(from:to));
tempavg=round(avg*1000);

for i=1:length(ct);
	if ct(i)*1000<=tempavg;
		time=i;
		return;
	end;
end

return;

function [slope, avg]=slopebyasympt(ct);

slope=-1;
clength=length(ct);
from=round(0.25*clength);
to=round(0.75*clength);
avg=mean(ct(from:to));
tempavg=round(avg*1000);

for i=1:length(ct);
	if ct(i)*1000<=tempavg;
		slope=(1-ct(i))/i;
		return;
	end;
end

return;
