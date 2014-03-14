function out=biased_gard(p, f, biasedmethod)
% out=biased_gard_v10(p, f, biasedmethod)
% Will run a biased GARD, towards the most frequent compotype.
% Defaults: f=1.1, biasedmethod=1.
% 16/06/2011 GARD10, by Omer Markovitch
% 19/09/11 some debugging

if ~exist('p', 'var') || isempty(p); p=tgs_agard; end;
if ~exist('f', 'var') || isempty(f); f=1.1; end;
if ~exist('biasedmethod', 'var') || isempty(biasedmethod); biasedmethod=1; end;

regular=tgs_agard_v10(p, 1);
[n, x]=hist(regular.tags, [0:1:size(regular.comps,2)]);
n(1)=[]; x(1)=[]; %don't include the drifts
it=x(n==max(n)); %index of the most frequent compotype, the target
trgt=regular.comps(:,it); trgt=trgt/sum(trgt)*p.NG*p.splitsize;
targetfreq=size(find(regular.tags==it), 1); %original frequency of the target
drift=size(find(regular.tags==0), 1); %original freq of the drift

biased=tgs_agard_v10(p, 1, f, biasedmethod, trgt);
h=tgs_H(trgt, biased.comps);
itbiased=find(h==max(h));
targetfreqbiased=size(find(biased.tags==itbiased), 1);
driftbiased=size(find(biased.tags==0), 1);

out=[];
out.p=p;
out.factor=f;
out.biasedmethod=biasedmethod;
out.regular=regular;
out.it=it;
out.target=trgt;
out.targetfreq=zerosih(targetfreq);
out.drift=zerosih(drift);
out.biased=biased;
out.h=h;
out.itbiased=itbiased;
out.targetfreqbiased=zerosih(targetfreqbiased);
out.driftbiased=zerosih(driftbiased);

out=orderfields(out);
return;

function z=zerosih(value)

if value==0; z=value+1e-8;
else z=value; end;
	
return;