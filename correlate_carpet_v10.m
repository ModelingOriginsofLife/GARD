function [deltat, ctnorm, ct, hits]=correlate_carpet_v10(trace, t0, disp, comment);
% [deltat, ctnorm, ct, hits]=correlate_carpet_v10(trace, t0, disp, comment);
% This funtion will calculate the autocorrelation for the H, based on the CARPET.
% Becuase the CARPET is already the correlation between generations then all that is left is to calculate the average H (<H>).
% 16/06/2011 GARD10, by Omer Markovitch

%With t0 we can select specific generations to calculate for
if ~exist('t0','var') | isempty(t0); t0=[1:size(trace,2)]; end

if ~exist('disp','var') | isempty(disp); disp=0; end
if ~exist('comment','var') | isempty(comment); comment='Autocorrelation'; end

carpet=tgs_carpet(trace,'none');
newcarpet=tril(carpet);
numgens=size(newcarpet,2);
newcarpet=flipud(newcarpet);
vec=[1:numgens]';
mat=repmat(vec,1,numgens);
temp=tril(mat-mat'+1);
newc=temp+triu(numgens+2-temp')-eye(numgens)*(numgens+1);
newc=newc+repmat([0:numgens-1],numgens,1)*numgens;
newcarpet=flipud(newcarpet(newc(:,t0)));
ct=sum(newcarpet,2);
hits=sum(newcarpet>0,2);
ctnorm=ct./(hits+0.00000000001);
deltat=[0:size(carpet,1)-1]';

if disp~=0
  plot(deltat,ctnorm);
  xlabel('Generation difference');
  ylabel('c(t)');
  title(comment);
end

return

