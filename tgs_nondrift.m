function [indx,angles] = tgs_nondrift(trace, HThreshhold,driftsize);

%function [indx,angles] = tgs_nondrift(trace, HThreshhold,driftsize);
% separates the drift composomes from nondrift
%inputs:
%trace - NGxgens trace matrix
%Hthreshold - the H value threshold for determining whether a composome is a drifter
%driftsize - if empty then technique 1 is used if not then technique 2 is used
%technique 1 uses avg H value between previous and next gen
%technique 2 requires no H value less than the threshold between consecutive driftsize gens
%outputs
%indx - logical vector specifying whether each composome is drift (0) or non-drift (1)
%angles - the H value between each generation and the one that preceeded it

%use technique 2?
if nargin==3;
    [indx,angles] = tgs_marktrace2(trace, HThreshhold,driftsize);
    return;
end

if nargin < 1,
    error('Missing parameters');
end
if nargin < 2,
    HThrehhold = .9;
end

nspecies = size(trace,1);  % nspecies = number of generations  - the size of trace

trace = full(trace);      % wtf
osize=size(trace,2); q=sum(trace); q=min(find(q==0));
if ~isempty(q);
    trace=trace(:,1:q-1); %removed generations after the micelle died
end
trace = trace ./ (ones(nspecies,1) * (sum(trace.*trace).^.5)); 

avgH = sum(trace(:,1:end-1).*trace(:,2:end));
angles = [avgH(1) avgH]';
avgH = ([avgH(1) avgH] + [avgH avgH(end)])/2;
indx = (avgH > HThreshhold);
if ~isempty(q)
    indx(q:osize)=0;
    angles(q:osize)=0;
end
return

function [indx,angles] = tgs_marktrace2(trace, HThreshhold,driftsize);
if nargin < 1,
    error('Missing parameters');
end
if nargin < 2,
    HThrehhold = .9;
end

nspecies = size(trace,1);  % nspecies = number of generations  - the size of trace

trace = full(trace);      % wtf
osize=size(trace,2); q=sum(trace); q=min(find(q==0));
if ~isempty(q);
    trace=trace(:,1:q-1); %removed generations after the micelle died
end
trace = trace ./ (ones(nspecies,1) * (sum(trace.*trace).^.5)); 

avgH = sum(trace(:,1:end-1).*trace(:,2:end));
avgH = ([avgH(1) avgH]);
angles=avgH';
a = (avgH > HThreshhold);
b=[0 a(1:end-1); a; a(2:end) 0];
c=[find(b(1,:)==0 & b(2,:)==1) ;find(b(2,:)==1 & b(3,:)==0)];
d=c(2,:)-c(1,:)+1;
e=find(d>=driftsize);
indx=zeros(1,size(trace,2));
for x=1:length(e);
    indx(c(1,e(x))-1:c(2,e(x))-1)=1;
end

if ~isempty(q)
    indx(q:osize)=0;
    angles(q:osize)=0;
end
indx=logical(indx);

return