function [out,stats] = tgs_cluster(trace,p, clusternorm);

%function [tags, comps,angles,silhvec,stats] = tgs_clust(trace,p);
%clusters a trace to extract the compotypes
%inputs:
%trace - NGxgens trace matrix
%p - parameter structure (tgs_parameters)
%outputs:
%a structure containing the following fields
%tags - a length gens vector specifying to which cluster each generation belongs
%comps - the compotypes/clusters
%stats - be careful if this is requested then execution time is greatly increased.
%stats provides the sum of distances in columns 1:k and the mean silhouettes in column k+1:2k for all replicas.
%---these are not output currently...
%angles - the H value between each generation and the one that preceeded it
%silvec - a vector of length # of k values to try indicating the sillhouette value for that k
%counts - the size in generations of each compotype

if ~exist('clusternorm', 'var') || isempty(clusternorm); clusternorm='cosine'; end; %omer 27/05/10

%get parameters
if nargin<2; p=[]; end
if ischar(p); p=load(p); end
driftsize=[]; if isfield(p,'driftsize'); driftsize=p.driftsize; end
s1=rand('state');s2=randn('state');

if nargout>1; stats=zeros(p.replicas,length(p.ks)*2); end

samples=trace;
optclst = p.ks;%optclst(optclst>p.NG)=[];
mink=p.mink; if isempty(mink) | mink==0; mink=Inf; end
hthreshold = p.hthresh;
replicas=p.replicas;
sils=zeros(length(p.ks));

if p.seed(3)~=-1;rand('seed',p.seed(3)); randn('seed',p.seed(3));end

%identify drift / non-drift
if isempty(driftsize);
    [indx,angles] = tgs_nondrift(samples, hthreshold); 
else
    [indx,angles] = tgs_nondrift(samples, hthreshold,p.driftsize); 
end    

tagmat = zeros(size(samples,2),length(optclst));
streak=0;
bestsil=-inf;
bestsilEuq=+inf; %omer

%test all requested values of k (# clusters)
for i=1:length(optclst),
    silh=-inf;
    flag=1;
    if sum(indx)==0; %all drift and no non-drift. single compotype will be center of gravity
        [tags comps silh] = tgs_kmeans(samples, 1,1, 'off', clusternorm); %omer 27/05/10
        tags=[];
        flag=0;
    end        
    if optclst(i) <= sum(indx), %don't request more clusters then points
        flag=0;
        if nargout>1;
            [tags, comps, silh,stat1] = tgs_kmeans(samples(:,indx), optclst(i),replicas, 'off', clusternorm); %omer 27/05/10
            stats(1:replicas,i)=stat1(:,1);stats(1:replicas,length(optclst)+i)=stat1(:,2);
        else
            [tags comps silh] = tgs_kmeans(samples(:,indx), optclst(i),replicas, 'off', clusternorm); %omer 27/05/10
        end
        if isempty(tags) flag=1; end %clustering was not succesful so don't crash!
    end
    if flag==1;%no clustering performed
        tagmat(indx,i) = NaN;
        compcell{i} = [];
        silhvec(i) = NaN;
    else
        %receives only the samples marked as composomes and clusters them. tag==the group.
        tagmat(indx,i) = tags;
        compcell{i} = comps';
        silhvec(i) = silh;
    end
    streak=streak+1;
	if silh>=bestsil; bestsil=silh; streak=0; end;
    if streak>=mink; break; end %don't go more then mink without seeing any improvement
end

%select number of clusters with best silhouette

[dummy compnum] = max(silhvec);
tags = tagmat(:,compnum);
comps = compcell{compnum};
if compnum == 1; comps = comps'; end

[counts,d]=hist(tags,max(tags+1)); counts=counts(2:end);

out.tags=tags;
out.comps=comps;
out.parameters=p; 
rand('state',s1);randn('state',s2);
return