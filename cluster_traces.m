function [comps, compnum, compcounts] = cluster_traces(trace)
% %% cluster_traces - clusters a group of vectors using k-means, using
% different k's and selecting the best result by the silhouette score.
% returns all the clustering results 
%      
%
%   Input parameters
%   ----------------+------------+------------------------------------------------------------------------------------------------------
%  | Name           |   Type     | Description
%   ----------------+------------+------------------------------------------------------------------------------------------------------
%   trace
%
%   variables returned: 
%   ----------------+------------+------------------------------------------------------------------------------------------------------
%  | Name           | Type       | Description
%   ----------------+------------+------------------------------------------------------------------------------------------------------
%   comps
%   compnum 
%   compcounts
%
%   functions used
%   ----------------+------------+------------------------------------------------------------------------------------------------------
%  | Name           | Type       | Description
%   ----------------+------------+------------------------------------------------------------------------------------------------------
%   tgs_clust
%
% (c) Aia Oz 4/2006
%  modified and documented 5/2007
if nargin<1
    index=0;
end

%comps - the center of the cluster (the average composition).
%compnum is the number of clusters chosen
%compcounts is how many splits in each of these clusters

    [tags, comps]= tgs_clust(trace, [1,2,3,4,5,6,7,8,9,10], 0.95); %receives the trace, the number of groups to be divided into and the cutoff for H - homology.
    
    compnum = size(comps, 2)
    
    compcounts = (hist(tags(tags >0), compnum))';
  
 return; 
    
    