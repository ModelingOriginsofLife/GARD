function [tags, comps] = tgs_clust(samples, optclst, hthreshold, clusternorm);
% TGS_CLUST - A trace clusterer based on k-means.
%

%   Copyright 2004-2005 Barak Shenhav & Weizmann Institute of Science
%
%   20050801 v0.001 - First incorporation The GARD Suite. Based on
%                     standalone code from 2004
%   20050803 v0.003 - Fix output in case of a single composome from row to
%                     colomn
%   20050907 v0.005 - Fix to avoid clustering with number of clusters larger
%                     then number of samples that path the H threshold

%samples==trace
%optclst==number of clusters - a vector containing a number of groups to be
%tested [2,3,4,5...]- divides into 2, 3 , 4 , 5  groups.
%hthreshold == the cutoff for drift/nondrift

if ~exist('clusternorm', 'var') || isempty(clusternorm); clusternorm='cosine'; end; %omer 27/05/10

p.hthresh=hthreshold;
p.ks=optclst;
[out] = tgs_acluster(samples,p, clusternorm);
tags=out.tags;
comps=out.comps;