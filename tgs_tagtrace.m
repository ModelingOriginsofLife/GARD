function [optidx, optcentroids, meansilh] = tgs_tagtrace(samples, clustnum, displaymode);
% TGS_TAGTRACE

replicas=8;
[optidx, optcentroids, meansilh] = tgs_kmeans(samples, clustnum,replicas,displaymode);
