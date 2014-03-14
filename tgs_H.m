function [h]=tgs_H(set1, set2);

%function [h]=tgs_H(set1, set2);
%calculates the H value
%(cosine of angle between two normalized vectors)
%higher values indicate more similarity
%inputs:
%set1, set2 - may be matrices of vectors 
%if set2 is omitted then set1 is compared with itself
%outputs:
%h
%the H value for every element in set one compared with every element in set two.

eps=0.0000001;
norm1 = set1 ./ (ones(size(set1,1),1) * max(sqrt(sum(set1.^2)),eps));
if nargin==1; 
    norm2=norm1;
else
norm2 = set2 ./ (ones(size(set2,1),1) * max(sqrt(sum(set2.^2)),eps));
end
if size(norm1,1)==1; norm1=norm1'; end
if size(norm2,1)==1; norm2=norm2'; end
h=norm1'*norm2;

h=max(h,0);
h=min(h,1);