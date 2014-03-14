function [childA, childB, parent]=tgs_split_v10(parent, p)
% [childA, childB, parent]=tgs_split_v10(parent, p)
% Will stochastically split the parent assembly into two children
% 16/06/2011 GARD10, by Omer Markovitch
% 19/09/11 some debugging

if ~exist('parent','var') || isempty(parent) || isempty(find(parent>0, 1)); error('Problem with parent'); end;

ng=length(parent);
childA=zeros(ng, 1);

for rs=1:floor(sum(parent)/2);
	species=tgs_rndpdf(parent);
	childA(species)=childA(species)+1;
	parent(species)=parent(species)-1;
end;

while (sum(childA)<sum(parent));
	species=tgs_rndpdf(parent);
	parent(species)=parent(species)-1;
end;

childB=parent;

return;
