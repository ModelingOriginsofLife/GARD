function [index, times, freq]=getcomposometime_v10(tags)
% [index, times]=getcomposometime_v10(tags)
% times(:,1)=from, times(:,2)=to, times(:,3)=duration
% Will focus only on the most frequent compotype.
% 19/09/2011 GARD10, by Omer Markovitch

%Make sure the input is OK
if ~exist('tags','var') || isempty(tags) ; error('tags'); end;

gens=length(tags);
n=hist(tags, [0:1:10]);
n(1)=[];
index=find(n==max(n));
index=index(1);
freq=n(index);

tags(find(tags~=index))=0;
tags(end)=0;
times=[];
manytimes=0;

for from=1:gens;
	
	if (tags(from)>0);
		
		for to=from+1:gens;
			
			if tags(to)==0;
				manytimes=manytimes+1;
				times(manytimes,1)=from;
				times(manytimes,2)=to-1;
				times(manytimes,3)=to-from;
				from=to+1;
				break;
			end;
			
		end;
		
	end;
	
end;

return;


