function [newassembly, asshist, dt]=tgs_grow_v10(oldassembly, p);
% [newassembly, species, dt, history]=grow_v10(oldassembly, p)
% Will grow an assembly until it reaches the split size. species holds a list of species joined (or left) until splitsize reached.
% 16/06/2011 GARD10, by Omer Markovitch

nmax=ceil(p.NG*p.splitsize);
newassembly=oldassembly;
sum_n=sum(newassembly);
kfrho=p.Kf*p.rho;
dt=0;
asshist=[];

while (sum_n<nmax);
	bn=(1+1/sum_n*(p.Beta*newassembly));
	rates=[ (kfrho*sum_n).*bn ; (p.Kb*newassembly).*bn ]; %this is the basic gard equation
	[mu sum_rates]=tgs_rndpdf(rates);
	dt=dt+sum_rates;
    
	if mu <= p.NG, %join
        newassembly(mu) = newassembly(mu) + 1;
        sum_n=sum_n+1;
		asshist(end+1)=mu;
	else %leave
        mu = mu - p.NG;
        newassembly(mu) = newassembly(mu) - 1;
        sum_n=sum_n-1;
		asshist(end+1)=-1*mu;
        if newassembly(mu)<0; error('Invalid reaction, program will be terminated'); end;
	end;
	
	if (sum_n==0); return; end; %omer 22/05/2010 assembly has died

end;

return;

