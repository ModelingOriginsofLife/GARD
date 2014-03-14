function out=population_gard_nmin_v10(p, numsplits, trgt)
% out=population_gard_nmin_v10(p, numsplits, trgt)
% Will perform a population dynamics on the population of assemblies, based on Moran process.
% numsplits is for how many splits to run the dynamics
% 19/09/2011 GARD10, by Omer Markovitch

if ~exist('p', 'var') || isempty(p); error('Must input parameters structure'); end;
if ~exist('numsplits', 'var') || isempty(numsplits) || numsplits<1; numsplits=p.gen; end;
if length(p.seed)==1; p.seed(2:3)=p.seed(1); end;
if size(p.seed,2)>1; rand('seed', p.seed(2)); randn('seed', p.seed(2)); end;
if isempty(p.Beta); p.Beta=tgs_newbeta_v10(p); end;
if isempty(p.rho); p.rho=ones(p.NG, 1)/p.NG; end;

popsize=p.gen;
initialpop=zeros(p.NG, popsize);
currpop=zeros(p.NG, popsize);
Nmax=ceil(p.splitsize*p.NG);
gb=zeros(popsize,1);
splitorder=zeros(numsplits,2); %[,1]=the index of the parent who split, [,2]=the index of the second assembly who was replaced
s=0;

% 13/07/11 Rafi suggested that a GARD trace is already biased towards composomes
for i=1:popsize;
	n = histc(rand(Nmax ,1)*p.NG, 0:p.NG);
	if length(n)~=p.NG; n=n(1:p.NG); end;
	if size(n,2)>1; n=n'; end;
	initialpop(:,i)=n;
end; %for i - generate initial random population

currsize=sum(initialpop);
i=find(currsize>=Nmax);
for j=1:max(size(i));
	childA=tgs_split_v10(initialpop(:,i(j)), p);
	currpop(:, i(j))=childA;
end; %assemblies need splitting before we begin

rates=zeros(2*p.NG, popsize); %is times 2 for join & leave
tmpbeta=zeros(p.NG, p.NG);
for i=1:popsize;
	if exist('trgt', 'var') && ~isempty(trgt);
		tmpbeta=p.Beta;
		gb(i)=tgs_H(currpop(:,i), trgt)*1.1;
		p.Beta=modifybeta(p.Beta, currpop(:,i), gb(i), 1);
		rates=modifyrates(p, rates, i, currpop(:,i));
		p.Beta=tmpbeta;
	else
		rates=modifyrates(p, rates, i, currpop(:,i));
	end;
end;

while s<numsplits; 
	
	currsize=sum(currpop);
	i=find(currsize>=Nmax);
	
 	if length(i)>1; error('more than 1 assembly reached split size?'); end;
	
	if ~isempty(i);
		j=randi([1 popsize], 2, 1);
		if (j(1)==i); j=j(2); else j=j(1); end;
		[childA, childB]=tgs_split_v10(currpop(:,i), p);
		currpop(:, i)=childA;
		currpop(:, j)=childB;
				
		if exist('trgt', 'var') && ~isempty(trgt);
			tmpbeta=p.Beta;
			gb(i)=tgs_H(childA, trgt)*1.1;
			p.Beta=modifybeta(p.Beta, childA, gb(i), 1);
			rates=modifyrates(p, rates, i, childA);
			p.Beta=tmpbeta;
			gb(j)=tgs_H(childB, trgt)*1.1;
			p.Beta=modifybeta(p.Beta, childB, gb(j), 1);
			rates=modifyrates(p, rates, j, childB);
			p.Beta=tmpbeta;
		else
			rates=modifyrates(p, rates, i, childA);
			rates=modifyrates(p, rates, j, childB);
		end;
		
		s=s+1; %if ( (s/100)==round(s/100) ); s, end;
		splitorder(s,:)=[i,j];
		if (s>=numsplits); break; end;
	else %if - split
		mu=tgs_rndpdf(rates(:));
		muass=floor((mu-1)/(2*p.NG));
		mung=mu-muass*2*p.NG;
		muass=muass+1;
		
		if (mung<=p.NG);
			currpop(mung, muass)=currpop(mung, muass)+1; %join
		else

			currpop(mung-p.NG, muass)=currpop(mung-p.NG, muass)-1; %leave
			
			if ( sum( currpop(:,muass) ) )<=0;
				n = histc(rand(Nmax ,1)*p.NG, 0:p.NG);
				if length(n)~=p.NG; n=n(1:p.NG); end;
				if size(n,2)>1; n=n'; end;
				currpop(:,muass)=n;
			end; %if - an assembly died
			
		end;
		
		if exist('trgt', 'var') && ~isempty(trgt);
			tmpbeta=p.Beta;
% 			gb(muass)=tgs_H(currpop(:, muass), trgt)*1.1; %
			p.Beta=modifybeta(p.Beta, currpop(:, muass), gb(muass), 1);
			rates=modifyrates(p, rates, muass, currpop(:, muass));
			p.Beta=tmpbeta;
		else
			rates=modifyrates(p, rates, muass, currpop(:, muass));
		end;
		
	end;
	
end; %while s - splits

out.p=p;
out.trace=currpop;
out.splitorder=splitorder;
if exist('trgt', 'var') && ~isempty(trgt); out.target=trgt; end;
out=orderfields(out);
return;

function newrates=modifyrates(p, oldrates, index, assembly)
% Will modify the 'rates' matrix at location 'index' according to the composition of 'assembly'

newrates=oldrates;
kfrho=p.Kf*p.rho;
sum_n=sum(assembly);
bn=( 1+1/sum_n*( p.Beta*assembly ) );
newrates(:,index)=[ (kfrho*sum_n).*bn ; ( p.Kb*assembly ).*bn ]; %this is the basic gard equation

return;

