function out=tgs_agard_v10(p, alsocluster, biasfactor, biasmethod, biascomposition)
% out=tgs_agard_v10(p, alsocluster, biasfactor, biasmethod, biascomposition)
% 02/06/2011 [dd/mm/yy] A-GARD (Amphiphile GARD) version 10, by Omer Markovitch
% a very short example for running: p=tgs_parameters_v10; p.seed=1; o=tgs_agard_v10(p, 1); c=tgs_carpet_v10(o.trace);
% contact information: omermar@gmail.com , http://sites.google.com/site/omermar
% 19/09/11 some debugging

if ~exist('p', 'var') || isempty(p); p=tgs_parameters; end;
if ~exist('alsocluster', 'var') || isempty(alsocluster); alsocluster=0; end;
if ~exist('biasmethod', 'var') || isempty(biasmethod); biasmethod=0; end;

if exist('biasfactor', 'var');
	if isempty(biasfactor); biasfactor=1.1; end;
	if ~exist('biascomposition', 'var'); error('Must input a target composition'); end;
	if ~exist('biasmethod', 'var'); biasmethod=1; end;
end;

n=p.n; %the assembly
NG=p.NG; %repertior size
gen=p.gen; %how long to run a simulation
nmax=ceil(p.splitsize*NG); %the size at which the assembly splits
nmin=floor(nmax/2);

if length(p.rho)~=NG;
	rho=ones(NG, 1)/NG;
	p.rho=rho;
else
	rho=p.rho;
end

%generate a beta matrix if not given
if isempty(p.Beta) || length(p.Beta)==1;
	Beta = tgs_newbeta_v10(p);
	p.Beta=Beta;
else
	Beta=p.Beta;
end;

%reset the random number generator
if length(p.seed)==1; p.seed(2:3)=p.seed(1); end;
state = randn('state');
randn('state', p.seed(2));

%initial random assembly
if isempty(n) || sum(n)<=0 || length(n)~=p.NG
    n = histc(rand(nmin ,1)*NG, 0:NG);
	if length(n)~=NG; n=n(1:NG); end;
	if size(n,2)>1; n=n'; end;
	p.n=n;
    if p.seed(2)~=-1; rand('seed', p.seed(2)); randn('seed', p.seed(2)); end;
end;

trace=zeros(NG, gen); %a matrix holding the compositions of all the generations
counter=1;
step=0;
sum_n=sum(n); %current size of assembly
newbeta=zeros(NG, NG); %for biased gard

while counter <= gen,
	
	if sum_n >= nmax, %if we are at the end of the growth cycle (Nmax)
        trace(:, counter)=n;
		time(counter)=1/dt;
        counter=counter+1;
        step=0;
		n=tgs_split_v10(n, p);
        sum_n=nmin;
		
		if biasmethod~=0; %for selection. note that H is calculated at the begining of the growth cycle (Nmin)
			newbeta=modifybeta(Beta, n, tgs_H(n, biascomposition)*biasfactor, biasmethod);
			p.Beta=newbeta;
		end;
		
        continue;
    elseif sum_n==0 %assembly can die, if all of its molecules left
        disp(['the assembly has died after ' num2str(counter) ' generations, R.I.P.']);
		trace(:, counter)=-1;
        break;
	end; %if - splitting

    step=step+1;
	[n, asshist{counter}, dt]=tgs_grow_v10(n, p);
	sum_n=sum(n);
	
end; %while - for each generation

out=[];
clusternorm='cosine'; if alsocluster==2; clusternorm='sqEuclidean'; end;
if alsocluster; out=tgs_acluster(trace,p, clusternorm); end; %do the clustering and find compotypes
out.n=n;
out.rho=rho;
out.Beta=Beta;
out.trace=trace;
out.parameters=p;
out.time=time;
out.clusternorm=clusternorm;
out.version=10.0;

if biasmethod~=0;
	out.bias=biasmethod;
	out.target=biascomposition;
	out.biasfactor=biasfactor;
end;

out=orderfields(out);
return;

function newbeta=modifybeta(oldbeta, composition, factor, method)
% Will modify the BETA matrix, according to the species that exists in the composition.
% method=1 (default) - only enhances Bij that both i & j are inside the assembly
% Bij is how i is catalyzed by j. Bij=beta[i,j]
% 02/06/2011 GARD10, by Omer Markovitch

newbeta=oldbeta;
species=find(composition>0);

switch (method),
	case 1, %enhance only Bij that both i&j are in the assembly
		
		for i=1:max(size(species));
			s=species(i);
			for j=1:max(size(species)); newbeta(s, species(j))=newbeta(s, species(j))*factor; end;
		end; %for i
				
	otherwise
		error('Unknown method');
end;

return;
