function [optidx, optcentroids, meansilh,stats] = tgs_kmeans(trace, knum, replicas,display, clusternorm);

%function [optidx, optcentroids, meansilh,stats] = tgs_kmeans(samples, knum, displaymode,replicas);
%the wrapper function for performing kmeans clustering
%inputs:
%trace - the NGxgens trace matrix
%knum - the number of clusters for this particular run
%replicas - number of replicas to perform
%display - optional to display k-means stats (not recommended but here for backwards compatibility)
%outputs:
%optidx - the cluster/compotype each generation/composome belongs to
%optcentroids - the clusters / compotypes
%meansilh - the mean silhouette value for the chosen clustering
%stats - be careful if this is requested then execution time is greatly increased.
%stats provides the sum of distances in column 1 and the mean silhouettes in column 2 for all replicas.

%if then next value is set to 1 then each replica calculates a silhouette and the max silhouette replica is chosen
%otherwise the min sum of distances to each composomes cluster is used to determine best replica (much faster)
usesil=0;

warning off MATLAB:divideByZero
warning off stats:kmeans:EmptyCluster
warning off stats:kmeans:FailedToConverge

if usesil;val=-Inf; else; val=Inf; end
if nargout==4; stats=zeros(replicas,2);end;
if nargin<4; display='off'; end
if ~exist('clusternorm', 'var') || isempty(clusternorm); clusternorm='cosine'; end; %omer 27/05/10

if knum > 1,
    repeats = replicas; 
    
    if nargin < 3,
        displaymode = 'off';
    end
    if nargin < 2,
        error('Missing parameters');
    end

    nspecies = size(trace,1);

    optd = val;
    optidx=[];
    optcentroids=[];
    for repindx=1:repeats,       % The program will run the function kmeans a predetermined number of times, looking for 
        trace = full(trace); % the most fitting clustering solution of the sample set.
        try
           %[idx centroids sumd] = kmeans(trace',knum, 'Distance', 'cosine', 'Display', displaymode, 'EmptyAction', 'singleton');
           [idx centroids sumd] = kmeans(trace',knum, 'Distance', clusternorm,  'EmptyAction', 'singleton','Replicates',1,'Display',display);
        catch %Error handling for the "Zero centroid" problem
            [errmsg, msg_id] = lasterr;
            switch (msg_id)
                case 'stats:kmeans:ZeroCentroid'
                    error_output= ['Zero centroid created when trying for ' num2str(knum) ' clusters.'];
                    %disp (error_output);
                    sumd = val; %Giving sumd an infinite value makes sure this iteration won't be chosen as the best fit
                otherwise
                    %error(msg_id,errmsg);
                    sumd = val;
            end    
        end

        %do some extra work if stats are requested
        if nargout==4;   stats(repindx,1)=mean(sumd);end
        if nargout==4 | usesil; sumd2 = mean(silhouette(trace',idx,'cosine'));end
        if nargout==4 & usesil;   stats(repindx,2)=sumd2;end
        if usesil; sumd=sumd2; end

        curd = sum(sumd);
        %calculate silhouette for each replica or use minimum distance for best replica?
        if usesil
            if curd > optd,               % If the current iteration had a higher sillhouette than all the previous ones,
                optd = curd;              % its parameters are kept as the best current fit.
                optcentroids = centroids;
                optidx = idx;
            end
        else
            if curd < optd,               % If the current iteration had a lower sumd than all the previous ones,
                optd = curd;              % its parameters are kept as the best current fit.
                optcentroids = centroids;
                optidx = idx;
            end
        end
    end

%special rules when k=1
elseif knum == 1,
    optcentroids = sum(trace, 2)/size(trace,2);
    %now let us normalize...
    optcentroids=    optcentroids ./ (ones(size(optcentroids,1),1) * sqrt(sum(optcentroids.^2)));
    %optcentroids=optcentroids/norm(optcentroids);
    %ie optcentroids=optcentroids/norm(optcentroids);
    optidx = ones(size(trace,2),1);
else
    error('Invalid number of clusters');
end

%calculate the sillhouette for the chosen replica
if nargout >= 3,
    silhindx = 1:size(trace,2);    
    if knum == 1,
        
		if strcmp(clusternorm, 'sqEuclidean')
			meansilh = tgs_carpet_metric(trace(:,silhindx), 'none'); %omer with aron 05/07/2010
		else
			meansilh = tgs_carpet(trace(:,silhindx),'none');
		end;
		
        meansilh = mean(meansilh(:));
% knum
% meansilh
%         if nargout==4; 
% 			
% 			if strcmp(clusternorm, 'sqEuclidean')
%
% 			else
% 				sumd=mean(tgs_H(trace,optcentroids));
% 			end;
% 			
%             if usesil;
% 				stats=[sumd*ones(replicas,1) meansilh*ones(replicas,1)];
% 			else
% 				stats=[sumd*ones(replicas,1) meansilh*zeros(replicas,1)];
% 			end
%         end
	else
        if usesil;
			meansilh=optd;
        else
            traceindex=~(sum(trace,2)==0); %remove all 0 lines as this crashes silhouette!
            if ~isempty(optidx) & ~isempty(silhindx);
                meansilh = mean(silhouette(trace(traceindex,silhindx)',optidx(silhindx), clusternorm)); %omer
% knum
% meansilh
            else
                meansilh=-1;
            end
        end
    end
end
warning on MATLAB:divideByZero
warning on stats:kmeans:EmptyCluster
warning on stats:kmeans:FailedToConverge
