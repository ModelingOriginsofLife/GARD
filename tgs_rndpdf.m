function [retval probsum dt]=tgs_rndpdf(probvec); %omer 26/07/12

%function [retval probsum dt]=tgs_rndpdf(probvec);
%genrate a random integer value given discrete weights
%inputs:
%probvec - a vector different numbers serving as weights 
%k - the number of elements to return (default=1)
%outputs:
%retval - an index on the probability vector
%probsum - sum of input vector
%see bottom of function code for more explanations

%if nargin>1; [retval probsum]=tgs_rndpdfk(probvec,k); return; end

if size(probvec,2)==1,  % if the vector is not oriented properly, "vertical"
    probvec = probvec'; %    rotate the vector to a proper alignment
end
probvec2 = [0 cumsum(probvec)]; %omer 26/07/12  % stores the sum of the values until now
                                % for each index
probsum = probvec2(end);         % the sum of probabilities 
                                % (merhav hapitaron) is the sum of all the 
                                % values in the vector- the end of the sums
                                % vector

randval = rand;         % get a random value from 0 to 1
while (randval == 0),   % this ensures that there is always a choice from 
    randval = rand;     % one of the vector cells - no zero value
end

randval = randval * probsum;    % multiplying the random number by the sum
                                % gives us a usable random choice

indx = find (probvec2 < randval);    % indx will be a vector holding the 
                                    % indexes of all the values in the sums
                                    % vector that are smaller than the 
                                    % randomized fraction.

retval = indx(end);             % return the last of the values in the 
                                % index vector. This means that the index
                                % of the largest value in the edge is
                                % returned. 
dt=probvec(retval); %omer 26/07/12
return


function [retval probsum]=tgs_rndpdfk(probvec,k)
if size(probvec,2)==1; probvec = probvec'; end
probvec = [cumsum(probvec)];                                  
probsum = probvec(end);         
retval=zeros(length(probvec),1);
kr=randsample(probsum(end),k,false);
for x=1:k; add=find(probvec>=kr(x));add=add(end);retval(add)=retval(add)+1; end

% 	DESCRIPTION
%		This function receives a vector of probabilities as input. It then
%		stores the sums of those probabilities in the vector, which means
%		that for the expected positive values, the cell values increase
%		with the index. 
%       As an example consider the following vector: 
%
%       example = [ 1 0 0 0 2 5 3 1 ] 
%
%       the vector of accumulated sums will be: 
%       [0     1     1     1     1     3     8    11    12]
%
%       the last value for this vector, called probsum, is 12.
%
%       We'll take a nonzero randomly generated value, 0.2311.
%       multiplied by the probsum(12), this gives 2.7732.
%       we now build a vector of the indexes of all the values smaller than
%       this number:
%       
%       [1     2     3     4     5]
%
%       the function returns the last of these values - the index 5, which
%       points to the number 2 in our original vector. 
%     
%       the larger weights will be picked more often this way - notice the
%       big gap between the values 3 and 8 in the second vector, which
%       means for the range of 0.25 to 0.6666 out of the possibilities we
%       get the number 5, which is a range of 0.4160. the probability of 5
%       out of 12 is 0.4167, which is pretty close. So each cell is
%       represented according to its weight