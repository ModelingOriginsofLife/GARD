function [carpet_mat,trace,trace2] = tgs_carpet(varargin)
%
%function carpet_mat = tgs_carpet(trace,display,analysis_mode, color); 
%function carpet_mat = tgs_carpet(trace,tags,display,analysis_mode, color); 
%function carpet_mat = tgs_carpet(trace1, trace2,display,analysis_mode, color); 
%function carpet_mat = tgs_carpet(comps1,tags1,comps2,tags2,display,analysis_mode, color); 
%creates a compositional carpets, PNAS style.
%inputs:
%trace - the trace matrix or matrices. if 2 are given then they are compared with each other
%tags - if given then the trace is first sorted by tags
%comps - two sets of compotypes and generation tags input will result in the compotype carpet
%display - graphic output format. 'default'or 'd' - carpet with colorbar, 'naked' - plain carpet, 'none' or 'n' - no carpet, just computation.
%analysis_mode - 'all' - analyses all the trace, 
%   'divisionx' - analyses each gen/xth generation. - DEFAULT 
%   'averagex' - each gen/x generations is averaged into a single bin
%   'smoothx' - a moving average of window size x
%color - color option. 'reduced_jet' use a color scale w/o the edge dark red and dark,  blue, called tgs_reduced_jet - DEFAULT
%'jet' - default matlab color scheme,  'bw' - black and white
%outputs:
%carpet_mat - the H value comparing every generation (column) in the trace with each other
%trace,trace2 - the traces after the analysis mode has been performed

%   20050731 v0.001 based on Daniel Segre code and function from 2002
%   modified version by Aia Oz feb 2006 limits the sizes of carpets to
%   1000*1000 by using every NG/1000 vector
%   modified december 2006 to narrow down the color scale.
%  10/2007 - imagesc scale limits are now [0 1]

maxgens=8000;
trace2=[];
%check which form the input is
if  nargin>=4 & ~ischar(varargin{4})
    %tgs_carpet(comps1,tags1,comps2,tags2,display,analysis_mode, color); 
    ostart=5;   
    similarity=0.7;
    [trace,trace2]=carpet_compotypes(varargin{1},varargin{2},varargin{3},varargin{4},similarity);
elseif nargin >=2 & ~ischar(varargin{2}) & min(size(varargin{2}))==1
    %tgs_carpet(trace,tags,display,analysis_mode, color); 
    ostart=3;
    trace=varargin{1};
    trace=tgs_sorttrace(varargin{1},varargin{2});
elseif nargin >=2 & ~ischar(varargin{2})
    %tgs_carpet(trace1, trace2,display,analysis_mode, color); 
    ostart=3;
    trace=varargin{1};
    trace2=varargin{2};
else
    %tgs_carpet(trace,display,analysis_mode, color); 
    ostart=2;
    trace=varargin{1};
end

display = 'default' ;
analysis_mode = 'division';
color = 'reduced_jet';
if nargin>=ostart; display=varargin{ostart}; end
if nargin>=ostart+1; analysis_mode=varargin{ostart+1}; end
if nargin>=ostart+2; color=varargin{ostart+1}; end

[NG,gen]=size(trace);
gen2=0;if ~isempty(trace2);[NG2,gen2]=size(trace2); end

%%%%%%%%%%%%%%%%%%
if strncmp(analysis_mode,'division',8)
    if length(analysis_mode)>8; maxgens=str2num(analysis_mode(9:end)); end
    if gen>maxgens
        W=int16(gen/maxgens);
        sample_mat1=zeros(NG,maxgens);
        m=1;
        for j=1:W:gen
            sample_mat1(:,m)=trace(:,j);
            m=m+1;
        end
        trace=sample_mat1;
    end
    if gen2>maxgens
        W=int16(gen2/maxgens);
        sample_mat2=zeros(NG2,maxgens);
        m=1;
        for j=1:W:gen2
            sample_mat2(:,m)=trace2(:,j);
            m=m+1;
        end
        trace2=sample_mat2;
    end
end

if strncmp(analysis_mode,'average',7)
    if length(analysis_mode)>7; maxgens=str2num(analysis_mode(8:end)); end
    if gen>maxgens
        W=int16(gen/maxgens);
        sample_mat1=zeros(NG,maxgens);
        m=1; j1=1;
        for j=W:W:gen
            sample_mat1(:,m)=mean(trace(:,j1:j),2);
            m=m+1; j1=j+1;
        end
        trace=sample_mat1;
    end
    if gen2>maxgens
        W=int16(gen2/maxgens);
        sample_mat2=zeros(NG2,maxgens);
        m=1; j1=1;
        for j=W:W:gen2
            sample_mat2(:,m)=mean(trace2(:,j1:j),2);
            m=m+1;j1=j+1;
        end
        trace2=sample_mat2;
    end
end

if strncmp(analysis_mode,'smooth',6)
    if length(analysis_mode)>6; maxgens=str2num(analysis_mode(7:end)); else; maxgens=2; end
    sample_mat1=zeros(size(trace));
    for x=1:size(trace,1);
        sample_mat1(x,:) = filter((ones(1,maxgens)/maxgens),1,trace(x,:));
    end
    trace=sample_mat1;
    sample_mat2=zeros(size(trace2));
    for x=1:size(trace2,1);
        sample_mat2(x,:) = filter((ones(1,maxgens)/maxgens),1,trace2(x,:));
    end
    trace2=sample_mat2;
end

%normalize
epsi=0.00000001;
trace = trace ./ (ones(NG,1) * sqrt(sum(trace.^2))+epsi);
if ~isempty(trace2); trace2 = trace2 ./ (ones(NG2,1) * sqrt(sum(trace2.^2))+epsi); end

if ~isempty(trace2)
    carpet_mat = trace' * trace2 ;
else
    carpet_mat = trace' * trace ;
end

%%%%%%%%%%%%%%%%%%%%%

switch display,
    case {'default','d'}
        imagesc(carpet_mat,[0 1]);
        set(gca,'YDir','normal');
        colorbar;

    case 'naked'
        imagesc(carpet_mat, [0 1]);
        set(gca,'YDir','normal');

    case {'none','n'}
       color='jet' ;
    otherwise
        error('Invalid display method');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch color,
    case 'reduced_jet'
        tgs_colors('jet');
    case 'jet'
        ;
    case 'bw'
        set(gcf, 'Colormap', gray);
end

        return;

  %-------------------------------------------------------------------------------------------------------------------------------------------      
function [trace,trace2]=carpet_compotypes(comps1,tags1,comps2,tags2,similarity)
if nargin<5; similarity=0.7; end
%sort first set by size
[comps1,tags1,count1]=tgs_aligncomps(comps1,tags1);
%align second set by similairty with first set
if ~isempty(comps2);
[comps1,tags1,count1,comps2,tags2,count2,angles]=tgs_aligncomps(comps1,tags1,comps2,tags2,similarity,0);
end
%handle the sizes of the compotypes in each set
set1=[];
for x=1:size(comps1,2)
    set1=[set1 repmat(comps1(:,x),1,count1(x))];
end
set2=[];
if ~isempty(comps2)
for x=1:size(comps2,2)
    set2=[set2 repmat(comps2(:,x),1,count2(x))];
end
end
trace=set1; trace2=set2; 