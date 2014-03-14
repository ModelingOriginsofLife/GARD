function [exp1, g1]=omer_curvefit_v10(t, ct, display)
% [exp1, g1]=omer_curvefit_v10(t, ct, display)
% 08/08/2010 Omer, generated using 'cftool' & file->export M file
% display=none to turn off display, =1 to display only exp1
% exp is a data structure holding the fitted parameters.
% g1 is goodness of fit data structure.
% exp1: (1-h0)*exp(-x/tau1)+h0
% 16/06/2011 GARD10, by Omer Markovitch

if ~exist('display', 'var') | isempty(display); display=1;
elseif display<1 | display=='none'; display=0;
end;

if size(t, 2)~=1; t=t'; end;
if size(ct, 2)~=1; ct=ct'; end;

exp1=[];
g1=[];

% Set up figure to receive datasets and fits
if (display>0);
	f_ = clf;
	figure(f_);
	set(f_,'Units','Pixels','Position',[318 115 680 481]);
	legh_ = []; legt_ = {};   % handles and text for legend
	xlim_ = [Inf -Inf];       % limits of x axis
	ax_ = axes;
	set(ax_,'Units','normalized','OuterPosition',[0 0 1 1]);
	set(ax_,'Box','on');
	axes(ax_); hold on;

	% --- Plot data originally in dataset "ct vs. t"
	t = t(:);
	ct = ct(:);
	h_ = line(t,ct,'Parent',ax_,'Color',[0.333333 0 0.666667],...
		'LineStyle','none', 'LineWidth',1,...
		'Marker','.', 'MarkerSize',12);
	xlim_(1) = min(xlim_(1), min(t));
	xlim_(2) = max(xlim_(2), 500); %,max(t));
	legh_(end+1) = h_;
	legt_{end+1} = 'ct vs. t';

	% Nudge axis limits beyond data limits
	if all(isfinite(xlim_))
		xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
		set(ax_,'XLim',xlim_)
	else
		set(ax_, 'XLim',[-1, 2000]);
	end
	
end;
	
% --- Create fit "exp1"
fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[0 0.001],'Upper',[1 Inf], 'MaxFunEvals', 750, 'Tolfun', 1e-5); %omer
ok_ = isfinite(t) & isfinite(ct);
if ~all( ok_ )
	warning( 'GenerateMFile:IgnoringNansAndInfs', ...
		'Ignoring NaNs and Infs in data' );
end
st_ = [0.4 1.0];
set(fo_,'Startpoint',st_);
ft_ = fittype('(1-h0)*exp(-x/tau1)+h0',...
	'dependent',{'y'},'independent',{'x'},...
	'coefficients',{'h0', 'tau1'});
[cf_ g1]= fit(t(ok_), ct(ok_), ft_, fo_); % Fit this model using new data
exp1=cf_;

% Plot this fit
if display>0;
	h_ = plot(cf_,'fit',0.95);
	legend off; % turn off legend from plot method call
	set(h_(1),'Color',[1 0 0],...
		'LineStyle','-', 'LineWidth',2,...
		'Marker','none', 'MarkerSize',6);
	legh_(end+1) = h_(1);
	legt_{end+1} = '(1-h0)*exp(-x/tau1)+h0';

	hold off;
	leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'};
	h_ = legend(ax_,legh_,legt_,leginfo_{:}); % create legend
	set(h_,'Interpreter','none');
	xlabel(ax_,''); % remove x label
	ylabel(ax_,''); % remove y label
end;

return;
