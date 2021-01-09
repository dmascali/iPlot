function iplot(varargin)

if nargin > 1
    %check dimension consistency among matrices
    for l = 2:nargin
        if ~isequal(size(varargin{l-1}),size(varargin{l}))
            in1 = inputname(l-1); in2 = inputname(l);
            if ~isempty(in1) && ~isempty(in2)
                error(['Matrix ''',in1,''' has a different size from ''',in2,'''.']); 
            else
                error('Input matrices have different dimensions.');
            end
        end
    end 
end
% Number of input matrices
nVariable = nargin;

%Find out input variable names:
labels = cell(nVariable,1);
for l = 1:nVariable
    labels{l} = inputname(l);
    if isempty(labels{l})
        labels{l} = ['Input',num2str(l)];
    end
end

h = figure;
%initialize config structure
cfg.showing = '';
cfg.showingLegend = false;
cfg.type = 'raw'; %raw or fft
cfg.mode = 'seq'; %seq, random, std
cfg.indx = 0;
cfg.indx_max = size(varargin{1},2);
cfg.nVariable = nVariable;
cfg.labels = labels;
%----------------------------

% print help screen
cfg = help_screen(cfg);

%wait for key press
guidata(h,cfg);
set(h,'KeyPressFcn',{@switcher,varargin{:}});


return
end

function switcher(src,event,varargin)

%retrive config structure
cfg = guidata(src);

switch event.Key
  
    case {'a'}
        cfg = update_cfg_a(cfg);
        plot_column(cfg,varargin{:});        
    case {'d'}
        cfg = update_cfg_d(cfg);
        plot_column(cfg,varargin{:});
    case {'f'} %switch bettwen raw data plotting and spectrum
        switch cfg.type
            case {'raw'}
                cfg.type = 'fft';
                plot_column(ts(:,cfg.indx),cfg);
            case {'fft'}
                cfg.type = 'raw';
                plot_column(ts(:,cfg.indx),cfg);
        end
        
    case {'l'}
        switch cfg.showingLegend
            case {true} %then toggle off
                legend off;
                cfg.showingLegend = false;
            case {false}
                legend(cfg.labels,'location','best');
                cfg.showingLegend = true;
        end
    case {'h'}
        switch cfg.showing
            case {'help'}
                %if it was already on help, go back to what previously
                %plotted
                plot_column(cfg,varargin{:});  
                cfg.showing = 'data';
            otherwise
                help_screen;
                cfg.showing = 'help';
        end
        
    case {'q'}
        close(gcf); 
        return
end

%update config structure
guidata(src,cfg);

return
end

function cfg = update_cfg_d(cfg)
if cfg.indx ~= cfg.indx_max
    cfg.indx = cfg.indx +1; 
end
cfg.showing = 'data';
return
end

function cfg = update_cfg_a(cfg)
if cfg.indx ~= 1
    cfg.indx = cfg.indx -1; 
end
cfg.showing = 'data';
return
end


function cfg = help_screen(cfg)
cfg.showing = 'help';
clf;
text(0.5, 0.5,'HELP');
set(gca,'Xtick',[],'Ytick',[], 'XtickLabel',{},'YtickLabel',{},'box','on');
return
end

function plot_column(cfg,varargin)
clf;
switch cfg.type
    case {'raw'}
        hold on;
        for l = 1:cfg.nVariable
            plot(varargin{l}(:,cfg.indx));
        end
        hold off;
    case {'fft'}
        %TODO
end
ylabel(['Column ',num2str(cfg.indx)],'Fontweight','bold');

return
end


 

