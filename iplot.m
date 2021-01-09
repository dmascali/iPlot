function iplot(varargin)
global IPLOT_FIG_POS 

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
s = size(varargin{1});

%Find out input variable names:
labels = cell(1,nVariable);
title_str = [];
for l = 1:nVariable
    labels{l} = inputname(l);
    if isempty(labels{l})
        labels{l} = ['Input',num2str(l)];
    end
    title_str = [title_str,', ',labels{l}];
end
title_str(1:2) = [];

figure_title = sprintf(['iPlot - variable(s): ',title_str]);
if isempty(IPLOT_FIG_POS)
    h = figure('Name',figure_title,'NumberTitle','off','MenuBar', 'None');
else
    h = figure('Name',figure_title,'NumberTitle','off','MenuBar', 'None','Position',IPLOT_FIG_POS);
end

%initialize config structure
cfg.showing = '';
cfg.showingLegend = {false,true};
cfg.type = 'raw'; %raw or fft
cfg.mode = {'seq','std+','std-','random'};%,'stdA','stdD'}; %seq, random, std
cfg.columns = 1:1:s(2);
cfg.indx = 0;
cfg.indx_max = s(2); %number of columns
cfg.row_number = s(1); %number of rows
cfg.nVariable = nVariable;
cfg.labels = labels;
cfg.lnwidths = [0.5 0.6 0.7 0.8 0.9 1 1.2 1.4 1.6 1.8 2 2.5 3 4];
cfg.lnwidthsIndx = find(cfg.lnwidths==1);
cfg.ylim_mode = 'auto'; %auto or lock
cfg.ylim = [];
cfg.fft = [];
cfg.Fs = 1;
%----------------------------

% print help screen
cfg = help_screen(cfg);

%wait for key press
guidata(h,cfg);
set(h,'KeyPressFcn',{@switcher,varargin{:}});

return
end

function switcher(src,event,varargin)
global IPLOT_FIG_POS
%retrive config structure
cfg = guidata(src);

switch event.Key
    %--------------------------navigate------------------------------------
    case {'a'} %backward
        cfg = update_cfg_a(cfg);
        plot_column(cfg,varargin{:});        
    case {'d'} %forward
        cfg = update_cfg_d(cfg);
        plot_column(cfg,varargin{:});
    %----------------------------raw/fft-----------------------------------
    case {'f'} %switch bettwen raw data plotting and spectrum
        switch cfg.type
            case {'raw'}
                cfg.type = 'fft';
                %if it's the first call to fft, let's create the spectrum
                if isempty(cfg.fft)
                   cfg = calculate_spectrum(cfg,varargin{:}); 
                end
                plot_column(cfg,varargin{:});
            case {'fft'}
                cfg.type = 'raw';
                plot_column(cfg,varargin{:});
        end
    %----------------------------modality----------------------------------
    case {'r'} %switch bettwen seq, random, std...
        cfg.mode = circshift(cfg.mode,-1);
        switch cfg.mode{1}
            case {'seq'}
                cfg.columns = 1:1:cfg.indx_max;
                %resest index
                cfg.indx = 0;
                print_title(cfg)
            case {'random'}
                cfg.columns = randperm(cfg.indx_max);
                cfg.indx = 0;
                print_title(cfg)
            case {'std+'}
                stds = std(varargin{1});
                [~,cfg.columns] = sort(stds,'descend');
                cfg.indx = 0;
                print_title(cfg)
            case {'std-'}
                stds = std(varargin{1});
                [~,cfg.columns] = sort(stds,'ascend');
                cfg.indx = 0;
                print_title(cfg)
        end
    %--------------------------legend--------------------------------------    
    case {'l'}
        cfg.showingLegend = circshift(cfg.showingLegend,-1);
        switch cfg.showingLegend{1}
            case {true} 
                legend(cfg.labels,'location','best');
            case {false}
                legend off;
        end
    %----------------------lock ylim---------------------------------------    
    case {'e'}
        switch cfg.ylim_mode
            case {'auto'} 
               cfg.ylim_mode = 'lock';
               cfg.ylim = get(gca,'ylim');
               print_title(cfg);
            case {'lock'}
               cfg.ylim_mode = 'auto';
               print_title(cfg);
        end
    %----------------------linewidth---------------------------------------
    case {'equal'}
        cfg = update_lnwd_plus(cfg);
        %update current
        lines = findobj(gcf,'Type','Line');
        for l = 1:numel(lines); lines(l).LineWidth = cfg.lnwidths(cfg.lnwidthsIndx); end
    case {'hyphen'}
        cfg = update_lnwd_minus(cfg);   
        lines = findobj(gcf,'Type','Line');
        for l = 1:numel(lines); lines(l).LineWidth = cfg.lnwidths(cfg.lnwidthsIndx); end
    %----------------------set param---------------------------------------
    case {'s'}
        prompt = {'Enter sampling frequency (Hz):'};
        dlgtitle = 'Set Parameters';
        dims = [1 35];
        definput = {'1'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        cfg.Fs = str2num(answer{1});
        % update freq:
        n = 2^nextpow2(cfg.row_number);
        cfg.freq = 0:(cfg.Fs/n):(cfg.Fs/2);
        % update plot:
        plot_column(cfg,varargin{:});   
    %--------------------------help----------------------------------------
    case {'h'}
        switch cfg.showing
            case {'help'}
                %if it was already on help, go back to what previously
                %plotted
                cfg.showing = 'data';
                plot_column(cfg,varargin{:});  
            otherwise
                help_screen;
                cfg.showing = 'help';
        end
    %--------------------------quit----------------------------------------    
    case {'q'}
        % before exiting save 
        IPLOT_FIG_POS = src.Position;
        close(gcf); 
        return
end

%update config structure
guidata(src,cfg);

return
end

function print_title(cfg)
if strcmp(cfg.showing,'help')
    % do nothing
    return
end
title(['Modality: \bf',cfg.type,'     \rmColumn ordering: \bf',cfg.mode{1},'     \rmYlim: \bf',cfg.ylim_mode],'fontweight','normal'); 
try % available in the 
    ax = gca;
    ax.TitleHorizontalAlignment = 'left'; ax.TitleHorizontalAlignment = 'left';
end
return
end

function cfg = update_lnwd_plus(cfg)
if cfg.lnwidthsIndx ~= length(cfg.lnwidths)
    cfg.lnwidthsIndx = cfg.lnwidthsIndx +1; 
end
cfg.showing = 'data';
return
end

function cfg = update_lnwd_minus(cfg)
if cfg.lnwidthsIndx ~= 1
    cfg.lnwidthsIndx = cfg.lnwidthsIndx -1; 
end
cfg.showing = 'data';
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
if cfg.indx > 1
    cfg.indx = cfg.indx -1; 
end
if cfg.indx == 0
    cfg.indx = 1; 
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
if cfg.indx <= 0 %take care of odd situations where indx = 0
    cfg.indx = 1;
end
switch cfg.type
    case {'raw'}
        hold on;
        for l = 1:cfg.nVariable
            plot(varargin{l}(:,cfg.columns(cfg.indx)),'linewidth',cfg.lnwidths(cfg.lnwidthsIndx));
        end
        hold off;
        ylabel(['Column \bf',num2str(cfg.columns(cfg.indx)),'/',num2str(cfg.indx_max)],'Fontweight','normal');
    case {'fft'}
        hold on;
        for l = 1:cfg.nVariable
            plot(cfg.freq,cfg.fft{l}(:,cfg.columns(cfg.indx)),'linewidth',cfg.lnwidths(cfg.lnwidthsIndx));
        end
        hold off;
        xlabel('Frequency');
        ylabel(['Spectrum, column \bf',num2str(cfg.columns(cfg.indx)),'/',num2str(cfg.indx_max)],'Fontweight','normal');
end

switch cfg.ylim_mode
    case {'lock'}
        ylim(cfg.ylim);
end
print_title(cfg);
box on;

switch cfg.showingLegend{1}
    case {true} 
        legend(cfg.labels,'location','best');
    case {false}
        legend off;
end

return
end

function cfg = calculate_spectrum(cfg,varargin)
%concantenate matrixes
y = [];
for l = 1:cfg.nVariable
    y = cat(2,y,varargin{l});
end

n = 2^nextpow2(cfg.row_number);

Y = fft(y,n);
Y = abs(Y/cfg.row_number);
Y = Y(1:n/2+1,:);
Y(2:end-1,:) = 2*Y(2:end-1,:);

cfg.freq = 0:(cfg.Fs/n):(cfg.Fs/2);

%split back data
for l = 1:cfg.nVariable
    cfg.fft{l} = Y(:,(1:cfg.indx_max));
    Y(:,1:cfg.indx_max) = [];
end

return
end

 

