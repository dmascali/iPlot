function iplot(varargin)
%IPLOT(X) plots the columns in X (or their spectral amplitude) interactively,
% one at a time. 
%IPLOT(X,Y,...) plots the columns from X,Y..., one above the other for easy 
% comparison.
%
%IPLOT functionalities are triggered by pressing the following keys:
%
% Navigate through columns:
%  D : plot the next column
%  A : plot the previous column
%  R : toggle between different column ordering:
%      sequential: from 1 to N (number of columns) {default}
%      std+      : column sorted by variance (of the first input), descending 
%      std-      : column sorted by variance (of the first input), ascending 
%      random    : random ordering 
%
% Plotting mode:
%  F : toogle between plotting modality:
%      raw : plot raw column {default}
%      fft : plot the spectral amplitude of the column
%
% Appearance:
%  E :  toogle between y-limit modality:
%       auto : automatically adjust limits {default}
%       lock : lock the current y-limits (you can also specify limits at
%              the matlab command window: ylim([a b]) and then locking them)
%  L :  show legend 
%  +/-: adjust linewidth
%
% Miscellaneous:
%  S : open a dialog box for specifying the sampling frequency (so that 
%      spectra reflect real frequencies)
%  Q : quit
%  H : show help
%__________________________________________________________________________
% Daniele Mascali
% ITAB, Chieti, 2021 
% danielemascali@gmail.com

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
    defpos = get(groot, 'DefaultFigurePosition');
    h = figure('Name',figure_title,'NumberTitle','off','MenuBar', 'None','Position', [defpos(1) defpos(2)-275 920 695]);
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
    case {'a','leftarrow'} %backward
        cfg = update_cfg_a(cfg);
        cfg = plot_column(cfg,varargin{:});        
    case {'d','rightarrow'} %forward
        cfg = update_cfg_d(cfg);
        cfg = plot_column(cfg,varargin{:});
    %----------------------------raw/fft-----------------------------------
    case {'f'} %switch bettwen raw data plotting and spectrum
        switch cfg.type
            case {'raw'}
                cfg.type = 'fft';
                %if it's the first call to fft, let's create the spectrum
                if isempty(cfg.fft)
                   cfg = calculate_spectrum(cfg,varargin{:}); 
                end
                cfg.showing = 'data';
                cfg = plot_column(cfg,varargin{:});
            case {'fft'}
                cfg.type = 'raw';
                cfg = plot_column(cfg,varargin{:});
        end
    %----------------------------modality----------------------------------
    case {'r'} %switch bettwen seq, random, std...
        if strcmp(cfg.showing,'help')
            return
        end
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
        if strcmp(cfg.showing,'help')
            return
        end
        cfg.showingLegend = circshift(cfg.showingLegend,-1);
        switch cfg.showingLegend{1}
            case {true} 
                legend(cfg.labels,'location','best');
            case {false}
                legend off;
        end
    %----------------------lock ylim---------------------------------------    
    case {'e'}
        if strcmp(cfg.showing,'help')
            return
        end
        switch cfg.ylim_mode
            case {'auto'}
                cfg.ylim_mode = 'lock';
                cfg.ylim = get(gca,'ylim');
            case {'lock'}
                cfg.ylim_mode = 'auto';
        end
        print_title(cfg);
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
        if ~isempty(answer)
            cfg.Fs = str2num(answer{1});
            % update freq:
            n = 2^nextpow2(cfg.row_number);
            cfg.freq = 0:(cfg.Fs/n):(cfg.Fs/2);
            % update plot:
            if strcmp(cfg.showing,'data')
                cfg = plot_column(cfg,varargin{:});   
            end
        end
    %--------------------------help----------------------------------------
    case {'h'}
        switch cfg.showing
            case {'help'}
                %if it was already on help, go back to what previously
                %plotted
                cfg.showing = 'data';
                cfg = plot_column(cfg,varargin{:});  
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
% plot functions in backgroun
x = 0:0.1:4*pi;
y = [4*sin(x);3*sin(2*x);2*sin(3*x);1*sin(4*x)]';
plot(x,y,'linewidth',1.5);
ylim([-13, 4.5]);
xlim([0, x(end)]);
%----------------------------
%title
text(0.5, 0.8,'\bfi\itPlot','FontSize',50,'HorizontalAlignment','center', 'Units', 'Normalized');
text(0.5, 0.70,'Interactive Plot','FontSize',25,'HorizontalAlignment','center', 'Units', 'Normalized');

%legend help
text(0.5, 0.47,'\bfKeyboard Map','FontSize',12,'HorizontalAlignment','center', 'Units', 'Normalized');

%left_colulm
lf_entry_x = 0.05+0.1;
lf_entry_title_x = lf_entry_x -0.015; 

rh_entry_x = 0.55+0.1;
rh_entry_title_x = rh_entry_x -0.015; 

%common
deltay = 0.05; %interline
entry_start_y = 0.4;
entry_fontsize = 12;

%left column
text(lf_entry_title_x, entry_start_y-0.*deltay,'\bf\itNavigate through columns','FontSize',entry_fontsize,'HorizontalAlignment','left', 'Units', 'Normalized');
text(lf_entry_x, entry_start_y-1.*deltay,'D - move forward','FontSize',entry_fontsize,'HorizontalAlignment','left', 'Units', 'Normalized');
text(lf_entry_x, entry_start_y-2.*deltay,'A - move backward','FontSize',entry_fontsize,'HorizontalAlignment','left', 'Units', 'Normalized');
text(lf_entry_x, entry_start_y-3.*deltay,'R - column ordering:','FontSize',entry_fontsize,'HorizontalAlignment','left', 'Units', 'Normalized');
text(lf_entry_x, entry_start_y-4.*deltay,'     [sequential/std+/std-/random]','FontSize',entry_fontsize,'HorizontalAlignment','left', 'Units', 'Normalized');
text(lf_entry_title_x, entry_start_y-5.*deltay,'\bf\itPlotting mode','FontSize',entry_fontsize,'HorizontalAlignment','left', 'Units', 'Normalized');
text(lf_entry_x, entry_start_y-6.*deltay,'F - raw or spectral amplitude','FontSize',entry_fontsize,'HorizontalAlignment','left', 'Units', 'Normalized');
text(lf_entry_x, entry_start_y-7.*deltay,'     [raw/fft]','FontSize',entry_fontsize,'HorizontalAlignment','left', 'Units', 'Normalized');

%right column
text(rh_entry_title_x, entry_start_y-0.*deltay,'\bf\itAppearance','FontSize',entry_fontsize,'HorizontalAlignment','left', 'Units', 'Normalized');
text(rh_entry_x, entry_start_y-1.*deltay,'E   - ylim [auto/lock]','FontSize',entry_fontsize,'HorizontalAlignment','left', 'Units', 'Normalized');
text(rh_entry_x, entry_start_y-2.*deltay,'L   - legend [on/off]','FontSize',entry_fontsize,'HorizontalAlignment','left', 'Units', 'Normalized');
text(rh_entry_x, entry_start_y-3.*deltay,'+/- - adjust linewidth','FontSize',entry_fontsize,'HorizontalAlignment','left', 'Units', 'Normalized');
text(rh_entry_title_x, entry_start_y-4.*deltay,'\bf\itMiscellaneous','FontSize',entry_fontsize,'HorizontalAlignment','left', 'Units', 'Normalized');
text(rh_entry_x, entry_start_y-5.*deltay,'S - set frequency for fft','FontSize',entry_fontsize,'HorizontalAlignment','left', 'Units', 'Normalized');
text(rh_entry_x, entry_start_y-6.*deltay,'Q - quit','FontSize',entry_fontsize,'HorizontalAlignment','left', 'Units', 'Normalized');
text(rh_entry_x, entry_start_y-7.*deltay,'H - show this help','FontSize',entry_fontsize,'HorizontalAlignment','left', 'Units', 'Normalized');


set(gca,'Xtick',[],'Ytick',[], 'XtickLabel',{},'YtickLabel',{},'box','on');
return
end

function cfg = plot_column(cfg,varargin)
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
        xlabel('row');
        xlim([0 cfg.row_number]);
    case {'fft'}
        hold on;
        for l = 1:cfg.nVariable
            plot(cfg.freq,cfg.fft{l}(:,cfg.columns(cfg.indx)),'linewidth',cfg.lnwidths(cfg.lnwidthsIndx));
        end
        hold off;
        xlabel(['Frequency (Hz) [sampling freq=',num2str(cfg.Fs),'Hz]']);
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

 

