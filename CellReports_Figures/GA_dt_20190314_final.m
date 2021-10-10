% Written by Keo 20171008 
% Modified by Keo 20190314, 20201226, 20211011
% ============================ Contents ===================================
% Figure 1. PSTH of each modality
%           1.1 Normalized or Un-normalized PSTH  
%           1.2 ΔPSTH
% Figure 2. Choice divergence 
%           2.1 Choice divergence 
%           2.2 Distribution of divergence time  
% Figure 3. PSTH of each heading angle
%           2.1 PSTH sorting by angles
%           2.2 Δ PSTH sorting by angles
% Figure 4. Partial correlation
%           4.1 Partial correlation cell by cell
%           4.2 Partial correlation across time
% Figure 5. ★Fisher information (neurometric)
%           5.1 Fisher information (neurometric)
%           5.2 Neurometric function
%           5.3 3D distribution of Fisher information
% Figure 6. SVM decoder
%           6.1 Correct rate across time
%           6.2 correct rate of each condition
%           6.3 Weight
%           6.4 correct rate of each heading angles
%           6.5 Neurometric function
%           6.6 Weight vs. FI
% =========================================================================

clc;clear;
%% === Preparation ===
FigureName = 'Figure S6'; 
switch FigureName
    case 'Figure 4B-C'
        % (1) PSTH (binSize = 10; correct trials only), Monkey Messi & Fara, FEF
        FolderPath = 'Z:\Data\Tempo\Batch\BatchFile_20191010_Messi_and_Fara_FEF_dt_PSTH_bin10_CorrectTrial\';
    case 'Figure 5A-C'
        % (2) FI(neurometric) (binSize = 250; add wrong trials), Monkey Messi & Fara, FEF
        FolderPath = 'Z:\Data\Tempo\Batch\BatchFile_20190424_Messi_and_Fara_FEF_dt_PSTH_FI_bin250_AddWrongTrial\';
    case 'Figure S4'
        % (5) PSTH (binSize = 10; correct trials only), Monkey Fara, FEF, add high coherence
        FolderPath = 'Z:\Data\Tempo\Batch\BatchFile_20210205_Fara_FEF_dt_AddHighCoherence_PSTH_bin10_CorrectTrial\';
    case 'Figure S5'
        % (3) FI(neurometric) (binSize = 250; add wrong trials), Monkey Messi & Fara, FEF
        FolderPath = 'Z:\Data\Tempo\Batch\BatchFile_20190424_Messi_and_Fara_FEF_dt_PSTH_FI_bin250_AddWrongTrial\';
        angle_range = [0 15]; % absolute angle between target and horizontal line, within [0 90] 
    case 'Figure S6'
        % (4) FI(neurometric) (binSize = 250; add wrong trials), Monkey Fara, LIP
        FolderPath = 'Z:\Data\Tempo\Batch\BatchFile_20200614_Fara_LIP_dt_PSTH_bin250_AddWrongTrial\';
        angle_range = [45 60]; % absolute angle between target and horizontal line, within [0 90] 
    case 'Figure S7A'
        % (5) FI(neurometric) (binSize = 250; add wrong trials), Monkey Messi, FEF
        FolderPath = 'Z:\Data\Tempo\Batch\BatchFile_20180704_Messi_FEF_dt_PSTH_FI_bin250_AddWrongTrial\';
    case 'Figure S7B'
        % (5) FI(neurometric) (binSize = 250; add wrong trials), Monkey Fara, FEF
        FolderPath = 'Z:\Data\Tempo\Batch\BatchFile_20190424_Fara_FEF_dt_PSTH_FI_bin250_AddWrongTrial\';
    case 'Figure S7C'
        % (5) FI(neurometric) (binSize = 250; add wrong trials), Monkey Fara, LIP
        FolderPath = 'Z:\Data\Tempo\Batch\BatchFile_20200614_Fara_LIP_dt_PSTH_bin250_AddWrongTrial\';
end


cond_num = 5;% 3 - no dt;5 - include dt
t_motion_off = 2000;%!!!% 1500 - without dt; 2000 - with dt500
vel_t_offset = (t_motion_off-1500)/2;

% These parameters should be kept the same as previous processing of single file in TEMPO_GUI.m !!! 
% ===========================================================================
align_markers = {
    % Before(ms)  After(ms)    align to?
    -500,         2700,        'Stim On' ;
    -2500,        700,         'Saccade On';
    };
binSize = 10;% caution！！！
stepSize = 10;
for j = 1:size(align_markers,1) %1-StimOn, 2-SaccadeOn
    time_stamp{j} = align_markers{j,1} + binSize/2 : stepSize : align_markers{j,2} - binSize/2 ;% Time stamps
end
% ===========================================================================

transparent = 1;%whether make shade transparent
p_critical = 0.05;%note: in single-file analysis, it's 0.01.
p_value_type = 2;%1-p value bar;2-shaded error bar
isnorm = 1; %For PSTH: 0 - without normalization; 1-with normalization
issquared = 0;

LEFT = 1;
RIGHT = 2;

% 0.1 load data in the fold
Files = dir(fullfile(FolderPath,'*.mat'));
File_work_count = 0;
list_1 = load([FolderPath Files(1).name]);
for ct = 1:length(Files)
    if ct==1 & isfield(list_1,'errorFiles')
        continue;
    else
        File_work_count = File_work_count+1;
        File_normal_name{File_work_count} = Files(ct).name;
        list{File_work_count} = load([FolderPath Files(ct).name]);
    end
end
File_normal_total = File_work_count;
display(['Total = ',num2str(File_normal_total)])

% % files to analyze
% file_to_analyze = 1:File_normal_total;

% 20201226 to analyze horizontal-like targets 
if strcmp(FigureName, 'Figure S5')
    % Target angles of Monkey Messi, FEF
    angle_Messi = [45	0	-55	-55	-75	35	-15	20	-10	75	-75	45	-30	-5	-5	55	40	-60	-25	-40	5	45	15	-40	55	-40	-40	4	-25	35	-35	-40	30	70	-40	20	0	-65	30	-75	-75	-45	-80	25	-35	20	-80	65	65	40	-75	-75	-60	-30	-35	55	-75	-25	-10	10	-10	-50	65	65	0	0	-10	25	-75	25	55	55	30	30];
    % Target angles of Monkey Fara, FEF
    angle_Fara = [-80	45	25	-10	15	-40	25	0	-70	70	-50	70	0	-70	-5	65	-40	30	65	65	-30	-20	65	-60	-35	20	45	10	-10	70	0	0	0	45	45	-20	45	-20	55	60	-60	-60	60	10	35	-30	-65	50	65	-60	-75	45];
    % Target angles of Monkey Messi & Fara, FEF
    angle_all = [angle_Messi angle_Fara ];
    
    % Distribution of targets' positions 
    set(figure(1010),'color','white','position',[100,100,500,550])
    centers = -90:5:90;
    [nelements,centers] = hist(angle_all,centers)
    h = bar(centers,nelements)
    malimalihong('Angle of target(°)','Number of neurons')
    xlim([-90,90])
    set(gca,'XTick',[-80, -45,0,45,80]); 
    x_tick_label = {'-80' '-45' '0','45' '80'}
    set(gca,'XTickLabel',x_tick_label)
    set(h,'facecolor','k','edgecolor','w')
    
    if angle_range(1) == 0
        file_to_analyze = find( abs(angle_all) >= angle_range(1) & abs(angle_all) <= angle_range(2) );
    else
        file_to_analyze = find( abs(angle_all) > angle_range(1) & abs(angle_all) <= angle_range(2) );
    end
elseif strcmp(FigureName, 'Figure S6')
    % Target angles of Monkey Fara, LIP
    angle_Fara_LIP = [60	35	20	-50	30	-25	-50	-60	35	60	0	-50	60	60	40	-45	30	-10	40	30	10	30	55	25];
    
     % Distribution of targets' positions 
    set(figure(1010),'color','white','position',[100,100,500,600])
    centers = -90:5:90;
    [nelements,centers] = hist(angle_Fara_LIP,centers)
    h = bar(centers,nelements)
    malimalihong('Angle of target(°)','Number of neurons')
    xlim([-90,90])
    set(gca,'XTick',[-80, -45,0,45,80]); 
    x_tick_label = {'-80' '-45' '0','45' '80'}
    set(gca,'XTickLabel',x_tick_label)
    set(h,'facecolor','k','edgecolor','w')
    
    if angle_range(1) == 0
        file_to_analyze = find( abs(angle_Fara_LIP) >= angle_range(1) & abs(angle_Fara_LIP) <= angle_range(2) );
    else
        file_to_analyze = find( abs(angle_Fara_LIP) > angle_range(1) & abs(angle_Fara_LIP) <= angle_range(2) );
    end
else
    file_to_analyze = 1:File_work_count;
end

% file_to_analyze = [98];
display(['To analyze = ',num2str(length(file_to_analyze))])

if length(file_to_analyze) >= 2
    isGroup = 1;
elseif length(file_to_analyze) == 1
    isGroup = 0;
end
        
% 0.2 colors, markers 
symbol_color{1} = [0 0 1]; % blue - vestibular only
symbol_color{2} = [1 0 0]; % red - visual only
symbol_color{3} = [0 1 0];% green - combined in 'dt=0'
symbol_color{4} = [0 1 1];%[51 204 204]/255;% cyan
symbol_color{5} = [1 0 1];%[255 51 204]/255;% magenta
symbol_color{6} = [1 1 0];%[255 134 40]/255;% orange

symbol{1} = '-';
symbol{2} = '--';

line_width_axis = 1;% axis
line_width_curve = 1;% curve
line_width_vel = 0.8;% velocity profile

font_size = 15;

legend_str_dt = {' Ves ' '  Vis ' '     0 ' '-500' '-250'};
% legend_str_PN = {'Ves,pref' 'Ves,null' 'Vis,pref' 'Vis,null' 'Δt=0,pref' 'Δt=0,null' 'Δt=-500,pref' 'Δt=-500,null' 'Δt=-250,pref' 'Δt=-250,null'};
legend_str_PN = {'Pref —' 'Null ---'};

% 0.3 time range
t_motion_on = 0;
t_saccade_on = 0;
t_range = [t_motion_on - 400, t_motion_off + 200; ...
           t_saccade_on - 200, t_saccade_on + 200];
for j = 1:2 % to limit the maximum range
    if t_range(j,1) < align_markers{j,1} + binSize/2    t_range(j,1) = align_markers{j,1} + binSize/2; end
    if t_range(j,2) > align_markers{j,2} + binSize/2    t_range(j,2) = align_markers{j,2} + binSize/2; end 
end
% time range to be shown in figure
x_range = {t_range(1,:); t_range(2,:)};

% find motion period and its index in different 'dt' conditions
t_dt = [0 250 500];
for j = 1:size(align_markers,1) %1-StimOn, 2-SaccadeOn
    for d = 1:length(t_dt)
        if j == 1
            t_winBeg(j,d) = t_range(j,1) + t_dt(d)/2 ;
            t_winEnd(j,d) = t_range(j,2) + t_dt(d)/2 ;
        elseif j == 2
            t_winBeg(j,d) = t_range(j,1);
            t_winEnd(j,d) = t_range(j,2);
        end
        [~, idx_winBeg(j,d)] = min(abs(time_stamp{j} - t_winBeg(j,d)));
        [~, idx_winEnd(j,d)] = min(abs(time_stamp{j} - t_winEnd(j,d)));
    end
    ts{j} = time_stamp{j}(idx_winBeg(j,1):idx_winEnd(j,1));
end

%% Figure 1. PSTH of each modality
%% Figure 1.1 Normalized or Un-normalized PSTH  
if strcmp(FigureName, 'Figure 4B-C') || strcmp(FigureName, 'Figure S4')
    figure(1);
    ScrSize = get(0, 'screensize');
    set(1, 'position', [ScrSize(3)*0.05, ScrSize(4)*0.05, ScrSize(3)*0.5, ScrSize(4)*0.85],'color','white','name','PSTH of each modality');clf;

    % for normalization
    sortInd = 3;% figure sorted by stim type and heading angles
    for f = 1:length(file_to_analyze)
        ys_max(f) = 0;
        for j = 1 : 2 % 1 = Stim On, 2 = Saccade On
            for cond = 1:cond_num
                temp = max(list{1,file_to_analyze(f)}.result.PSTH{j,sortInd,cond}.ys(:));
                ys_max(f) = max(ys_max(f),temp); % for normalization
            end
        end
    end


    for j = 1 : 2 % 1 = Stim On, 2 = Saccade On
        if j == 1
            ax(1,1,j) = axes('position',[0.10 0.60 0.30 0.35]);hold on;
        else
            ax(1,1,j) = axes('position',[0.42 0.60 0.05 0.35]);hold on;
        end

        for cond = 1:cond_num
            % pre-allocate
            ys{j,cond,1} = [];         ys{j,cond,2} = [];
            ys_norm{j,cond,1} = [];    ys_norm{j,cond,2} = [];
            ps{j,cond} = [];
            ChoiceDivergence{j,cond} = [];

            if isGroup == 1 %group analysis
                for f = 1:length(file_to_analyze)
                    if cond_num == 3 % pick motion period [0,1500] regardless of dt
                        st = length(list{1,file_to_analyze(f)}.result.unique_stim_type)-2;% number of stim types - 2
                    elseif cond_num == 5
                        st = 1;
                    end



                    for PN = 1:2 % 1-pref, 2-null
                        % (1) PSTH 'without' normalization
                        ys{j,cond,PN} = [ys{j,cond,PN}; list{1,f}.result.PSTH{j,1,1}.ys(2*cond-mod(PN,2),idx_winBeg(j,st):idx_winEnd(j,st))];
                        % (2) PSTH 'with' normalization
                        ys_norm{j,cond,PN} = [ys_norm{j,cond,PN}; list{1,f}.result.PSTH{j,1,1}.ys(2*cond-mod(PN,2),idx_winBeg(j,st):idx_winEnd(j,st))/ys_max(f)];
                    end

                    % (1) ΔPSTH 'without' normalization
                    ys_delta{j,cond} = ys{j,cond,1} - ys{j,cond,2};%ΔPSTH
                    % (2) ΔPSTH 'with' normalization
                    ys_norm_delta{j,cond} = ys_norm{j,cond,1} - ys_norm{j,cond,2};%ΔPSTH

                    % choice divergence
                    try
                        ChoiceDivergence{j,cond} = [ChoiceDivergence{j,cond}; list{1,f}.result.ChoiceDivergence_ALL{1,j}(cond,idx_winBeg(j,st):idx_winEnd(j,st))];
                    end
                end

                % mean and ste of PSTH with normalization
                for PN = 1:2 % 1-pref, 2-null
                    % (1) mean and ste of PSTH 'without' normalization
                    ys_mean{j,cond,PN} = mean(ys{j,cond,PN});% mean
                    ys_ste{j,cond,PN} = std(ys{j,cond,PN},0,1)/sqrt(size(ys{j,cond,PN},1));% ste
                    ys_error{j,cond,PN} = ys_ste{j,cond,PN}*1.96;% 95% CI

                    % (2) mean and ste of PSTH 'with' normalization
                    ys_norm_mean{j,cond,PN} = mean(ys_norm{j,cond,PN},1);% mean
                    ys_norm_ste{j,cond,PN} = std(ys_norm{j,cond,PN},0,1)/sqrt(size(ys_norm{j,cond,PN},1));% ste
                    ys_norm_error{j,cond,PN} = ys_norm_ste{j,cond,PN}*1.96;% 95% CI

                    if isnorm == 0
                        % (1) plot PSTH 'without' normalization
                        if p_value_type == 1
                            % type 1: p_value bar
                            plot(ts{j},ys_mean{j,cond,PN},'Color',symbol_color{cond},'LineStyle',symbol{PN},'linewidth',line_width_curve);hold on;
                        else
                            % type 2: shaded error bar
                            h = shadedErrorBar(ts{j},ys_mean{j,cond,PN},ys_ste{j,cond,PN},{'Color',symbol_color{cond},'LineStyle',symbol{PN}},transparent);
                            set(h.mainLine,'LineWidth',line_width_curve);  hold on;
                        end
                    elseif isnorm == 1
                        % (2) plot PSTH 'with' normalization
                        if p_value_type == 1
                            % type 1: p_value bar
                            plot(ts{j},ys_norm_mean{j,cond,PN},'Color',symbol_color{cond},'LineStyle',symbol{PN},'linewidth',line_width_curve);hold on;
                        else
                            % type 2: shaded error bar
                            h = shadedErrorBar(ts{j},ys_norm_mean{j,cond,PN},ys_norm_ste{j,cond,PN},{'Color',symbol_color{cond},'LineStyle',symbol{PN}},transparent);
                            set(h.mainLine,'LineWidth',line_width_curve);  hold on;
                        end
                    end
                end

                % for plotting p-value bar 
                ys_mean_delta{j,cond} = ys_mean{j,cond,1} - ys_mean{j,cond,2}; % '>0'→ pref,  '<0'→ null 
                for t = 1:length(ts{j})                
                    p_value{j,cond}(t) = malegetest(ys_norm{j,cond,1}(:,t),ys_norm{j,cond,2}(:,t));
                end

            elseif isGroup == 0 % one-cell analysis
                f = file_to_analyze
                if cond_num == 3 % pick motion period [0,1500] regardless of dt
                    st = length(list{1,f}.result.unique_stim_type)-2;% number of stim types - 2
                elseif cond_num == 5
                    st = 1;
                end
                for PN = 1:2 % 1-pref, 2-null
                    % (1) PSTH 'without' normalization
                    ys{j,cond,PN} = list{1,f}.result.PSTH{j,1,1}.raw{2*cond-mod(PN,2)}(:,idx_winBeg(j,st):idx_winEnd(j,st));
                end
                ps{j,cond} = list{1,f}.result.PSTH{j,1,1}.ps(cond,idx_winBeg(j,st):idx_winEnd(j,st));

                for PN = 1:2 % 1-pref, 2-null
                    ys_mean{j,cond,PN} = mean(ys{j,cond,PN});% mean
                    ys_ste{j,cond,PN} = std(ys{j,cond,PN},0,1)/sqrt(size(ys{j,cond,PN},1));% ste
                    h = shadedErrorBar(ts{j},ys_mean{j,cond,PN},ys_ste{j,cond,PN},{'Color',symbol_color{cond},'LineStyle',symbol{PN}},transparent);
                    set(h.mainLine,'LineWidth',line_width_curve);  hold on;
                end
                % (1) ΔPSTH 'without' normalization
                ys_delta{j,cond} = ys_mean{j,cond,1} - ys_mean{j,cond,2};%ΔPSTH
            end
        end

    end
    % p-value indicator
    if isGroup == 0
        y_range = [min(min(get(h_sp_1(1),'ylim')),min(get(h_sp_1(2),'ylim'))) max(max(get(h_sp_1(1),'ylim')),max(get(h_sp_1(2),'ylim')))];
        for j = 1:2;%1-StimOn, 2-SaccadeOn
            axes(ax(1, 1, j))
            for cond = 1:cond_num
                % p-value indicator            
                indicator_pos = zeros(size(ps{j,cond}));
                indicator_pos (ys_delta{j,cond} >= 0) = indicator_pos (ys_delta{j,cond} >= 0) + y_range(2)-(cond-1)*y_range(2)/60; % right -> up
                indicator_pos (ys_delta{j,cond} < 0) = indicator_pos (ys_delta{j,cond} < 0) + (cond-1)*y_range(2)/60;% left -> down

                plot(ts{j}(ps{j,cond}<p_critical),indicator_pos(ps{j,cond}<p_critical),'s','Color',symbol_color{cond},'MarkerFaceColor',symbol_color{cond},'MarkerSize',5);hold on;
            end
        end
        y_range = [0 y_range(2)+3];
    end

    % ---------- PSTH makers ---------
    % uimenufcn(h_f(1),'EditCopyFigure')
    % y_range = [min(min(get(h_sp_1(1),'ylim')),min(get(h_sp_1(2),'ylim'))) max(max(get(h_sp_1(1),'ylim')),max(get(h_sp_1(2),'ylim')))];
    y_range = [0.1 0.7];

    % To plot p-value bar if needed
    if p_value_type == 1
        for j = 1 : 2; % 1 = Stim On, 2 = Saccade On
            if j == 1
                axes(ax(1,1,j));
            else
                axes(ax(1,1,j));
            end

            for cond = 1:cond_num            
                % p-value indicator
                indicator_pos = zeros(size(p_value{j,cond}));
                indicator_pos (ys_mean_delta{j,cond} >= 0) =  y_range(2)-(cond-1)*y_range(2)/60; % right -> up
                indicator_pos (ys_mean_delta{j,cond} < 0) = (cond-1)*y_range(2)/60;% left -> down

                plot(ts{j}(p_value{j,cond}<p_critical),indicator_pos(p_value{j,cond}<p_critical),'s','Color',symbol_color{cond},'MarkerFaceColor',symbol_color{cond},'MarkerSize',5);hold on;

            end
        end
    end

    for j = 1:2;%1-StimOn, 2-StimOff
        axes(ax(1,1,j))
        if j == 1
            if isnorm == 1 && cond_num == 5
                malimalihong('Time to motion onset (ms)','Normalized firing rate','fontsize',font_size)
            elseif isnorm == 0 && cond_num == 5
                malimalihong('Time to motion onset (ms)','Firing rate(Hz)','fontsize',font_size)
            elseif isnorm == 1 && cond_num == 3
                malimalihong('Time to stimulus onset (ms)','Normalized firing rate','fontsize',font_size)
            elseif isnorm == 0 && cond_num ==5
                malimalihong('Time to stimulus onset (ms)','Firing rate(Hz)','fontsize',font_size)
            end
            xlim(x_range{j})
            ylim(y_range)   

    %         set(gca,'tickdir','out')
            set(gca,'ticklength',[0.025 0.025]);

            % vertical lines
            plot([0 0 ] ,[y_range(1) y_range(2) ],'k','linewidth',line_width_axis)
            plot([t_motion_off t_motion_off] ,[y_range(1) y_range(2)],'k','linewidth',line_width_axis)
            % axes lines
            set(gca,'linewidth',line_width_axis);

            % add velocity profile
            load('Z:\Data\Acceleration\system4\sigma3.5dis0.2dur1.5.mat');
            time = -0.5:0.01:2.0;
            acc = V20160121_AccelerometerTest_sigma3_Ch1.values';
            for i = 1:250
                vel(i) = -trapz(time(1):0.01:time(i+1),acc(1:i+1));
            end
            maxy = max(vel);
            miny = min(vel);
            vel_norm_1 = ([vel(1) vel]-miny)/(maxy-miny);%normalized to [0,1]
            vel_norm_2 = vel_norm_1 * (y_range(2)-y_range(1));%normalized to [y_range(1),y_range(2)]
            amplitude = 1/4;
            vel_y_offset = y_range(1)*1.0;
            vel_new = vel_norm_2 * amplitude + vel_y_offset*0.80;
            plot(time*1000 + vel_t_offset,vel_new ,'color' ,[0.7 0.7 0.7],'linewidth',1,'linestyle','-.');

            % acceleration profile
            maxy = max(acc);
            miny = min(acc);
            acc_norm_1 = (-acc(1:251)-miny)/(maxy-miny);%normalized to [0,1]
            acc_norm_2 = acc_norm_1 * (y_range(2)-y_range(1));%normalized to [y_range(1),y_range(2)]
            amplitude = 1/4.8;
            acc_y_offset = y_range(1)*1.18;
            acc_new = acc_norm_2 * amplitude + acc_y_offset;
            index_0 = find(time==(0-vel_t_offset)/1000);
            index_2000 = find(time==(2000-vel_t_offset)/1000);
            plot(time(index_0:index_2000)*1000 + vel_t_offset,acc_new(index_0:index_2000) ,'color' ,[0.7 0.7 0.7],'linewidth',1,'linestyle','-');
        else
            malimalihong('','','fontsize',font_size)
            xlim(x_range{j})
            ylim(y_range)
            set(gca,'xtick',[0])
            set(gca,'XTickLabel',{'Sac'})        
            set(gca,'ytick',[])
    %         set(gca,'tickdir','out')
            % vertical lines
            set(gca,'YColor','white');% make y axis disappear
            line([0 0],y_range,'color','k','linestyle','-','linewidth',line_width_axis)
            % axes lines
            set(gca,'linewidth',line_width_axis);
        end
    end
    title(ax(1,1,1),'PSTH')
    axes(ax(1,1,1))
    for cond = 1:cond_num
        text(100, 0.64 - cond*0.04, [legend_str_dt{cond}],'color',symbol_color{cond},'fontsize',font_size);
    end
    for PN = 1:2
        if cond_num == 3 
            text(500, 0.64 - PN*0.04, [legend_str_PN{PN} ],'color','k','fontsize',font_size);
        else
            text(800, 0.64 - PN*0.04, [legend_str_PN{PN} ],'color','k','fontsize',font_size);
        end
    end
    %% Figure 1.2 ΔPSTH  
    figure(1);
    for j = 1:2;%1-StimOn, 2-SaccadeOn   
        if j == 1
            ax(1,2,j) = axes('position',[0.60 0.60 0.30 0.35]);hold on;
        else
            ax(1,2,j) = axes('position',[0.92 0.60 0.05 0.35]);hold on;
        end

        for  cond = 1:cond_num
            if isnorm == 0
                % Without normalization
                ys_delta_mean{j,cond} = mean(ys_delta{j,cond},1);
                ys_delta_ste{j,cond} = std(ys_delta{j,cond},1)/sqrt(size(ys_delta{j,cond},1));
                ys_delta_error{j,cond} = ys_delta_ste{j,cond}*1.96;
                if p_value_type == 1
                    % type 1: p_value bar
                    plot(ts{j},ys_delta_mean{j,cond},'Color',symbol_color{cond},'LineStyle',symbol{1},'linewidth',line_width_curve);hold on;
                else
                    % type 2: shaded error bar with 95% CI
                    h=shadedErrorBar(ts{j},ys_delta_mean{j,cond},ys_delta_ste{j,cond},{'Color',symbol_color{cond},'LineStyle',symbol{1}},transparent);
                    set(h.mainLine,'LineWidth',line_width_curve);  hold on;
                end
            else
                % With normalization
                ys_norm_delta_mean{j,cond} = mean(ys_norm_delta{j,cond},1);
                ys_norm_delta_ste{j,cond} = std(ys_norm_delta{j,cond},1)/sqrt(size(ys_norm_delta{j,cond},1));
                ys_norm_delta_error{j,cond} = ys_norm_delta_ste{j,cond}*1.96;
                if p_value_type == 1
                    % type 1: p_value bar
                    plot(ts{j},ys_norm_delta_mean{j,cond},'Color',symbol_color{cond},'LineStyle',symbol{1},'linewidth',line_width_curve);hold on;
                else
                    % type 2: shaded error bar with 95% CI
                    h=shadedErrorBar(ts{j},ys_norm_delta_mean{j,cond},ys_norm_delta_ste{j,cond},{'Color',symbol_color{cond},'LineStyle',symbol{1}},transparent);
                    set(h.mainLine,'LineWidth',line_width_curve);  hold on;
                end
            end
        end

    end

    % ---------- ΔPSTH markers  ----------
    % y_range = [min(min(get(h_sp_22(1),'ylim')),min(get(h_sp_22(2),'ylim'))) max(max(get(h_sp_22(1),'ylim')),max(get(h_sp_22(2),'ylim')))];
    y_range = [-0.1 0.4];

    % To plot p-value bar if needed
    if p_value_type == 1
        for j = 1 : 2; % 1 = Stim On, 2 = Saccade On
            if j == 1
                axes(ax(1,2,j));
            else
                axes(ax(1,2,j));
            end

            for cond = 1:cond_num            
                % p-value indicator
                indicator_pos = zeros(size(p_value{j,cond}));
                indicator_pos (ys_mean_delta{j,cond} >= 0) =  y_range(2)-(cond-1)*y_range(2)/40; % right -> up
                indicator_pos (ys_mean_delta{j,cond} < 0) = (cond-1)*y_range(2)/40;% left -> down

                plot(ts{j}(p_value{j,cond}<p_critical),indicator_pos(p_value{j,cond}<p_critical),'s','Color',symbol_color{cond},'MarkerFaceColor',symbol_color{cond},'MarkerSize',5);hold on;

            end
        end
    end

    for j = 1:2; % 1 = StimOn, 2 = StimOff
        axes(ax(1,2,j))
        if j == 1
            if isnorm == 1
                malimalihong('Time to motion onset (ms)',' Δ Normalized PSTH','fontsize',15)
            else
                malimalihong('Time to motion onset (ms)','ΔPSTH','fontsize',15)
            end
            xlim(x_range{j})
            ylim(y_range)
    %         set(gca,'tickdir','out')
            set(gca,'ticklength',[0.025 0.025]);

            % vertical lines
            plot([0 0 ] ,[y_range(1) y_range(2) ],'k','linewidth',line_width_axis)
            plot([t_motion_off t_motion_off] ,[y_range(1) y_range(2)],'k','linewidth',line_width_axis)
            % axes lines
            set(gca,'linewidth',line_width_axis);

            % add velocity profile
            load('Z:\Data\Acceleration\system4\sigma3.5dis0.2dur1.5.mat');
            time = -0.5:0.01:2.0;
            acc = V20160121_AccelerometerTest_sigma3_Ch1.values';
            for i = 1:250
                vel(i) = -trapz(time(1):0.01:time(i+1),acc(1:i+1));
            end
           maxy = max(vel);
            miny = min(vel);
            vel_norm_1 = ([vel(1) vel]-miny)/(maxy-miny);%normalized to [0,1]
            vel_norm_2 = vel_norm_1 * (y_range(2)-y_range(1));%normalized to [y_range(1),y_range(2)]
            amplitude = 1/4;
            vel_y_offset = y_range(1);
            vel_new = vel_norm_2 * amplitude + vel_y_offset*1.2;
            plot(time*1000 + vel_t_offset,vel_new ,'color' ,[0.7 0.7 0.7],'linewidth',1,'linestyle','-.');

            % acceleration profile
            maxy = max(acc);
            miny = min(acc);
            acc_norm_1 = (-acc(1:251)-miny)/(maxy-miny);%normalized to [0,1]
            acc_norm_2 = acc_norm_1 * (y_range(2)-y_range(1));%normalized to [y_range(1),y_range(2)]
            amplitude = 1/5;
            acc_y_offset = y_range(1)*0.85;
            acc_new = acc_norm_2 * amplitude + acc_y_offset;        
            index_0 = find(time==(0-vel_t_offset)/1000);
            index_2000 = find(time==(2000-vel_t_offset)/1000);
            plot(time(index_0:index_2000)*1000 + vel_t_offset,acc_new(index_0:index_2000) ,'color' ,[0.7 0.7 0.7],'linewidth',1,'linestyle','-');


        else
            malimalihong('','','fontsize',15)
            xlim(x_range{j})
            ylim(y_range)
            set(gca,'xtick',[0])
            set(gca,'XTickLabel',{'Sac'})
            set(gca,'ytick',[])        
            set(gca,'ticklength',[0.025 0.025]);
    %         set(gca,'tickdir','out')
            % vertical lines
            set(gca,'YColor','white');% make y axis disappear
            line([0 0],y_range,'color','k','linestyle','-','linewidth',line_width_axis)
            % axes lines
            set(gca,'linewidth',line_width_axis);
        end
    end
    axes(ax(1,2,1))
    for cond = 1:cond_num
        text(200, 0.35- cond*0.03, [legend_str_dt{cond}],'color',symbol_color{cond},'fontsize',font_size);
    end
    title(ax(1,2,1),'ΔPSTH')
end

%% Figure 2. Choice divergence 
if strcmp(FigureName, 'Figure ???')
%% Figure 2.1 Choice divergence 
    try % try if Choice divergence exist
        figure(1);
        for j = 1:2;%1-StimOn, 2-SaccadeOn
            if j == 1
                ax(1,3,j) = axes('position',[0.10 0.10 0.30 0.35]);hold on;
            else
                ax(1,3,j) = axes('position',[0.42 0.10 0.05 0.35]);hold on;
            end

            for  cond = 1:cond_num
                ChoiceDivergence_mean{j,cond} = mean(ChoiceDivergence{j,cond},1);
                ChoiceDivergence_ste{j,cond} = std(ChoiceDivergence{j,cond},1)/sqrt(size(ChoiceDivergence{j,cond},1));
                ChoiceDivergence_error{j,cond} = ChoiceDivergence_ste{j,cond}*1.96;
                if p_value_type == 1
                    % type 1: p_value bar
                    plot(ts{j},ChoiceDivergence_mean{j,cond},'Color',symbol_color{cond},'LineStyle',symbol{1},'linewidth',line_width_curve);hold on;
                else
                    % type 2: shaded error bar with 95% CI
                    h = shadedErrorBar(ts{j},ChoiceDivergence_mean{j,cond},ChoiceDivergence_ste{j,cond},{'Color',symbol_color{cond},'LineStyle',symbol{1}},transparent);
                    set(h.mainLine,'LineWidth',line_width_curve);  hold on;
                end
            end

        end

        % ---------- Choice divergence markers ----------
        axes(ax(1,3,1));
        % legend(str,'location','northwest')
        % y_range = [min(min(get(ax(1,3,1),'ylim')),min(get(ax(1,3,2),'ylim'))) max(max(get(ax(1,3,1),'ylim')),max(get(ax(1,3,2),'ylim')))];
        y_range = [-0.1 0.32];

        % To plot p-value bar if needed
        if p_value_type == 1
            for j = 1 : 2; % 1 = Stim On, 2 = Saccade On
                if j == 1
                    axes(ax(1,3,j));
                else
                    axes(ax(1,3,j));
                end

                for cond = 1:cond_num
                    % p-value indicator
                    indicator_pos = zeros(size(p_value{j,cond}));
                    indicator_pos (ys_mean_delta{j,cond} >= 0) =  y_range(2)-(cond-1)*y_range(2)/40; % right -> up
                    indicator_pos (ys_mean_delta{j,cond} < 0) = (cond-1)*y_range(2)/40;% left -> down

                    plot(ts{j}(p_value{j,cond}<p_critical),indicator_pos(p_value{j,cond}<p_critical),'s','Color',symbol_color{cond},'MarkerFaceColor',symbol_color{cond},'MarkerSize',5);hold on;

                end
            end
        end

        for j = 1:2;%1-StimOn, 2-StimOff
            axes(ax(1,3,j))
            if j == 1
                if cond_num == 5
                    malimalihong('Time to motion onset (ms)','Choice divergence','fontsize',15)
                else
                    malimalihong('Time to stimulus onset (ms)','Choice divergence','fontsize',15)
                end
                xlim(x_range{j})
                ylim(y_range)
                plot([0 0 ] ,[y_range(1) y_range(2) ],'color','k','linewidth',line_width_axis)
                plot([t_motion_off t_motion_off] ,[y_range(1) y_range(2)],'color','k','linewidth',line_width_axis)
                set(gca,'linewidth',line_width_axis);
                set(gca,'ytick',0:0.1:0.3)
    %             set(gca,'tickdir','out')
                set(gca,'ticklength',[0.025 0.025]);

                % add velocity profile
                load('Z:\Data\Acceleration\system4\sigma3.5dis0.2dur1.5.mat');
                time = -0.5:0.01:2.0;
                acc = V20160121_AccelerometerTest_sigma3_Ch1.values';
                for i = 1:250
                    vel(i) = -trapz(time(1):0.01:time(i+1),acc(1:i+1));
                end
                maxy = max(vel);
                miny = min(vel);
                vel_norm_1 = ([vel(1) vel]-miny)/(maxy-miny);%normalized to [0,1]
                vel_norm_2 = vel_norm_1 * (y_range(2)-y_range(1));%normalized to [y_range(1),y_range(2)]
                amplitude = 1/4;
                vel_y_offset = y_range(1);
                vel_new = vel_norm_2 * amplitude + vel_y_offset*1.15;
                plot(time*1000 + vel_t_offset,vel_new ,'color' ,[0.7 0.7 0.7],'linewidth',2,'linestyle','-.');

            else
                malimalihong('','','fontsize',15)
                xlim(x_range{j})
                ylim(y_range)
                plot([0 0 ] ,[y_range(1) y_range(2) ],'k','linewidth',line_width_axis)
                set(gca,'xtick',[0])
                set(gca,'XTickLabel',{'Sac'})
                set(gca,'linewidth',line_width_axis);
                set(gca,'ytick',[])
                ylim(y_range)
                set(gca,'YColor','white');% make y axis disappear
                line([0 0],y_range,'color','k','linestyle','-')
    %             set(gca,'tickdir','out')
                set(gca,'ticklength',[0.025 0.025]);
            end
        end
        axes(ax(1,3,1))
        for cond = 1:cond_num
            text(100, 0.30 - cond*0.03, [legend_str_dt{cond}],'color',symbol_color{cond},'fontsize',font_size);
        end
        title(ax(1,3,1),'Choice divergence','fontsize',font_size)
    end
    
    %% Figure 2.2 Distribution of divergence time
    % % note: Divergence time is defined by the first occurence of a 250-ms window 
    % % in which the IN firing are consistently larger than the OUT firing (p < 0.05). 
    % t_range_new = [t_motion_on + 200 , t_motion_off];% limit the selection range 
    % for d = 1:length(t_dt)
    %     t_winBeg_new(d) = t_range_new(1) + t_dt(d)/2 ;
    %     t_winEnd_new(d) = t_range_new(2) + t_dt(d)/2 ;
    %     
    %     [~, idx_winBeg_new(d)] = min(abs(time_stamp{1} - t_winBeg_new(d)));
    %     [~, idx_winEnd_new(d)] = min(abs(time_stamp{1} - t_winEnd_new(d)));
    % end
    % ts_new = time_stamp{1}(idx_winBeg_new(1):idx_winEnd_new(1));
    % 
    % j = 1;
    % sortInd = 1;
    % for cond = 1:cond_num
    %     ps{cond} = [];
    %     
    %     for f = 1:File_normal_total
    %         if cond_num == 3 % pick motion period [0,1500] regardless of dt
    %             st = length(list{1,f}.result.unique_stim_type)-2;% number of stim types - 2
    %         elseif cond_num == 5
    %             st = 1;
    %         end
    %         
    %         % p-value
    %         temp_ps = list{1,f}.result.PSTH{j,sortInd,1}.ps(cond,idx_winBeg_new(st):idx_winEnd_new(st));
    %         temp_ys_delta = list{1,f}.result.PSTH{j,sortInd,1}.ys(2*cond-1,idx_winBeg_new(j,st):idx_winEnd_new(j,st)) - ...
    %                         list{1,f}.result.PSTH{j,sortInd,1}.ys(2*cond,idx_winBeg_new(j,st):idx_winEnd_new(j,st));
    %         try
    %             % General method: define the time of the first continuous 'stepSize_num' stepSize with significant divergence as divergence time
    %             stepSize_num = 25;
    %             for s = 1:stepSize_num
    %                 if s == 1
    %                     str_temp = ['(temp_ps(1:end) < p_critical) & temp_ys_delta(1:end) > 0'];
    %                 else
    %                     str_temp = [str_temp ' & ([temp_ps(' num2str(s) ':end) ' num2str(ones(1,s-1)) '] < p_critical) '...
    %                                          ' & ([temp_ys_delta(' num2str(s) ':end) ' num2str(zeros(1,s-1)) '] > 0) '];
    %                 end
    %             end
    %             eval(['divergence_time_idx = find(' str_temp ',1);']);
    %             divergence_time(cond,f) = ts_new(divergence_time_idx);
    %         catch 
    %             divergence_time(cond,f) = NaN;
    %         end        
    %         
    %     end
    %     divergence_time_num(cond) = sum(~isnan(divergence_time(cond,:))==1)
    %     divergence_time_mean(cond) = nanmean(divergence_time(cond,:));
    %     divergence_time_ste(cond) = nanstd(divergence_time(cond,:))/sqrt(size(divergence_time,2));
    % end
    % 
    % % ----------  plot Distribution of divergence time ----------
    % figure(1);
    % ax(1,4,1) = axes('position',[0.60 0.05 0.40 0.35]);hold on;
    % 
    % plot(divergence_time,repmat(-[1; 2; 3],1,size(divergence_time,2)),'colo',[0.7 0.7 0.7]);hold on;
    % bar_space = 0.32;
    % bar_half_width = 0.15;
    % for cond = 1:cond_num
    %     plot(divergence_time(cond,:),ones(1,size(divergence_time,2))*(-cond),'colo',symbol_color{cond},'marker','o','markersize',5,'linestyle','none');
    %     
    %     x_left = divergence_time_mean(cond) - divergence_time_ste(cond);
    %     x_right = divergence_time_mean(cond) + divergence_time_ste(cond);
    %     y_mid = (-cond) - bar_space;
    %     y_up = ( - cond) - bar_space + bar_half_width;
    %     y_down = (-cond)  - bar_space - bar_half_width;    
    %     % mean
    %     plot(divergence_time_mean(cond),y_mid,'colo',symbol_color{cond},'marker','.')
    %     % vertical lines for error bar
    %     plot([x_left x_left],[y_down, y_up ],'colo',symbol_color{cond},'linewidth',line_width)
    %     plot([x_right x_right],[y_down, y_up],'colo',symbol_color{cond},'linewidth',line_width)
    %     % horizontal line for error bar
    %     plot([x_left x_right],[y_mid, y_mid],'colo',symbol_color{cond},'linewidth',line_width)
    % end
    % xlim([t_motion_on , t_motion_off])
    % y_range = [-3-1.5*bar_space -0.2];
    % ylim(y_range)
    % box off;
    % set(gca,'YColor','white') 
    % set(gca,'tickdir','out')
    % malimalihong('Time to stimulus onset(ms)','','fontsize',15)
    % 
    % % profile
    %  % add velocity profile
    %  load('Z:\Data\Acceleration\system4\sigma3.5dis0.2dur1.5.mat');
    %  time = -0.5:0.01:2.0;
    %  acc = V20160121_AccelerometerTest_sigma3_Ch1.values';
    %  for i = 1:250
    %      vel(i) = -trapz(time(1):0.01:time(i+1),acc(1:i+1));
    %  end
    %  % velocity profile
    %  maxy = max(vel);
    %  miny = min(vel);
    %  plot(time*1000 + vel_t_offset,([vel(1) vel]-miny)/(maxy-miny)*(y_range(2)-y_range(1))/1.1+y_range(1)-0.3,'color' ,[1 0 0 ],'linewidth',line_width,'linestyle','-.');
    %  % acceleration profile
    %  maxy = max(acc);
    %  miny = min(acc);
    %  plot(time*1000 + vel_t_offset,(-acc(1:251)-miny)/(maxy-miny)*(y_range(2)-y_range(1))/1.25+y_range(1)+0.35,'color' ,[0 0 1],'linewidth',line_width,'linestyle','-.');
    %  
    % % test 
    % for i = 1:cond_num
    %     for j = 1:cond_num
    %         divergence_time_p(i,j) = malegetest(divergence_time(i,:),divergence_time(j,:));
    %     end
    % end
    % disp(['p(Ves,Vis) =  ' num2str(divergence_time_p(1,2))]);
    % disp(['p(Vis,Comb) =  ' num2str(divergence_time_p(2,3))]);
    % disp(['p(Ves,Comb) =  ' num2str(divergence_time_p(1,3))]);
    % 
    % %plot p-value
    % h_sp_23(2) = axes('position',[0.85 0.15 0.10 0.9]);
    % % vertical
    % offset = 0.02;
    % plot([0 0],[-1 - offset, -2 + offset],'k-');hold on;
    % plot([0 0],[-2 - offset, -3 + offset],'k-');
    % plot([1 1],[-1 - offset, -3 + offset],'k-');
    % % horizontal
    % hori_wid = 0.3;
    % plot([-hori_wid 0],[-1 - offset, -1 - offset],'k-');
    % plot([-hori_wid 0],[-2 + offset, -2 + offset],'k-');
    % plot([-hori_wid 0],[-2 - offset, -2 - offset],'k-');
    % plot([-hori_wid 0],[-3 + offset, -3 + offset],'k-');
    % plot([-hori_wid 0] + 1,[-1 - offset, -1 - offset],'k-');
    % plot([-hori_wid 0] + 1,[-3 + offset, -3 + offset],'k-');
    % ylim(y_range)
    % xlim([-1 2])
    % box off
    % set(gca,'XColor','white') 
    % set(gca,'YColor','white') 
    % % text
    % text(0.6,-1.6,'***','rotation',90,'fontsize',15);
    % text(0.6,-2.6,'***','rotation',90,'fontsize',15);
    % text(1.6,-2,'**','rotation',90,'fontsize',15);
    % title(ax(1,4,1),'Distribution of divergence time','fontsize',font_size)
end

%% Figure 3. PSTH sorting by angles
if strcmp(FigureName, 'Figure ??')
%% Figure 3.1 PSTH sorting by angles
    figure(2);clf;
    ScrSize = get(0, 'screensize');
    set(2, 'position', [ScrSize(3)*0.05, ScrSize(4)*0.05, ScrSize(3)*0.9, ScrSize(4)*0.8],'color','white','name','PSTH sorting by angles');
    font_size_2 = 12;
    heading = [-8 -4 -2 -1 +1 +2 +4 +8];
    heading_str = {'- 8' '- 4' '- 2' '- 1' '+1' '+2' '+4' '+8'};

    sortInd = 3;% keep the same with single-file analysis

    p_value_type = 1;
    line_color_temp = colormap('jet');
    color_range = [0.05 0.45;0.55 0.95];
    line_color = line_color_temp([round(linspace(length(line_color_temp)*color_range(1,1),length(line_color_temp)*color_range(1,2),4)),...
        round(linspace(length(line_color_temp)*color_range(2,1),length(line_color_temp)*color_range(2,2),4))],:);

    for f = 1:length(file_to_analyze) % for normalization
        ys_max(f) = 0;
        for cond = 1:cond_num
            ys_max(f) = max(ys_max(f),max(max(list{1,file_to_analyze(f)}.result.PSTH{j,sortInd,cond}.ys))); % for normalization
        end
    end

    for cond = 1:cond_num
        for j = 1:2;% 1 = StimOn, 2 = SaccadeOn
            if j == 1
                ax(2,1,2*(cond-1)+j) = axes('position',[0.07+0.19*(cond-1), 0.55, 0.11, 0.30]);hold on;
            else
                ax(2,1,2*(cond-1)+j) = axes('position',[0.19+0.19*(cond-1), 0.55, 0.018, 0.30]);hold on;
            end

            for ang = 1:8 % heading angles
                ys_per_heading{j,cond,ang} = [];
                ys_per_heading_norm{j,cond,ang} = [];
            end


            for f = 1:File_normal_total
                if cond_num == 3 % pick motion period [0,1500] regardless of dt
                    st = length(list{1,f}.result.unique_stim_type)-2;% number of stim types - 2
                elseif cond_num == 5
                    st = 1;
                end

                for ang = 1:8 % heading angles
                    if list{1,f}.result.PREF == RIGHT % % Preferred RIGHT by default
                        % (1) PSTH 'without' normalization
                        ys_per_heading{j,cond,ang} = [ys_per_heading{j,cond,ang}; list{1,f}.result.PSTH{j,sortInd,cond}.ys(ang,idx_winBeg(j,st):idx_winEnd(j,st))];
                        % (2) PSTH 'with' normalization
                        ys_per_heading_norm{j,cond,ang} = [ys_per_heading_norm{j,cond,ang}; list{1,f}.result.PSTH{j,sortInd,cond}.ys(ang,idx_winBeg(j,st):idx_winEnd(j,st))/ys_max(f)];
                    elseif list{1,f}.result.PREF == LEFT
                        % (1) PSTH 'without' normalization
                        ys_per_heading{j,cond,ang} = [ys_per_heading{j,cond,ang}; list{1,f}.result.PSTH{j,sortInd,cond}.ys(9-ang,idx_winBeg(j,st):idx_winEnd(j,st))];
                        % (2) PSTH 'with' normalization
                        ys_per_heading_norm{j,cond,ang} = [ys_per_heading_norm{j,cond,ang}; list{1,f}.result.PSTH{j,sortInd,cond}.ys(9-ang,idx_winBeg(j,st):idx_winEnd(j,st))/ys_max(f)];
                    end

                    if ang >= 5
                        % (1) PSTH 'without' normalization
                        ys_per_heading_delta{j,cond,ang-4} = mean(ys_per_heading{j,cond,ang}) - mean(ys_per_heading{j,cond,9-ang});%1-smallest, 4-largest
                        % (2) PSTH 'with' normalization
                        ys_per_heading_norm_delta{j,cond,ang-4} = mean(ys_per_heading_norm{j,cond,ang}) - mean(ys_per_heading_norm{j,cond,9-ang});%1-smallest, 4-largest
                    end
                end
            end

            for ang = 1:8 % heading angles
                % (1) mean and ste of PSTH 'without' normalization
                ys_per_heading_mean{j,cond,ang} = mean(ys_per_heading{j,cond,ang});% mean
                ys_per_heading_ste{j,cond,ang} = std(ys_per_heading{j,cond,ang},0,1)/sqrt(size(ys_per_heading{j,cond,ang},1));% ste
                ys_per_heading_error{j,cond,ang} = ys_per_heading_ste{j,cond,ang}*1.96;% 95% CI

                % (2) mean and ste of PSTH 'with' normalization
                ys_per_heading_norm_mean{j,cond,ang} = mean(ys_per_heading_norm{j,cond,ang},1);% mean
                ys_per_heading_norm_ste{j,cond,ang} = std(ys_per_heading_norm{j,cond,ang},0,1)/sqrt(size(ys_per_heading_norm{j,cond,ang},1));% ste
                ys_per_heading_norm_error{j,cond,ang} = ys_per_heading_norm_ste{j,cond,ang}*1.96;% 95% CI

                if isnorm == 0
                    % (1) plot PSTH 'without' normalization
                    if p_value_type == 1
                        % type 1: p_value bar
                        plot(ts{j},ys_per_heading_mean{j,cond,ang},'Color',line_color(ang,:),'linewidth',line_width_curve);hold on;
                    else
                        % type 2: shaded error bar
                        h = shadedErrorBar(ts{j},ys_per_heading_mean{j,cond,ang},ys_per_heading_ste{j,cond,ang},{'Color',line_color(ang,:)},transparent);
                        set(h.mainLine,'LineWidth',line_width_curve);  hold on;
                    end
                else
                    % (2) plot PSTH 'with' normalization
                    if p_value_type == 1
                        % type 1: p_value bar
                        plot(ts{j},ys_per_heading_norm_mean{j,cond,ang},'Color',line_color(ang,:),'linewidth',line_width_curve);hold on;
                    else
                        % type 2: shaded error bar
                        h = shadedErrorBar(ts{j},ys_per_heading_norm_mean{j,cond,ang},ys_per_heading_norm_ste{j,cond,ang},{'Color',line_color(ang,:)},transparent);
                        set(h.mainLine,'LineWidth',line_width_curve);  hold on;
                    end
                end
            end
            end

    end

    % =====  markers =====

    for cond = 1:cond_num    
    %     y_range = [min(min(get(h_sp_4(1,cond),'ylim')),min(get(h_sp_4(2,cond),'ylim'))) max(max(get(h_sp_4(1,cond),'ylim')),max(get(h_sp_4(2,cond),'ylim')))];
    %     y_range = [min(get(ax(2,1,cond),'ylim')) max(get(ax(2,1,cond),'ylim'))];
        if isnorm == 1
            y_range = [0.2 0.6];
        else
            y_range = [10 30];
        end
        for j = 1:2;%1-StimOn, 2-StimOff
            axes(ax(2,1,2*(cond-1)+j));hold on;
            if j == 1
                if isnorm == 1
                    malimalihong('Time to motion onset (ms)','Normalized firing rate','fontsize',font_size_2)
                else
                    malimalihong('Time to motion onset (ms)','Firing rate(Hz)','fontsize',font_size_2)
                end
                xlim(x_range{j})
                ylim(y_range)
                set(gca,'tickdir','out')

                % vertical lines
                plot([0 0 ] ,[y_range(1) y_range(2) ],'k','linewidth',line_width_axis)
                plot([t_motion_off t_motion_off] ,[y_range(1) y_range(2)],'k','linewidth',line_width_axis)
                % axes lines
                set(gca,'linewidth',line_width_axis);

                % add velocity profile
                load('Z:\Data\Acceleration\system4\sigma3.5dis0.2dur1.5.mat');
                time = -0.5:0.01:2.0;
                acc = V20160121_AccelerometerTest_sigma3_Ch1.values';
                for i = 1:250
                    vel(i) = -trapz(time(1):0.01:time(i+1),acc(1:i+1));
                end
                maxy = max(vel);
                miny = min(vel);
                vel_norm_1 = ([vel(1) vel]-miny)/(maxy-miny);%normalized to [0,1]
                vel_norm_2 = vel_norm_1 * (y_range(2)-y_range(1));%normalized to [y_range(1),y_range(2)]
                amplitude = 1/6;
                vel_y_offset = y_range(1);
                vel_new = vel_norm_2 * amplitude + vel_y_offset*0.97;
                plot(time*1000 + vel_t_offset,vel_new ,'color' ,[0.7 0.7 0.7],'linewidth',line_width_vel,'linestyle','-.');
            else
                malimalihong('','','fontsize',font_size_2)
                xlim(x_range{j})
                ylim(y_range)
                set(gca,'xtick',[0])
                set(gca,'XTickLabel',{'Sac'})
                set(gca,'ytick',[])
                set(gca,'tickdir','out')

                % vertical lines
                set(gca,'YColor','white');% make y axis disappear
                line([0 0],y_range,'color','k','linestyle','-','linewidth',line_width_axis)
                % axes lines
                set(gca,'linewidth',line_width_axis);
            end

            if cond == 1 & j == 1
                for ang = 1:8
                    if isnorm == 1
                        text(200, 0.4+ang*0.02,heading_str(ang),'color',line_color(ang,:),'fontsize',font_size_2);
                    else
                        text(200, 22+ang*1.0,heading_str(ang),'color',line_color(ang,:),'fontsize',font_size_2);
                    end
                end
            end

        end
        title(ax(2,1,2*(cond-1)+1),legend_str_dt{cond},'fontsize',font_size_2)
    end

    %% Figure 3.2 Δ PSTH sorting by angles
    figure(2)
    heading_dif_str = {'|1°|' '|2°|' '|4°|' '|8°|'};
    for ang = 1:4 % 1-smallest, 4-largest
        for j = 1:2
            if j == 1
                ax(2,2,2*(ang-1)+j) = axes('position',[0.10+0.23*(ang-1), 0.10, 0.14, 0.30]);hold on;
            else
                ax(2,2,2*(ang-1)+j) = axes('position',[0.25+0.23*(ang-1), 0.10, 0.023, 0.30]);hold on;
            end
            for cond = 1:cond_num
                if isnorm == 0
                    plot(ts{j},ys_per_heading_delta{j,cond,ang},'Color',symbol_color{cond},'linewidth',line_width_curve);hold on;
                else
                    plot(ts{j},ys_per_heading_norm_delta{j,cond,ang},'Color',symbol_color{cond},'linewidth',line_width_curve);hold on;
                end
            end
        end
    end

    for ang = 1:4
        for j = 1:2
            axes(ax(2,2,2*(ang-1)+j));
            %     y_range = [min(min(get(h_sp_4(1,5+ii),'ylim')),min(get(h_sp_4(2,5+ii),'ylim'))) max(max(get(h_sp_4(1,5+ii),'ylim')),max(get(h_sp_4(2,5+ii),'ylim')))];
    %             y_range = [min(get(h_sp_4(1,5+ang),'ylim')) max(get(h_sp_4(1,5+ang),'ylim'))];
            if isnorm == 1
                y_range = [-0.05 0.3];
            else
                y_range = [-1 15];
            end
            if j == 1
                if isnorm == 1
                    malimalihong('Time to motion onset (ms)','ΔNormalized PSTH', 'fontsize',font_size_2)
                else
                    malimalihong('Time to motion onset (ms)','ΔPSTH','fontsize',font_size_2)
                end
                xlim(x_range{j})
                ylim(y_range)
                set(gca,'tickdir','out')

                % vertical lines
                plot([0 0 ] ,[y_range(1) y_range(2) ],'k','linewidth',line_width_axis)
                plot([t_motion_off t_motion_off] ,[y_range(1) y_range(2)],'k','linewidth',line_width_axis)
                % axes lines
                set(gca,'linewidth',line_width_axis);

                % add velocity profile
                load('Z:\Data\Acceleration\system4\sigma3.5dis0.2dur1.5.mat');
                time = -0.5:0.01:2.0;
                acc = V20160121_AccelerometerTest_sigma3_Ch1.values';
                for i = 1:250
                    vel(i) = -trapz(time(1):0.01:time(i+1),acc(1:i+1));
                end
                maxy = max(vel);
                miny = min(vel);
                vel_norm_1 = ([vel(1) vel]-miny)/(maxy-miny);%normalized to [0,1]
                vel_norm_2 = vel_norm_1 * (y_range(2)-y_range(1));%normalized to [y_range(1),y_range(2)]
                amplitude = 1/6;
                vel_y_offset = y_range(1);
                vel_new = vel_norm_2 * amplitude + vel_y_offset*1.1;
                plot(time*1000 + vel_t_offset,vel_new ,'color' ,[0.7 0.7 0.7],'linewidth',line_width_vel,'linestyle','-.');
            else
                malimalihong('','','fontsize',font_size_2)
                xlim(x_range{j})
                ylim(y_range)
                set(gca,'xtick',[0])
                set(gca,'XTickLabel',{'Sac'})
                set(gca,'ytick',[])
                set(gca,'tickdir','out')
                % vertical lines
                set(gca,'YColor','white');% make y axis disappear
                line([0 0],y_range,'color','k','linestyle','-','linewidth',line_width_axis)
                % axes lines
                set(gca,'linewidth',line_width_axis);
            end
            if ang == 1 & j == 1
                for cond = 1:cond_num
                    if isnorm == 1
                        text(200,y_range(2)-cond*0.025,[legend_str_dt{cond} '—'],'color',symbol_color{cond},'fontsize',font_size_2)
                    else
                        text(200,y_range(2)-cond*1,[legend_str_dt{cond} '—'],'color',symbol_color{cond},'fontsize',font_size_2)
                    end
                    end
            end
        end
        title(ax(2,2,2*(ang-1)+1),heading_dif_str(ang),'fontsize',font_size_2)

    end


    % %% Figure 3. Fisher information (linear)
    % %% Figure 3.1.1 Linear Fisher information
    % % FI has been calculated in single file and save in 'result', added by Keo 20180704
    % figure(3);clf;
    % ScrSize = get(0, 'screensize');
    % set(3, 'position', [ScrSize(3)*0.05, ScrSize(4)*0.05, ScrSize(3)*0.9, ScrSize(4)*0.8],'color','white','name','Fisher Information(linear)');
    % 
    % for j = 1:2;%1-StimOn, 2-SaccadeOn
    %         if j == 1
    %             ax(3,1,j) = axes('position',[0.07, 0.50, 0.3, 0.40]);hold on;
    %         else
    %             ax(3,1,j) = axes('position',[0.39, 0.50, 0.06, 0.40]);hold on;
    %         end
    %     for cond = 1:cond_num
    %         FI_linear{j,cond} = [];
    %         
    % %         FI_slope_squared{j,cond} = [];
    % %         FI_var{j,cond} = [];
    %         for f = 1:File_normal_total
    %             if cond_num == 3
    %                 st = length(list{1,f}.result.unique_stim_type)-2;% number of stim types - 2
    %             elseif cond_num == 5
    %                 st = 1;
    %             end
    %             FI_linear{j,cond} = [FI_linear{j,cond} ; list{1,f}.result.FI{1,j}(cond,idx_winBeg(j,st):idx_winEnd(j,st))];
    % 
    % %             FI_slope_squared{j,cond} = [FI_slope_squared{j,cond} ; (list{1,f}.result.FI_slope{1,j}(cond,idx_winBeg(j,st):idx_winEnd(j,st))).^2];
    % %             FI_var{j,cond} = [FI_var{j,cond} ; list{1,f}.result.FI_var{1,j}(cond,idx_winBeg(j,st):idx_winEnd(j,st))];
    %         end
    %         % calculate sum of FI with bootstrap sampling
    %         nboot = 1000;
    %         [bootstat,bootsam] = bootstrp(nboot,@sum,FI_linear{j,cond});
    %         FI_sum_mean{j,cond} = nanmean(bootstat);
    %         FI_sum_ste{j,cond} = nanstd(bootstat,0,1);
    % 
    % %         plot(ts{j},FI_sum_mean{j,cond},'Color',symbol_color{cond},'LineStyle',symbol{1},'linewidth',line_width_0);hold on;
    %         h = shadedErrorBar(ts{j},FI_sum_mean{j,cond},FI_sum_ste{j,cond},{'Color',symbol_color{cond},'LineStyle',symbol{1}},transparent);
    %         set(h.mainLine,'LineWidth',line_width_curve);  hold on;
    %     end
    %     % prediction of FI of combined condition
    % %    plot(ts{j},FI_sum_mean{j,1} + FI_sum_mean{j,2},'Color',[1 0 1],'LineStyle',symbol{1},'linewidth',line_width_0-2);hold on;
    % end
    % 
    % % ---------- axes markers ----------
    % % y_range = [min(min(get(h_sp_3(1),'ylim')),min(get(h_sp_3(2),'ylim'))) max(max(get(h_sp_3(1),'ylim')),max(get(h_sp_3(2),'ylim')))];
    % % y_range = [min(get(ax(3,1,1),'ylim')) max(get(ax(3,1,1),'ylim'))];
    % y_range = [0 5300];
    % for j = 1:2;%1-StimOn, 2-StimOff
    %     axes(ax(3,1,j))
    %     if j == 1
    %         malimalihong('Time to motion onset (ms)','Fisher information(rad ^-^2)','fontsize',15)
    % %         malimalihong('Time to motion onset (ms)','Slope^2','fontsize',15)
    % %         malimalihong('Time to motion onset (ms)','Variance','fontsize',15)
    %         xlim(x_range{j})
    %         ylim(y_range)
    %         plot([0 0 ] ,[y_range(1) y_range(2) ],'color','k','linewidth',line_width_axis)
    %         plot([t_motion_off t_motion_off] ,[y_range(1) y_range(2)],'color','k','linewidth',line_width_axis)
    %         set(gca,'linewidth',line_width_axis);
    % %         set(gca,'ytick',0:0.1:0.3)
    %         set(gca,'tickdir','out')
    % 
    %         % add velocity profile
    %         load('Z:\Data\Acceleration\system4\sigma3.5dis0.2dur1.5.mat');
    %         time = -0.5:0.01:2.0;
    %         acc = V20160121_AccelerometerTest_sigma3_Ch1.values';
    %         for i = 1:250
    %             vel(i) = -trapz(time(1):0.01:time(i+1),acc(1:i+1));
    %         end
    %         maxy = max(vel);
    %         miny = min(vel);
    %         vel_norm_1 = ([vel(1) vel]-miny)/(maxy-miny);%normalized to [0,1]
    %         vel_norm_2 = vel_norm_1 * (y_range(2)-y_range(1));%normalized to [y_range(1),y_range(2)]
    %         amplitude = 1/6;
    %         vel_y_offset = y_range(1);
    %         vel_new = vel_norm_2 * amplitude + vel_y_offset*1.1;
    %         plot(time*1000 + vel_t_offset,vel_new ,'color' ,[0.7 0.7 0.7],'linewidth',line_width_vel,'linestyle','-.');
    %     else
    %         malimalihong('','','fontsize',15)
    %         xlim(x_range{j})
    %         ylim(y_range)
    %         plot([0 0 ] ,[y_range(1) y_range(2) ],'k','linewidth',line_width_axis)
    %         set(gca,'xtick',[0])
    %         set(gca,'XTickLabel',{'Sac'})
    %         set(gca,'linewidth',line_width_axis);
    %         set(gca,'ytick',[])
    %         set(gca,'tickdir','out')
    %         ylim(y_range)
    %         set(gca,'YColor','white');% make y axis disappear
    %     end
    % end
    % title(ax(3,1,1),'Linear FI')
    % 
    % %% Figure 3.1.2  Fisher information of comb-ves and comb-vis 
    % figure(3);
    % ax(3,2,1) = axes('position',[0.02, 0.15, 0.20, 0.20]);hold on;
    % time_stamp_FI = ts;
    % 
    % % find motion period(including pre- and post- 500ms) and its index
    % win_selected = [500 1000];% align to motion onset
    % t_winBeg_selected = win_selected(1) + t_dt/2;
    % t_winEnd_selected = win_selected(2) + t_dt/2;
    % j = 1;%1-StimOn, 2-SaccadeOn
    % for d = 1:length(t_dt)
    %     [~, idx_winBeg_selected(d)] = min(abs(time_stamp_FI{j}-t_winBeg_selected(d)));
    %     [~, idx_winEnd_selected(d)] = min(abs(time_stamp_FI{j}-t_winEnd_selected(d)));
    % end
    % f_selected = [];
    % for f = 1:File_normal_total
    %     d = length(list{1,f}.result.unique_stim_type)-2;% number of stim types - 2
    %     FI_ves(f) = nanmean(FI_linear{j,1}(f,idx_winBeg_selected(d):idx_winEnd_selected(d)));
    %     FI_vis(f) = nanmean(FI_linear{j,2}(f,idx_winBeg_selected(d):idx_winEnd_selected(d)));
    %     FI_comb(f) = nanmean(FI_linear{j,3}(f,idx_winBeg_selected(d):idx_winEnd_selected(d)));
    %     if (FI_comb(f) >= FI_ves(f)) & (FI_comb(f) >= FI_vis(f))
    %         f_selected = [f_selected f];
    %     end
    % end
    % num_f_selected = length(f_selected)
    % % plot and do linear fitting 
    % x = FI_comb - FI_ves; 
    % y = FI_comb - FI_vis; 
    % scatter(x,y,100,'fill','markeredgecolor',[0 0 0],'markerfacecolor',[0 1 0]);
    % scatter(x(f_selected),y(f_selected),100,'fill','markeredgecolor',[0 0 0],'markerfacecolor',[1 0 0]);
    % 
    % format long
    % X = [ones(length(x),1) x'];
    % b = X\y';
    % yCalc = X*b;
    % plot(x,yCalc,'--','linewidth',2,'color','k')
    % Rsq = 1 - sum((y' - yCalc).^2)/sum((y' - mean(y)).^2);
    % text(50,20,['R^2 = ' num2str(Rsq)])
    % 
    % axis equal
    % x_range_2 = get(gca,'xlim');
    % y_range_2 = get(gca,'ylim');
    % x_range_2 = [min(x_range_2(1),y_range_2(1)) max(x_range_2(2),y_range_2(2))];
    % y_range_2 = x_range_2;
    % line([x_range_2(1) x_range_2(2) ],[0 0],'linestyle','-.','color','k')
    % line([0 0],[y_range_2(1)  y_range_2(2)],'linestyle','-.','color','k')
    % line([x_range_2(1) x_range_2(2)],[y_range_2(1)  y_range_2(2)],'linestyle',':','color','k')
    % xlim(x_range_2)
    % ylim(y_range_2)
    % malimalihong('FI(comb - ves)','FI(comb - vis)','fontsize',font_size_2)
    % 
    % %% Figure 3.1.3  Fisher information selected ( comb > ves and comb > vis)  
    % figure(3);
    % ts_FI = ts;
    % for j = 1:size(align_markers,1) % Include two align methods
    %     if j == 1
    %         ax(3,2,2+j) = axes('position',[0.25, 0.15, 0.18, 0.20]);hold on;
    %     else
    %         ax(3,2,2+j) = axes('position',[0.44, 0.15, 0.03, 0.20]);hold on;
    %     end
    %     for cond = 1:cond_num
    %         % calculate sum of FI with bootstrap sampling
    %         nboot = 1000;        
    %         [bootstat,bootsam] = bootstrp(nboot,@sum,FI_linear{j,cond}(f_selected,:));
    %         FI_sum_mean{j,cond} = mean(bootstat);
    %         FI_sum_ste{j,cond} = std(bootstat,0,1);
    %         
    %         h = shadedErrorBar(ts_FI{j},FI_sum_mean{j,cond},FI_sum_ste{j,cond},{'Color',symbol_color{cond},'LineStyle',symbol{1}},transparent);
    %         set(h.mainLine,'LineWidth',line_width_curve);  hold on;
    %     end
    % end
    % % % y_range = [min(min(get(h_sp_3(1),'ylim')),min(get(h_sp_3(2),'ylim'))) max(max(get(h_sp_3(1),'ylim')),max(get(h_sp_3(2),'ylim')))];
    % % y_range = [min(get(ax(3,2,2+1),'ylim')) max(get(ax(3,2,2+1),'ylim'))];
    % y_range = [0 4100];
    % % axes markers
    % for j = 1:2;%1-StimOn, 2-StimOff
    %     axes(ax(3,2,2+j))
    %     if j == 1
    %         malimalihong('Time to motion onset (ms)','Fisher information(rad ^-^2)','fontsize',font_size_2)
    %         %         malimalihong('Time to motion onset (ms)','Slope^2','fontsize',15)
    %         %         malimalihong('Time to motion onset (ms)','Variance','fontsize',15)
    %         xlim(x_range{j})
    %         ylim(y_range)
    %         plot([0 0 ] ,[y_range(1) y_range(2) ],'color','k','linewidth',line_width_axis)
    %         plot([2000 2000] ,[y_range(1) y_range(2)],'color','k','linewidth',line_width_axis)
    %         set(gca,'linewidth',line_width_axis);
    %         %         set(gca,'ytick',0:0.1:0.3)
    %         set(gca,'tickdir','out')
    %         
    %         % add velocity profile
    %         load('Z:\Data\Acceleration\system4\sigma3.5dis0.2dur1.5.mat');
    %         time = -0.5:0.01:2.0;
    %         acc = V20160121_AccelerometerTest_sigma3_Ch1.values';
    %         for i = 1:250
    %             vel(i) = -trapz(time(1):0.01:time(i+1),acc(1:i+1));
    %         end
    %         maxy = max(vel);
    %         miny = min(vel);
    %         vel_norm_1 = ([vel(1) vel]-miny)/(maxy-miny);%normalized to [0,1]
    %         vel_norm_2 = vel_norm_1 * (y_range(2)-y_range(1));%normalized to [y_range(1),y_range(2)]
    %         amplitude = 1/6;
    %         vel_y_offset = y_range(1);
    %         vel_new = vel_norm_2 * amplitude + vel_y_offset*1.1;
    %         plot(time*1000 + vel_t_offset,vel_new ,'color' ,[0.7 0.7 0.7],'linewidth',line_width_vel,'linestyle','-.');        
    %     else
    %         malimalihong('','','fontsize',font_size)
    %         xlim(x_range{j})
    %         ylim(y_range)
    %         plot([0 0 ] ,[y_range(1) y_range(2) ],'k','linewidth',line_width_axis)
    %         set(gca,'xtick',[0])
    %         set(gca,'XTickLabel',{'Sac'})
    %         set(gca,'linewidth',line_width_axis);
    %         set(gca,'ytick',[])
    %         set(gca,'tickdir','out')
    %         ylim(y_range)
    %         set(gca,'YColor','white');% make y axis disappear
    %     end
    % end
end

%% Figure 4. Partial correlation
if strcmp(FigureName, 'Figure ???')
    %% Figure 4.1 Partial correlation cell by cell
    for f = 1:File_normal_total % file/cell
        for cond = 1:cond_num
            cond_idx{f,cond} = find( list{1,f}.result.stim_type_per_trial == cond); % find trial index of some stim type
            X_parr_cell{f,cond} = list{1,f}.result.heading_per_trial(cond_idx{f,cond})';% Heading
            Z_parr_cell{f,cond} = list{1,f}.result.choice_per_trial(cond_idx{f,cond}); % Choice
        end
    end

    time_idx = [0 t_motion_off] + 501;%the 501st element correponse to 'StimOn',1 - 1ms
    for f = 1:File_normal_total % file/cell
        for cond = 1:cond_num
            Y_parr_cell{f,cond} = (sum(list{1,f}.result.spike_aligned{1,1}(cond_idx{f,cond},time_idx(1):time_idx(2)),2))/(time_idx(2)-time_idx(1))*1000; %Response = (sum of spike)/time*1000(Hz)

            % two-way anova for unbalanced design
            p_anova{cond}(f,:) = anovan(Y_parr_cell{f,cond},{X_parr_cell{f,cond},Z_parr_cell{f,cond}},'display','off');

            % partial correlation
            [r_heading_cell(f,cond),p_heading_cell(f,cond)] = partialcorr(X_parr_cell{f,cond},Y_parr_cell{f,cond},Z_parr_cell{f,cond});%heading component
            [r_choice_cell(f,cond),p_choice_cell(f,cond)] = partialcorr(Z_parr_cell{f,cond},Y_parr_cell{f,cond},X_parr_cell{f,cond});%choice component
        end
    end
    h_f(500,7) = figure(500+7);
    set(h_f(500,7),'color','white','position',[10 50 300 400]);clf;
    for cond = 1:cond_num
    %     h_f(500,cond) = figure(500+cond);
    %     set(h_f(500,cond),'color','white','position',[10+270*(cond-1) 50 260 320]);clf;

    %     plot(r_heading_cell(:,cond),r_choice_cell(:,cond),'Color',symbol_color{cond},'marker','o','markersize',5,'linestyle','none');hold on;axis equal;

        % Only cells that showed a significant main effect of heading or choice were included (p < 0.05, two-way ANOVA)
        significance_idx = (p_anova{cond}(:,1)<p_critical) | (p_anova{cond}(:,2)<p_critical);
        display(['condition = ' num2str(cond) ' ,sigificant number = ' num2str(sum(significance_idx))]);
        plot(r_heading_cell(significance_idx,cond),r_choice_cell(significance_idx,cond),'Color',symbol_color{cond},'marker','.','markersize',15,'linestyle','none');hold on;

        malimalihong('Heading partial corr.(R)','Choice partial corr.(R)','fontsize',15);box off;axis equal;
        x_range_2 = [-0.6 0.5];
        y_range_2 = [-0.8 0.8];
        xlim(x_range_2);ylim(y_range_2)
        plot(x_range_2,[0 0],'k-','linewidth',line_width_axis);
        plot([0 0],y_range_2,'k-','linewidth',line_width_axis);
    %     text(x_range_2(1)+(x_range_2(2)-x_range_2(1))/1.6, y_range_2(2)-(y_range_2(2)-y_range_2(1))/20,['n = ', num2str(sum(significance_idx)) ])

    %     %linear fitting
    %     x = r_heading_cell(significance_idx,cond);
    %     y = r_choice_cell(significance_idx,cond);    
    %     format long
    %     X = [ones(length(x),1) x];
    %     b = X\y;
    %     yCalc = X*b;
    %     plot(x,yCalc,'-','linewidth',2,'color','k')
    %     Rsq = 1 - sum((y - yCalc).^2)/sum((y - mean(y)).^2);
    %     text(x_range_2(1)+(x_range_2(2)-x_range_2(1))/1.5,y_range_2(2)-(y_range_2(2)-y_range_2(1))/10,sprintf('R^2 = %1.3f', Rsq))
    end

    % Significant correlation proportion
    for cond = 1:cond_num
        significance_idx = (p_anova{cond}(:,1)<p_critical) | (p_anova{cond}(:,2)<p_critical);
        significance_both_idx{cond} = significance_idx & (p_heading_cell(:,cond) < p_critical) & (p_choice_cell(:,cond) < p_critical);
        significance_heading_idx{cond} = significance_idx & (p_heading_cell(:,cond) < p_critical) & (p_choice_cell(:,cond) >= p_critical);
        significance_choice_idx{cond} = significance_idx & (p_heading_cell(:,cond) >= p_critical) & (p_choice_cell(:,cond) < p_critical);
        bar_data(cond,1:3) = [sum(significance_both_idx{cond}) sum(significance_heading_idx{cond}) sum(significance_choice_idx{cond})]/ File_normal_total;
    end
    h_f(500,6) = figure(500+6);
    set(h_f(500,6),'color','white','position',[310 50 300 400]);clf;
    hArray = bar(bar_data,'stacked');
    set(hArray(1),'Linewidth',2,'EdgeColor','k','Facecolor','k');
    set(hArray(2),'Linewidth',2,'EdgeColor','k','Facecolor',[0.3 0.3 0.3]);
    set(hArray(3),'Linewidth',2,'EdgeColor','k','Facecolor',[0.7 0.7 0.7]);
    if cond_num == 5
        set(gca,'xticklabel',legend_str_dt);
    elseif cond_num == 3
        set(gca,'xticklabel',{'Ves', 'Vis', 'Comb'});
    end
    malimalihong('','Significant correlation proportion','fontsize',15)
    box off;

    %% Figure 4.2 Partial correlation across time
    for f = 1:File_normal_total % file/cell
        for cond = 1:cond_num
            cond_idx{f,cond} = find( list{1,f}.result.stim_type_per_trial == cond); % find trial index of some stim type
            % select specific angles
    %         cond_idx{f,cond} = find( list{1,f}.result.stim_type_per_trial == cond & abs(list{1,f}.result.heading_per_trial)<=2 ); 

            % ================ ★★★★★ ====================
            % Note that heading angles: LEFT is '-', and RIGHT is '+'
            % Thus heading/choice partial correlation: '-' means preferred LEFT,'+' means preferred RIGHT 
            % Therefore, if we simply add them together, on average, the choice partial correlations across cells cancel out

    %         if r_choice_cell(f,cond) > 0 % method 1.2 -- based on overall choice partial correlation
          if r_heading_cell(f,cond) > 0  % method 1.1 -- based on overall heading partial correlation
    %       if list{1,f}.result.PREF == RIGHT % method 2 -- based on CD PREF:

                % Preferred direction is RIGHT, make RIGHT positive
                X_parr{f,cond} = list{1,f}.result.heading_per_trial(cond_idx{f,cond})'; % Heading: keep signs of heading angles the same
                Z_parr{f,cond} = 2 * list{1,f}.result.choice_per_trial(cond_idx{f,cond}) - 3; % Choice :change (L,R) from (1,2) to (-1,1)
            else
                % Preferred direction is LEFT, make LEFT positive
                X_parr{f,cond} = - list{1,f}.result.heading_per_trial(cond_idx{f,cond})'; % change signs of heading angles
                Z_parr{f,cond} = - (2 * list{1,f}.result.choice_per_trial(cond_idx{f,cond}) - 3); % Choice :change (L,R) from (1,2) to (1,-1)
            end
            % ================ ★★★★★ ====================

        end
    end

    for j = 1:2 % align type: 1-StimOn, 2-StimOff
        for f = 1:File_normal_total % file/cell
            for cond = 1:cond_num
                for thresh = 1:length(ts{j})            
                    [~, tt_idx] = min(abs(time_stamp{j}-ts{j}(thresh)));
                    Y_parr{j,cond,f,thresh} = list{1,f}.result.spike_hist{j}(cond_idx{f,cond},tt_idx); %Response 
                    [r_heading{j,cond}(f,thresh),p_heading{j,cond}(f,thresh)] = partialcorr(X_parr{f,cond},Y_parr{j,cond,f,thresh},Z_parr{f,cond});%heading component 
                    [r_choice{j,cond}(f,thresh),p_choice{j,cond}(f,thresh)] = partialcorr(Z_parr{f,cond},Y_parr{j,cond,f,thresh},X_parr{f,cond});%choice component 
                end
            end
        end
    end
    % calulate mean and ste
    for j = 1:2 % align type: 1-StimOn, 2-StimOff
        for cond = 1:cond_num
            for thresh = 1:length(ts{j})
                if issquared == 0
                    % 1. non-squared
                    % (1) heading component
                    r_heading_mean{j,cond}(thresh) = nanmean(r_heading{j,cond}(:,thresh),1);
                    r_heading_ste{j,cond}(thresh) = nanstd(r_heading{j,cond}(:,thresh),0,1)/sqrt(size(r_heading{j,cond}(:,thresh),1));
                    % (2) choice component
                    r_choice_mean{j,cond}(thresh) = nanmean(r_choice{j,cond}(:,thresh),1);
                    r_choice_ste{j,cond}(thresh) = nanstd(r_choice{j,cond}(:,thresh),0,1)/sqrt(size(r_choice{j,cond}(:,thresh),1));
                elseif issquared == 1
                    % 2. squared
                    % (1) heading component
                    r_heading_mean{j,cond}(thresh) = nanmean(r_heading{j,cond}(:,thresh).^2,1);
                    r_heading_ste{j,cond}(thresh) = nanstd(r_heading{j,cond}(:,thresh).^2,0,1)/sqrt(size(r_heading{j,cond}(:,thresh),1));
                    % (2) choice component
                    r_choice_mean{j,cond}(thresh) = nanmean(r_choice{j,cond}(:,thresh).^2,1);
                    r_choice_ste{j,cond}(thresh) = nanstd(r_choice{j,cond}(:,thresh).^2,0,1)/sqrt(size(r_choice{j,cond}(:,thresh),1));
                end
            end
        end
    end

    % 5.3 plot
    for k = 1:2% 1-heading, 2-choice
        h_f(5,k) = figure(50+k);
        set(h_f(5,k),'color','white','position',[50+600*(k-1) 50 600 400]);clf;
        for j = 1:2 % align type: 1-StimOn, 2-StimOff
            if j == 1
                h_sp_5(k,j) = axes('position',[0.15 0.15 0.65 0.8]);hold on;
            else
                h_sp_5(k,j) = axes('position',[0.85 0.15 0.10 0.8]);hold on;
            end
            for cond = 1:cond_num
                if k == 1
                    h = shadedErrorBar(ts{j},r_heading_mean{j,cond},r_heading_ste{j,cond},{'Color',symbol_color{cond}},transparent);
                    set(h.mainLine,'LineWidth',line_width_curve);  hold on;
                elseif k == 2
                    h = shadedErrorBar(ts{j},r_choice_mean{j,cond},r_choice_ste{j,cond},{'Color',symbol_color{cond}},transparent);
                    set(h.mainLine,'LineWidth',line_width_curve);  hold on;
                end
            end
        end
    end

    % 5.4 axes markers
    % y_range = [min([min(get(h_sp_5(1,1),'ylim')),min(get(h_sp_5(1,2),'ylim')),min(get(h_sp_5(2,1),'ylim')),min(get(h_sp_5(2,2),'ylim'))]) ...
    %     max([max(get(h_sp_5(1,1),'ylim')),max(get(h_sp_5(1,2),'ylim')),max(get(h_sp_5(2,1),'ylim')),max(get(h_sp_5(2,2),'ylim'))])];
    y_range = [min([min(get(h_sp_5(1,1),'ylim')),min(get(h_sp_5(2,1),'ylim'))]) ...
        max([max(get(h_sp_5(1,1),'ylim')),max(get(h_sp_5(2,1),'ylim'))])];
    for k = 1:2
        for j = 1:2;%1-StimOn, 2-StimOff
            axes(h_sp_5(k,j))
            if j == 1
                if issquared == 0 
                    if k == 1
                        malimalihong('Time to motion onset (ms)','Heading partial corr.(R)','fontsize',15);
                    else
                        malimalihong('Time to motion onset (ms)','Choice partial corr.(R)','fontsize',15);
                    end
                elseif issquared == 1
                    if k == 1
                        malimalihong('Time to motion onset (ms)','Square of heading partial corr.(R^2)','fontsize',15);
                    else
                        malimalihong('Time to motion onset (ms)','Square of choice partial corr.(R^2)','fontsize',15);
                    end
                end

                xlim(x_range{j})
                ylim(y_range)
                plot([0 0 ] ,[y_range(1) y_range(2) ],'color','k','linewidth',line_width_axis)
                plot([t_motion_off t_motion_off] ,[y_range(1) y_range(2)],'color','k','linewidth',line_width_axis)
                set(gca,'linewidth',line_width_axis+1);
                %         set(gca,'ytick',0:0.1:0.3)
                set(gca,'tickdir','out')

                % add velocity profile
                load('Z:\Data\Acceleration\system4\sigma3.5dis0.2dur1.5.mat');
                time = -0.5:0.01:2.0;
                acc = V20160121_AccelerometerTest_sigma3_Ch1.values';
                for i = 1:250
                    vel(i) = -trapz(time(1):0.01:time(i+1),acc(1:i+1));
                end

                if k == 2
                    % velocity profile
                    maxy = max(vel);
                    miny = min(vel);
                    vel_norm_1 = ([vel(1) vel]-miny)/(maxy-miny);%normalized to [0,1]
                    vel_norm_2 = vel_norm_1 * (y_range(2)-y_range(1));%normalized to [y_range(1),y_range(2)]
                    amplitude = 1/6;
                    vel_y_offset = y_range(1);
                    vel_new = vel_norm_2 * amplitude + vel_y_offset;
                    plot(time*1000 + vel_t_offset,vel_new ,'color' ,[0.7 0.7 0.7],'linewidth',2,'linestyle','-.');
                else
                    % velocity profile
                    maxy = max(vel);
                    miny = min(vel);
                    vel_norm_1 = ([vel(1) vel]-miny)/(maxy-miny);%normalized to [0,1]
                    vel_norm_2 = vel_norm_1 * (y_range(2)-y_range(1));%normalized to [y_range(1),y_range(2)]
                    amplitude = 1/6;
                    vel_y_offset = y_range(1);
                    vel_new = vel_norm_2 * amplitude + vel_y_offset;
                    plot(time*1000 + vel_t_offset,vel_new ,'color' ,[0.7 0.7 0.7],'linewidth',2,'linestyle','-.');

                    % acceleration profile
                    maxy = max(acc);
                    miny = min(acc);
                    acc_norm_1 = (-acc(1:251)-miny)/(maxy-miny);%normalized to [0,1]
                    acc_norm_2 = acc_norm_1 * (y_range(2)-y_range(1));%normalized to [y_range(1),y_range(2)]
                    amplitude = 1/6;
                    acc_y_offset = y_range(1);
                    acc_new = acc_norm_2 * amplitude + acc_y_offset;
                    plot(time*1000 + vel_t_offset,acc_new ,'color' ,[0.7 0.7 0.7],'linewidth',2,'linestyle','-.');
                end
                plot(x_range{j},[0 0],'colo',[0 0 0],'linestyle','-','linewidth',line_width_axis)
            else
                malimalihong('','','fontsize',15)
                xlim(x_range{j})
                ylim(y_range)
                plot([0 0 ] ,[y_range(1) y_range(2) ],'k','linewidth',line_width_axis)
                set(gca,'xtick',[0])
                set(gca,'XTickLabel',{'Sac'})
                set(gca,'linewidth',line_width_axis+1);
                set(gca,'ytick',[])
                set(gca,'tickdir','out')
                ylim(y_range)
                set(gca,'YColor','white');% make y axis disappear
                plot(x_range{j},[0 0],'colo',[0 0 0],'linestyle','-','linewidth',line_width_axis)
            end
        end
    end
end

%% Figure 5. ★Fisher information (neurometric)
if strcmp(FigureName, 'Figure 5A-C') || strcmp(FigureName, 'Figure S5') || ...
        strcmp(FigureName, 'Figure S6') || strcmp(FigureName, 'Figure S7A') || ...
        strcmp(FigureName, 'Figure S7B') || strcmp(FigureName, 'Figure S7C')          

    %% Figure 5.1 Fisher information (neurometric)
    align_type = 1;% align to stim on
    sortInd = 3;% keep the same with single-file analysis

    for cond = 1:cond_num
        FI_neuro_allfile{cond} = [];
        FI_neuro_sum{cond} = [];
        for f = 1:length(file_to_analyze)
            switch list{1,file_to_analyze(f)}.result.FILE(1:3)
                case 'm13'
                    angle_list = [-6 -3 -1.5 -1 1 1.5 3 6]; % quicker
                case 'm10'
                    angle_list = [-8 -4 -2 -1 1 2 4 8];
            end
            FI_neuro_alltime{f,cond} = [];

            for t = 1:length(ts{align_type}) 
                neurometric{f,t,cond} = [];
                fit_data{f,t,cond} = [];
                bias{f,t,cond} = [];


                % extract raw firing rate data sorted by conditions and angles
                for ang = 1:length(angle_list)% heading angles
                    firing_rate{cond,ang} = list{1,file_to_analyze(f)}.result.PSTH{align_type,sortInd,cond}.raw{ang}(:,t)';%列向量转置
                end

                % calculate neurometric function
                for ang = 5:8
                    % make sure AUC > 0.5 by put the input of rocN at the right place
                    if list{1,file_to_analyze(f)}.result.PREF == RIGHT 
                        AUC(ang) = rocN( firing_rate{cond,ang},firing_rate{cond,9-ang}, 100 );% if prefer right, then put right first
                    else
                        AUC(ang) = rocN( firing_rate{cond,9-ang}, firing_rate{cond,ang}, 100 );% if prefer left, then put left first
                    end
                end
                neurometric{f,t,cond} = [1-AUC(8:-1:5) AUC(5:8)];

                % cummulative Gaussian fitting for neurometric function to get threshold
                fit_data{f,t,cond}(:,1) = angle_list';
                fit_data{f,t,cond}(:,2) = neurometric{f,t,cond};
                fit_data{f,t,cond}(:,3) = list{1,file_to_analyze(f)}.result.repetitionN;
                [~,thresh] = cum_gaussfit_max1(fit_data{f,t,cond});
                thresh = thresh * 2; % Britten(1992,JNS) compared with 0°,but we compared between two opposite directions, so multiply 2

                FI_neuro{f,t,cond} = 2/(thresh/180*pi)^2; % convert threshold of neurometric function to FI ,unit: radian (Angelaki, 2011, JNP)
                FI_neuro_alltime{f,cond} = [FI_neuro_alltime{f,cond}, FI_neuro{f,t,cond}];
            end
            FI_neuro_allfile{cond} = [FI_neuro_allfile{cond}; FI_neuro_alltime{f,cond}];
        end
    end
    % bootstrap 1000 times to sum picked FI
    FI_selected = setdiff([1:length(file_to_analyze)],[]);%[8 13 15 19 22 29 30 36 39  41 42]
    for cond = 1:cond_num    
        nboot = 1000;
        [bootstat,bootsam] = bootstrp(nboot,@sum,FI_neuro_allfile{cond}(FI_selected,:));
        FI_neuro_sum_mean{cond} = mean(bootstat);
        FI_neuro_sum_ste{cond} = std(bootstat,0,1);
    end

    % plot overall shuffled FI and at different dt
    figure(5);clf;
    ScrSize = get(0, 'screensize');
    set(5, 'position', [ScrSize(3)*0.05, ScrSize(4)*0.05, ScrSize(3)*0.9, ScrSize(4)*0.8],'color','white','name','Fisher Information(neurometric)');
    time_offset = [0 0; 250 -250; 125 -125;]; %row: comb 0/-500/-250; column: ves/vis
    % all
    ax(5,1,1) = axes('position',[0.08 0.60 0.18 0.35]);hold on;
    for cond = 1:cond_num
        h = shadedErrorBar(ts{align_type},FI_neuro_sum_mean{cond},FI_neuro_sum_ste{cond},{'Color',symbol_color{cond},'LineStyle',symbol{1}},transparent);hold on;
        set(h.mainLine,'LineWidth',line_width_curve);  hold on;
    end
    y_range = get(gca,'ylim');
    malimalihong('Time to stimulus onset(ms)','Fisher information(rad ^-^2)','fontsize',15);hold on;
    xlim(x_range{1});hold on;
    ylims = ylim;
    plot([0 0 ] ,[ylims(1) ylims(2) ],'color','k','linewidth',line_width_axis); hold on;
    plot([t_motion_off t_motion_off] ,[ylims(1) ylims(2)],'color','k','linewidth',line_width_axis);
    set(gca,'linewidth',line_width_axis);
    set(gca,'ticklength',[0.025 0.025]);
    axes(ax(5,1,1))
    for cond = 1:cond_num
        text(100, y_range(2) - cond*y_range(2)/10, [legend_str_dt{cond}],'color',symbol_color{cond},'fontsize',font_size);
    end

    % ves, vis & 0/-250/500
    count = 0;
    for cond = [3 5 4] 
        count = count + 1;
        ax(5,1,1+cond) = axes('position',[0.08+0.24*count 0.60 0.18 0.35]);hold on;
        % ves
        h = shadedErrorBar(ts{align_type}+time_offset(cond-2,1),FI_neuro_sum_mean{1},FI_neuro_sum_ste{1},{'Color',symbol_color{1},'LineStyle',symbol{1}},transparent);hold on;
        set(h.mainLine,'LineWidth',line_width_curve);  hold on;
        % vis
        h = shadedErrorBar(ts{align_type}+time_offset(cond-2,2),FI_neuro_sum_mean{2},FI_neuro_sum_ste{2},{'Color',symbol_color{2},'LineStyle',symbol{1}},transparent);hold on;
        set(h.mainLine,'LineWidth',line_width_curve);  hold on;
        % comb
        h = shadedErrorBar(ts{align_type},FI_neuro_sum_mean{cond},FI_neuro_sum_ste{cond},{'Color',symbol_color{cond},'LineStyle',symbol{1}},transparent);hold on;
        set(h.mainLine,'LineWidth',line_width_curve);  hold on; 

        t_window = [125 1875];
        [~, t_begin] = min(abs(ts{align_type}-t_window(1)));
        [~, t_end] = min(abs(ts{align_type}-t_window(2)));
        [~, ves_begin] = min(abs(ts{align_type}+time_offset(cond-2,1)-t_window(1)));
        [~, ves_end] = min(abs(ts{align_type}+time_offset(cond-2,1)-t_window(2)));
        [~, vis_begin] = min(abs(ts{align_type}+time_offset(cond-2,2)-t_window(1)));
        [~, vis_end] = min(abs(ts{align_type}+time_offset(cond-2,2)-t_window(2)));
        FI_ves_vis_sum = FI_neuro_sum_mean{1}(ves_begin:ves_end) + FI_neuro_sum_mean{2}(vis_begin:vis_end);
        plot(ts{align_type}(t_begin:t_end),FI_ves_vis_sum,'color',[0.5 0.5 0.5],'linewidth',line_width_curve)

        malimalihong('Time to stimulus onset(ms)','','fontsize',15);hold on;
        xlim(x_range{1});hold on;
        ylim(y_range)
        plot([0 0 ] ,[ylims(1) ylims(2) ],'color','k','linewidth',line_width_axis); hold on;
        plot([t_motion_off t_motion_off] ,[ylims(1) ylims(2)],'color','k','linewidth',line_width_axis);
        set(gca,'linewidth',line_width_axis);
        set(gca,'ticklength',[0.025 0.025]);

        text(100, y_range(2) - 1*y_range(2)/10, [legend_str_dt{1}],'color',symbol_color{1},'fontsize',font_size);
        text(100, y_range(2) - 2*y_range(2)/10, [legend_str_dt{2}],'color',symbol_color{2},'fontsize',font_size);
        text(100, y_range(2) - 3*y_range(2)/10, [legend_str_dt{cond}],'color',symbol_color{cond},'fontsize',font_size);
        text(100, y_range(2) - 4*y_range(2)/10, 'Ves+Vis','color',[0.5 0.5 0.5],'fontsize',font_size-5);

    end
end

if strcmp(FigureName, 'Figure ???') 
    %% Figure 5.2 Neurometric funtion
    figure(5)
    t_point = [500 1200 1500 2000];
    y_range = [0.25 0.75];
    for tt = 1:length(t_point)
        ax(5,2,tt) = axes('position',[0.08+0.24*(tt-1) 0.10 0.18 0.35]);hold on;
        title(['Time = ' num2str(t_point(tt)) ' (ms)'])
        [~, t_idx] = min(abs(ts{align_type}-t_point(tt)));
        for cond = 1:cond_num
            neurometric_allfile{tt,cond} = [];
            for f = 1:length(FI_selected)
                neurometric_allfile{tt,cond} = [neurometric_allfile{tt,cond}; neurometric{FI_selected(f),t_idx,cond}];
            end
            neurometric_allfile_mean{tt,cond} = mean(neurometric_allfile{tt,cond});
            neurometric_allfile_ste{tt,cond} = std(neurometric_allfile{tt,cond})/sqrt(length(neurometric_allfile{tt,cond}));
            if cond>0
            errorbar(angle_list,neurometric_allfile_mean{tt,cond},neurometric_allfile_ste{tt,cond},'color',symbol_color{cond},'linewidth',line_width_curve)
            end
        end
        ylim(y_range)
        xlabel('Heading angles(°)')
        ylabel('Proportion correct')
        set(findall(gca,'-property','fontsize'),'fontsize',font_size);
    end



    %% Figure 5.3 3D distribution of Fisher information(too complex to see the detail, so plot in multiple figures)
    for cond = 1:cond_num
        figure(90+cond);clf;
        ScrSize = get(0, 'screensize');
        set(gcf, 'position', [ScrSize(3)*0.05, ScrSize(4)*0.05, ScrSize(3)*0.9, ScrSize(4)*0.8],'color','white','name','Fisher Information(linear)');
        set(gcf,'name',legend_str_dt{cond})
        axes('position',[0.15 0.20 0.7 0.7])

        FI_neuro_allfile_mean{cond} = mean(FI_neuro_allfile{cond}(FI_selected,:),2);
        [~, order_idx] = sort(FI_neuro_allfile_mean{cond},'descend');
        [X,Y] = meshgrid(ts{align_type},1:length(order_idx));
        Z = FI_neuro_allfile{cond}(FI_selected,:);
        surf(X,Y,Z);

        axis ij
        colormap default
        view(gca,[-45 30])

        xlabel('Number of neuron','fontsize',font_size)
        ylabel('Time to stimulus onset(ms)','fontsize',font_size)
        zlabel('Fisher information(rad ^-^2)','fontsize',font_size)
        set(findall(gcf,'-property','fontsize'),'fontsize',20);
    end
end

%% Figure 6. SVM decoder
if strcmp(FigureName, 'Figure 5A-C') 
    % Choice decoder
    % Training 
    f_selected = [1:4 6:37 39:58 60:62 65:74]; % don't choose that less than 15 repetition: 5, 38, 59, 63, 64
    align_type = 1; % align to stim on
    heading = [-8 -4 -2 -1 1 2 4 8];
    train_win = [1000 1200];
    [~,train_win_idx(1)] = min(abs(ts{1}-train_win(1)));
    [~,train_win_idx(2)] = min(abs(ts{1}-train_win(2)));
    cond_selected = [1:5];
    run_num = 10;

    % count number for choosing 
    for cond = cond_selected
        for ha = 1:length(heading)
            for f = 1:length(f_selected)
                ff = f_selected(f);
                train_trial_idx{cond,ha,f} = find( list{1, ff}.result.stim_type_per_trial == cond & ...
                                                   list{1, ff}.result.heading_per_trial  == heading(ha));
                train_trial_num_temp{cond,ha}(f) = length(train_trial_idx{cond,ha,f});
            end
            train_trial_num(cond,ha) = min(train_trial_num_temp{cond,ha});
        end
    end

    for run = 1:run_num % run multiple simulation, and then take the average
        disp(['run = ' num2str(run)])

        % ========== SVM training ==========
        % choose trial in each condition for given times
        clear train_data
        clear train_label
        sample_count = 0;
        for cond = cond_selected
            for ha = 1:length(heading)
                if train_trial_num(cond,ha) > 0
                    for f = 1:length(f_selected)
                        % choose n trial in total N trial
                        temp = randperm(length(train_trial_idx{cond,ha,f}));
                        trial_chosen = train_trial_idx{cond,ha,f}(temp(1:train_trial_num(cond,ha)));
                        for tr = 1:length(trial_chosen)
                            train_data(tr + sample_count, f) = mean(list{1,f_selected(f)}.result.spike_hist{align_type,1}(trial_chosen(tr),train_win_idx(1):train_win_idx(2)))/ys_max(f_selected(f));
                            % training SVM with ground true, but not monkey's behavioral choices
                            if heading(ha) > 0
                                train_label(tr + sample_count) = RIGHT;
                            else
                                train_label(tr + sample_count) = LEFT;
                            end
                        end
                    end
                    sample_count = sample_count + train_trial_num(cond,ha);
                end
            end
        end

        SVMModel = fitcsvm(train_data,train_label);
        weight{run} = SVMModel.Beta;

        % ========== SVM test ==========
        step_size = 50;
        bin_size = 100;
        ts_test = ts{align_type}(1):step_size:ts{align_type}(end)-bin_size;
        for t = 1:length(ts_test)
            test_win = [ts_test(t) ts_test(t) + bin_size];
            [~,test_win_idx(1)] = min(abs(ts{1}-test_win(1)));
            [~,test_win_idx(2)] = min(abs(ts{1}-test_win(2)));

            clear test_data
            clear test_label

            for cond = cond_selected
                for ha = 1:length(heading)
                    sample_count(cond,ha) = 0;
                    for f = 1:length(f_selected)
                        % choose n trial in total N trial
                        temp = randperm(length(train_trial_idx{cond,ha,f}));
                        trial_chosen = train_trial_idx{cond,ha,f}(temp(1:train_trial_num(cond,ha)));
                        for tr = 1:length(trial_chosen)
                            test_data{cond,ha}(tr + sample_count(cond,ha), f) = mean(list{1,f_selected(f)}.result.spike_hist{align_type,1}(trial_chosen(tr),test_win_idx(1):test_win_idx(2)))/ys_max(f_selected(f));
                            if heading(ha) > 0
                                test_label{cond,ha}(tr + sample_count(cond,ha)) = RIGHT;
                            else
                                test_label{cond,ha}(tr + sample_count(cond,ha)) = LEFT;
                            end
                        end
                    end
                    sample_count(cond,ha) = sample_count(cond,ha) + train_trial_num(cond,ha);
                end
            end

            % calcuate correct rate
            correct_trial_count = 0;
            trial_count = 0;
            for cond = cond_selected
                for ha = 1:length(heading)
                    [pred_label_temp, score_temp] = predict(SVMModel, test_data{cond,ha});
                    pred_label{cond,ha} = pred_label_temp;
                    score{cond,ha} = score_temp;

                    % ATTENTION: correct rate here means the correct rate of predicting ground true directions 
                    % but not the correct rate of predicting monkey's choices
                    if heading(ha) > 0
                        CR(run,cond,ha,t) = sum(reshape(pred_label{cond,ha},1,length(pred_label{cond,ha})) == ones(1,length(pred_label{cond,ha}))*RIGHT)/length(pred_label{cond,ha});
                        correct_trial_count = correct_trial_count + sum(reshape(pred_label{cond,ha},1,length(pred_label{cond,ha})) == ones(1,length(pred_label{cond,ha}))*RIGHT);
                    else
                        CR(run,cond,ha,t) = sum(reshape(pred_label{cond,ha},1,length(pred_label{cond,ha})) == ones(1,length(pred_label{cond,ha}))*LEFT)/length(pred_label{cond,ha});
                        correct_trial_count = correct_trial_count + sum(reshape(pred_label{cond,ha},1,length(pred_label{cond,ha})) == ones(1,length(pred_label{cond,ha}))*LEFT);
                    end


                    % The correct rate of predicting monkey's choices is:
                    % CR(run,cond,ha,t) = sum(reshape(pred_label{cond,ha},1,length(pred_label{cond,ha})) == test_label{cond,ha})/length(pred_label{cond,ha});

                    trial_count = trial_count + length(pred_label{cond,ha});
                end
            end
            CR_overall(run,t) = correct_trial_count/trial_count;
        end
    end
    disp('Done!')

    %% Figure 6.1 Correct rate across time
    figure(6);clf;
    ScrSize = get(0, 'screensize');
    set(6, 'position', [ScrSize(3)*0.05, ScrSize(4)*0.05, ScrSize(3)*0.9, ScrSize(4)*0.8],'color','white' );

    ax(6,1,1) = axes('position',[0.05 0.70 0.3 0.23]);cla;
    ax(6,1,2) = axes('position',[0.05 0.40 0.3 0.23]);cla;
    ax(6,1,3) = axes('position',[0.05 0.10 0.3 0.23]);cla;

    y_range = [0.35 1.0];
    % 1.1 overall correct rate
    axes(ax(6,1,1))
    hold on;

    % average across run
    CR_overall_mean = mean(CR_overall);
    CR_overall_ste = std(CR_overall)/sqrt(size(CR_overall,1));
    h = shadedErrorBar(ts_test,CR_overall_mean,CR_overall_ste,{'Color','k'},transparent);
    set(h.mainLine,'LineWidth',line_width_curve);  hold on;

    plot([0 0],y_range,'color','k','linewidth',line_width_axis)
    plot([2000 2000],y_range,'color','k','linewidth',line_width_axis)
    ylim(y_range)
    xlabel('Time to stimulus onset(ms)')
    ylabel('Correct rate')

    %% Figure 6.2 correct rate of each condition
    axes(ax(6,1,2))
    hold on;
    y_range = [0.35 1.0];
    for cond = cond_selected
        CR_heading_temp1 = squeeze(CR(:,cond,:,:));
        CR_heading_temp2 = squeeze(mean(CR_heading_temp1,2)); % average across heading angle
        CR_heading_mean{cond} = mean(CR_heading_temp2,1);% average across run
        CR_heading_ste{cond} = std(CR_heading_temp2)/sqrt(size(CR_heading_temp2,1));

        h = shadedErrorBar(ts_test,CR_heading_mean{cond},CR_heading_ste{cond},{'Color',symbol_color{cond}},transparent);
        set(h.mainLine,'LineWidth',line_width_curve);  hold on;
        % plot(ts_test,CR_heading_mean{cond},'color',symbol_color{cond});
    end

    plot([0 0],y_range,'color','k','linewidth',line_width_axis)
    plot([2000 2000],y_range,'color','k','linewidth',line_width_axis)
    ylim(y_range)
    xlabel('Time to stimulus onset(ms)')
    ylabel('Correct rate')

    %% Figure 6.3 Weight
    axes(ax(6,1,3));cla;
    hold on;
    weight_mean = mean(cell2mat(weight),2)

    [~, idx] = sort(abs(weight_mean),'descend')
    plot(weight_mean(idx),'k-')
    xlabel('Number of neurons')
    ylabel('Weight')
    % for i = 1:10
    % temp = randn(1,690)/2.2;
    % [~, idx] = sort(abs(temp),'descend');
    % plot(0.1:0.1:69,temp(idx),'color',[0.8 0 1],'linestyle','none','marker','.')
    % end

    %% Figure 6.4 correct rate of each heading angles
    hold on;

    line_color_temp = colormap('hsv');
    line_color = line_color_temp([round(linspace(length(line_color_temp)*0.05,length(line_color_temp)*0.15,4)),...
                                  round(linspace(length(line_color_temp)*0.25,length(line_color_temp)*0.35,4))],:);

    for cond = cond_selected
        ax(6,2,ha)= axes('position',[0.40 0.85-(cond-1)*0.2 0.22 0.14]);cla;
        for ha = 1:length(heading)

            CR_cond_temp{cond,ha} = squeeze(CR(:,cond,ha,:));
            CR_cond_mean{cond,ha} = mean(CR_cond_temp{cond,ha},1);% average across run
            CR_cond_ste{cond,ha} = std(CR_cond_temp{cond,ha})/sqrt(size(CR_cond_temp{cond,ha},1));

            h = shadedErrorBar(ts_test,CR_cond_mean{cond,ha},CR_cond_ste{cond,ha},{'Color',line_color(ha,:)},transparent);
            set(h.mainLine,'LineWidth',line_width_curve);  hold on;
            %     plot(ts_test,CR_cond_mean{ha},'color',line_color(ha,:));

            if heading(ha)>=0
                text(1050+(9-ha)*200,0.4,num2str(heading(ha)),'color',line_color(ha,:),'fontsize',15)
            elseif heading(ha)<=-0
                text(1000+ha*200,0.25,num2str(heading(ha)),'color',line_color(ha,:),'fontsize',15)
            end
        end

        box off;
        y_range = [0.15 1.05];
        x_range = [-500 2300];
        plot([0 0],y_range,'color','k','linewidth',line_width_axis)
        plot([2000 2000],y_range,'color','k','linewidth',line_width_axis)
        xlim(x_range)
        ylim(y_range)
        xlabel('Time to stimulus onset(ms)')
        ylabel('Correct rate')
    end

    %% Figure 6.5 Neurometric function
    ax(6,2,2) = axes('position',[0.67 0.10 0.28 0.48]);cla;
    axes(ax(6,2,2));hold on;cla;
    font_size = 13;

    xi = min(heading):0.05:max(heading);
    dt = [0 0 0 -500 -250]
    RT_psy = [1033-400 1234-200 1011-330 1137-250 1064-250];% ves, vis, 0, -500, -250
    test_win_psy = [RT_psy'-100 RT_psy'];

    for cond = cond_selected  
        [~,test_win_idx(1)] = min(abs(ts_test-test_win_psy(cond,1)));
        [~,test_win_idx(2)] = min(abs(ts_test-test_win_psy(cond,2)));
        for ha = 1:length(heading)
            CR_psy_temp1 = squeeze(CR(:,cond, ha, test_win_idx(1):test_win_idx(2)));
            CR_psy_temp2 = mean(CR_psy_temp1,1); % average across run
            if heading(ha)>0                                 
                CR_psy(cond,ha) = mean(CR_psy_temp2)
            else
                CR_psy(cond,ha) = 1- mean(CR_psy_temp2);
            end

        end
        [bb(cond),tt(cond)] = cum_gaussfit_max1([heading' CR_psy(cond,:)' ones(length(heading),1)]);
        plot(xi, cum_gaussfit([bb(cond),tt(cond)], xi),'color',symbol_color{cond},'linewidth',line_width_curve)
        plot(heading, CR_psy(cond,:),'color',symbol_color{cond},'marker','.','markersize',30,'linestyle','none')
        text(-1, 0.5 - 0.07 * cond, sprintf('%7.0f   %5.2f   %5.2f  %5.1f%% %5.0f', dt(cond), bb(cond), tt(cond), 100*mean([1-CR_psy(cond,1:4) CR_psy(cond,5:8)]),RT_psy(cond)),'color',symbol_color{cond},'fontsize',font_size)
    end

    text(-1, 0.5, ['     Δt       μ        σ        CR     RT'],'fontsize',font_size)
    tt_pred = sqrt( (tt(1)^2 * tt(2)^2) / (tt(1)^2 + tt(2)^2) );
    bb_pred = ( bb(1)*tt(2)^2 + bb(2)*tt(1)^2 )/( tt(1)^2 + tt(2)^2 );
    text(-1, 0.5 - 0.07 * 6, sprintf('    Pred    %5.2f   %5.2f ', bb_pred, tt_pred),'color',[0.7 0.7 0.7],'fontsize',font_size)
    xlabel('Heading angel(°)')
    ylabel('Rightward choices proportion')

    axes(ax(6,1,2))
    for cond = 1:cond_num    
    %     plot([RT_psy(cond) RT_psy(cond)],y_range,'color',symbol_color{cond})
        text(RT_psy(cond),0.4,'▲','rotation',180,'color',symbol_color{cond})

    end

    %% Figure 6.6 Weight vs. FI
    ax(6,3,1) = axes('position',[0.720 0.70 0.25 0.23]);cla;
    for cond = cond_selected
        temp = FI_neuro_allfile_mean{cond};
        plot(temp(f_selected),abs(SVMModel.Beta),'color',symbol_color{cond},'linestyle','none','marker','o','markerface',symbol_color{cond})
    end
    xlabel('Averaged Fisher information(rad ^-^2)','fontsize',font_size)
    ylabel('Weight')
    box off
end