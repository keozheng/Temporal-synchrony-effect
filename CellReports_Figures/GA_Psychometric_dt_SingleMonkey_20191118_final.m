function result = GA_Psychometric_dt_SingleMonkey_20180518()
% ------------Keo 20170814--------------
% example:
FolderPath = 'Z:\Data\Tempo\Batch\BatchFile_20171217_Fara_dt_Psycho\';
% Messi: BatchFile_20180515_Messi_dt_Psycho
% Fara: BatchFile_20171217_Fara_dt_Psycho

%% load data in the fold
Files = dir(fullfile(FolderPath,'*.mat'));
File_work_count = 0;
list_1 = load([FolderPath Files(1).name]);

% find monkey's number
find_c = strfind(list_1.result.FILE,'c'); 
monkey_num = str2num(list_1.result.FILE(2:find_c-1));
switch monkey_num
    case 10
        monkey_name = 'Messi';
    case 13
        monkey_name = 'Fara';
    case 12
        monkey_name = 'ZhiZi';
    case 3
        monkey_name = 'ZhiMa';
    case 9
        monkey_name = 'CC';
end
        
% save data in 'list'
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
display(['# Session = ' num2str(File_normal_total)])

%% calculate mean and ste
% To find all of conflict time
conflict_time_list = [];
for ct = 1:File_normal_total
    conflict_time_list = union(conflict_time_list,[list{1,ct}.result.unique_conflict_time]);
end
conflict_time_list = sort(conflict_time_list);
if monkey_num == 9 
    conflict_time_list = conflict_time_list(1:end-1) % for monkey CC, ignor 'dt = +500'
end

for st = 1:3
    for ct = 1:length(conflict_time_list)
        Thresh_psy_raw{ct,st} = [];
        Bias_psy_raw{ct,st} = [];
        correct_rate_raw{ct,st} = [];
        for f = 1:File_normal_total
            ct_index = find(list{1,f}.result.unique_conflict_time==conflict_time_list(ct));
            if ~isempty(ct_index) & ~isempty(list{1,f}.result.Thresh_psy{ct_index,st})
                Thresh_psy_raw{ct,st} = [Thresh_psy_raw{ct,st}  [list{1,f}.result.Thresh_psy{ct_index,st}]];
                Bias_psy_raw{ct,st} = [Bias_psy_raw{ct,st}  [list{1,f}.result.Bias_psy{ct_index,st}]];
                correct_rate_raw{ct,st} = [correct_rate_raw{ct,st}  [list{1,f}.result.correct_rate{ct_index,st}]];
            else
                Thresh_psy_raw{ct,st} = [Thresh_psy_raw{ct,st}  NaN];
                Bias_psy_raw{ct,st} = [Bias_psy_raw{ct,st}  NaN];
                correct_rate_raw{ct,st} = [correct_rate_raw{ct,st}  NaN];
            end
        end
        
        if ~isempty(Thresh_psy_raw{ct,st})
            [index_not_nan{ct,st}] = ~isnan(Thresh_psy_raw{ct,st});
            % mean
            Thresh_psy_mean(ct,st) = mean(Thresh_psy_raw{ct,st}([index_not_nan{ct,st}]));
            Bias_psy_mean(ct,st) = mean(Bias_psy_raw{ct,st}([index_not_nan{ct,st}]));
            correct_rate_mean(ct,st) = mean(correct_rate_raw{ct,st}([index_not_nan{ct,st}]));
            % ste
            Thresh_psy_ste(ct,st) = std(Thresh_psy_raw{ct,st}([index_not_nan{ct,st}]))/sqrt(length([index_not_nan{ct,st}]));
            Bias_psy_ste(ct,st) = std(Bias_psy_raw{ct,st}([index_not_nan{ct,st}]))/sqrt(length([index_not_nan{ct,st}]));
            correct_rate_ste(ct,st) = std(correct_rate_raw{ct,st}([index_not_nan{ct,st}]))/sqrt(length([index_not_nan{ct,st}]));
        end
    end
end
%% predicted threshold and bias
Thresh_pred_raw = [];
Bias_pred_raw = [];
for f = 1:File_normal_total
    Thresh_pred_raw = [Thresh_pred_raw [list{1,f}.result.Thresh_pred]];
    Bias_pred_raw = [Bias_pred_raw [list{1,f}.result.Bias_pred]];
end
Thresh_pred_mean = mean(Thresh_pred_raw);
Bias_pred_mean = mean(Bias_pred_raw);
Thresh_pred_ste = std(Thresh_pred_raw)/sqrt(length(Thresh_pred_raw));
Bias_pred_ste = std(Bias_pred_raw)/sqrt(length(Bias_pred_raw));

%% ---------------------------------------Plot Figures----------------------------------
figure(666); clf;
set(666,'color','white','Position', [63 51 1280 600], 'Name', 'Heading Discrimination with Asynchronical Cues');

% colors and markers
COLOR{1} = [0 0 1]; % blue - vestibular only
COLOR{2} = [1 0 0]; % red - visual only
COLOR{3} = [0.6 0.6 0.6];% grey - predicted value
COLOR{4} = [255 192 0]/255;%[255 134 40]/255;% orange - combined in 'dt=-750'
COLOR{5} = [0 1 1];%[51 204 204]/255;% cyan - combined in 'dt=-500'
COLOR{6} = [1 0 1];%[255 51 204]/255;% magenta - combined in 'dt=-250'
COLOR{7} = [0 1 0];% green - combined in 'dt=0'
COLOR{8} = [0 0 0];% black - combined in 'dt=+250'
COLOR{9} = [153 0 153]/255;% purple


for ct = 1:5 % set the maximum number of conflict time is 5
    symbol_color{ct,1} = COLOR{1};% vestibular in 'dt=0' is always blue
    symbol_color{ct,2} = COLOR{2};% visual in 'dt=0' is always red
    symbol{ct,1} = '.';
    symbol{ct,2} = '.';
    symbol{ct,3} = '.';
end

for ct = 1:length(conflict_time_list)
	symbol_color{ct,3} = COLOR{3+ct};
end
symbol_color{ct+1,3} = COLOR{3};% grey - predicted value

marker_size = 6;
line_width = 3;
line_width_bar = 1;
line_width_errorbar = 1;
font_size = 15;

label_string = {'Ves' 'Vis' };
for ct = 1:length(conflict_time_list)
    label_string{ct+2} = num2str(conflict_time_list(ct));
end
label_string{2+length(conflict_time_list)+1} =  'Pred';

%% 1. performance across sessions
% (1) threshold across sessions
a11 = axes('position',[ 0.07 0.7 0.38 0.25]);
x_range = [0 length([Thresh_psy_raw{ct,st}])+1];
for st = 1:3
    for ct = 1:length(conflict_time_list)
        if ~isempty(Thresh_psy_raw{ct,st}([index_not_nan{ct,st}]))
            plot([Thresh_psy_raw{ct,st}],'marker',symbol{ct,st},'markersize',marker_size,'color',symbol_color{ct,st},'linewidth',line_width);
            hold on;
        end
    end
end
plot(Thresh_pred_raw,'marker','.','markersize',marker_size,'color',[0.8 0.8 0.8],'linewidth',line_width);
malimalihong('','Threshold(бу)')
% legend(gca,label_string,'location',[ 0.45 0.77 0.08 0.1],'fontsize',font_size-2)
xlim(x_range);
y_range = [min([Thresh_psy_raw{:}]) max([Thresh_psy_raw{:}])];
ylim(y_range);
text(x_range(1)+(x_range(2)-x_range(1))*0.35,y_range(1)+(y_range(2)-y_range(1))*0.9,['Monkey ' monkey_name ', n = ' num2str(File_normal_total)],'fontsize',font_size);

% (2) bias across sessions
a12 = axes('position',[ 0.07 0.4 0.38 0.25]);
for st = 1:3
    for ct = 1:length(conflict_time_list)
        if ~isempty(Thresh_psy_raw{ct,st}([index_not_nan{ct,st}]))
            plot([Bias_psy_raw{ct,st}],'marker',symbol{ct,st},'markersize',marker_size,'color',symbol_color{ct,st},'linewidth',line_width);
            hold on;
        end
    end
end
plot(Bias_pred_raw,'marker','.','markersize',marker_size,'color',[0.8 0.8 0.8],'linewidth',line_width);
y_range = [min([Bias_psy_raw{:}]) max([Bias_psy_raw{:}])];
ylim(y_range);
xlim(x_range);
malimalihong('','Bias(бу)')

% (3) Correct rate across sessions
a13 = axes('position',[ 0.07 0.1 0.38 0.25]);
for st = 1:3
    for ct = 1:length(conflict_time_list)
        if ~isempty(Thresh_psy_raw{ct,st}([index_not_nan{ct,st}]))
            plot([correct_rate_raw{ct,st}],'marker',symbol{ct,st},'markersize',marker_size,'color',symbol_color{ct,st},'linewidth',line_width);
            hold on;
        end
    end
end
malimalihong('#Sessions','Correct Rate');
xlim(x_range);
y_range = [min([correct_rate_raw{:}]) max([correct_rate_raw{:}])];
ylim(y_range);

%% summary of performance
% 1. bar
a21 = axes('position',[ 0.58 0.7 0.38 0.25]);
a22 = axes('position',[ 0.58 0.4 0.38 0.25]);
a23 = axes('position',[ 0.58 0.1 0.38 0.25]);

bar_count = 0;
Thresh_psy_mean_bar = [];Thresh_psy_ste_bar = [];Thresh_psy_raw_bar = [];
Bias_psy_mean_bar =[];Bias_psy_ste_bar =[];
correct_rate_mean_bar =[];correct_rate_ste_bar =[];

bar_count = 0;
for st = 1:3
    for ct = 1:length(conflict_time_list)
        if ~isempty(Thresh_psy_raw{ct,st}([index_not_nan{ct,st}]))
            bar_count = bar_count + 1;
            %(3) correct rate
            axes(a23);
            bar(bar_count,correct_rate_mean(ct,st),'facecolor',symbol_color{ct,st},'edgecolor','k','linewidth',line_width_bar);hold on;
%             plot(ones(1,length(correct_rate_raw{ct,st}))*bar_count,correct_rate_raw{ct,st},'linestyle','none','marker','.','markersize',15,'color',symbol_color{ct,st})
            correct_rate_mean_bar = [correct_rate_mean_bar correct_rate_mean(ct,st)];
            correct_rate_ste_bar = [correct_rate_ste_bar correct_rate_ste(ct,st)];
            correct_rate_raw_bar{bar_count} = [correct_rate_raw{ct,st}];
        end
    end
end

bar_count = 0;
for st = 1:3
    for ct = 1:length(conflict_time_list)
        if ~isempty(Thresh_psy_raw{ct,st}([index_not_nan{ct,st}]))
            bar_count = bar_count + 1;
            
            %(1) threshold
            axes(a21);            
            bar(bar_count,Thresh_psy_mean(ct,st),'facecolor',symbol_color{ct,st},'edgecolor','k','linewidth',line_width_bar);hold on;
%             plot(ones(1,length(Thresh_psy_raw{ct,st}))*bar_count,Thresh_psy_raw{ct,st},'linestyle','none','marker','.','markersize',15,'color',symbol_color{ct,st})
            Thresh_psy_mean_bar = [Thresh_psy_mean_bar Thresh_psy_mean(ct,st)];
            Thresh_psy_ste_bar = [Thresh_psy_ste_bar Thresh_psy_ste(ct,st)];
            Thresh_psy_raw_bar{bar_count} = [Thresh_psy_raw{ct,st}];
            
            %(2) bias
            axes(a22);
            bar(bar_count,Bias_psy_mean(ct,st),'facecolor',symbol_color{ct,st},'edgecolor','k','linewidth',line_width_bar);hold on;
%             plot(ones(1,length(Bias_psy_raw{ct,st}))*bar_count,Bias_psy_raw{ct,st},'linestyle','none','marker','.','markersize',15,'color',symbol_color{ct,st})
            Bias_psy_mean_bar = [Bias_psy_mean_bar Bias_psy_mean(ct,st)];
            Bias_psy_ste_bar = [Bias_psy_ste_bar Bias_psy_ste(ct,st)];
            Bias_psy_raw_bar{bar_count} = [Bias_psy_raw{ct,st}];
            
        end
        % put the 'pred' at the 3rd place
%         if bar_count == 2
%             % To add predicted threshold and bias
%             bar_count = bar_count + 1;
%             
%             axes(a21);  
%             bar(bar_count,Thresh_pred_mean,'facecolor',[0.8 0.8 0.8],'edgecolor','k','linewidth',line_width_bar);hold on;
%             Thresh_psy_mean_bar = [Thresh_psy_mean_bar Thresh_pred_mean];
%             Thresh_psy_ste_bar = [Thresh_psy_ste_bar Thresh_pred_ste];
%             Thresh_psy_raw_bar{bar_count} = Thresh_pred_raw;
%                       
%             axes(a22);
%             bar(bar_count,Bias_pred_mean,'facecolor',[0.8 0.8 0.8],'edgecolor','k','linewidth',line_width_bar);hold on;           
%             Bias_psy_mean_bar = [Bias_psy_mean_bar Bias_pred_mean];
%             Bias_psy_ste_bar = [Bias_psy_ste_bar Bias_pred_ste];
%             Bias_psy_raw_bar{bar_count} = Bias_pred_raw;
%         end
    end
end

% put the 'pred' at the last place
% To add predicted threshold and bias
bar_count = bar_count + 1;

axes(a21);
bar(bar_count,Thresh_pred_mean,'facecolor',[0.8 0.8 0.8],'edgecolor','k','linewidth',line_width_bar);hold on;
Thresh_psy_mean_bar = [Thresh_psy_mean_bar Thresh_pred_mean];
Thresh_psy_ste_bar = [Thresh_psy_ste_bar Thresh_pred_ste];
Thresh_psy_raw_bar{bar_count} = Thresh_pred_raw;

axes(a22);
bar(bar_count,Bias_pred_mean,'facecolor',[0.8 0.8 0.8],'edgecolor','k','linewidth',line_width_bar);hold on;
Bias_psy_mean_bar = [Bias_psy_mean_bar Bias_pred_mean];
Bias_psy_ste_bar = [Bias_psy_ste_bar Bias_pred_ste];
Bias_psy_raw_bar{bar_count} = Bias_pred_raw;


bar_total = bar_count ;

% p-value
for b1 = 1:bar_total
    for b2 = 1:bar_total
        Thresh_p(b1,b2) = malegetest(Thresh_psy_raw_bar{b1},Thresh_psy_raw_bar{b2});
        Bias_p(b1,b2) = malegetest(Bias_psy_raw_bar{b1},Bias_psy_raw_bar{b2});
        if b1 < bar_total & b2 < bar_total
            correct_rate_p(b1,b2) = malegetest(correct_rate_raw_bar{b1},correct_rate_raw_bar{b2});
        end
    end
end

% 2. error bar
%(1) threshold ------------------------------------------
axes(a21);
% bar(bar_count,Thresh_pred_mean,'facecolor','none','edgecolor',[0.8 0.8 0.8],'linewidth',line_width_bar);hold on;
errorbar(1:bar_total,Thresh_psy_mean_bar,Thresh_psy_ste_bar ,'color','k','linestyle','none','linewidth',line_width_errorbar);hold on;
set(gca,'ylim',y_range);
% for i = 3:bar_total-1
%     bar_pairs(i-2,1:2) = [bar_total-i+2,bar_total];% bar pair
% end
if bar_count == 4
    bar_pairs = [1 2;3 4;2 3;1 3];
elseif bar_count == 6
    bar_pairs = [1 2;5 6;4 6;3 6];
elseif bar_count == 8
%     bar_pairs = [1 2;5 6;4 5; 3 4; 3 5; 3 6; 3 7; 3 8]; % pred is the 3rd
    bar_pairs = [4 5; 3 4;5 6; 4 6;3 6; 6 7; 6 8; 6 2; 6 1]; % pred is the 3rd

end
star_bar(gca,Thresh_psy_mean_bar,Thresh_psy_ste_bar,Thresh_p,bar_pairs);
malimalihong('','Threshold(бу)')
set(gca,'xtick',1:bar_count)
set(gca,'xticklabel',label_string)
% legend(gca,label_string,'location',[ 0.9 0.8 0.08 0.1],'fontsize',font_size)
% ylim([1 4.5])

%(2) bias---------------------------------------------------------------
axes(a22);
% bar(bar_total,Bias_pred_mean,'facecolor','none','edgecolor',[0.8 0.8 0.8],'linewidth',line_width_bar);hold on;
errorbar(1:bar_total,Bias_psy_mean_bar,Bias_psy_ste_bar ,'color','k','linestyle','none','linewidth',line_width_errorbar);hold on;
set(gca,'ylim',y_range);
% for i = 3:bar_total-1
%     bp(i-2,1:2) = [bar_total-i+2,bar_total];% bar pair
% end
% bp = [1 2;bp];
if bar_count == 4
    bp = [1 2;3 4;2 3;1 3];
elseif bar_count == 8
    bp = [1 2;5 6;4 5;4 6; 3 4; 3 5; 3 6; 3 7; 3 8];
elseif bar_count == 6
    bp = [1 2;5 6;4 6;3 6];
end
star_bar(gca,Bias_psy_mean_bar, Bias_psy_ste_bar, Bias_p, bp)
% set(gca,'ytick',-0.8:0.4:0.4)
set(gca,'xtick',1:bar_total)
set(gca,'xticklabel',label_string)
malimalihong('','Bias(бу)')
% ylim([-1.5 1])

%(3) correct rate----------------------------------------------------------------
axes(a23);
errorbar(1:bar_total-1,correct_rate_mean_bar,correct_rate_ste_bar,'color','k','linestyle','none','linewidth',line_width_errorbar);hold on;
set(gca,'ylim',y_range);
clear bp
% for i = 3:bar_total-2
%     bp(i-2,1:2) = [bar_total-i+1,bar_total-1];% bar pair
% end
% bp = [1 2;3 4;bp];
if bar_count == 4
    bp = [1 2;2 3;1 3];
elseif bar_count == 8
    bp = [1 2;6 7;5 6;4 6;3 6;];
elseif bar_count == 6
    bp = [1 2;3 4;4 5;];
elseif bar_count == 6
    bp = [1 2; 3 6; 4 6; 5 6; 7 6;]
end

star_bar(gca,correct_rate_mean_bar, correct_rate_ste_bar, correct_rate_p,bp)
set(gca,'xtick',1:bar_total-1)
set(gca,'xticklabel',label_string)
malimalihong('Conditions','Correct Rate')
% ylim([0.7 1.5])
% uimenufcn(gcf,'EditCopyFigure')

%% output packed into 'result'
result.Thresh_psy_mean_bar = Thresh_psy_mean_bar;
result.Thresh_psy_ste_bar = Thresh_psy_ste_bar;
result.Thresh_psy_raw_bar = Thresh_psy_raw_bar;
result.Bias_psy_mean_bar = Bias_psy_mean_bar;
result.Bias_psy_ste_bar = Bias_psy_ste_bar;
result.Bias_psy_raw_bar = Bias_psy_raw_bar;
result.correct_rate_mean_bar = correct_rate_mean_bar;
result.correct_rate_ste_bar = correct_rate_ste_bar;
result.correct_rate_raw_bar = correct_rate_raw_bar;


% %% single figure
% bar_count = 0;
% figure(91)
% set(91,'color','white','position',[100 100 1000 500])
% gain = 1.5;
% offset =0.5;
% Thresh_psy_mean_bar_2 = [];
% Thresh_psy_ste_bar_2 = [];
% symbol_color_1{1} = [0 0 1];
% symbol_color_1{2} = [1 0 0];
% symbol_color_1{3} = [145 72 200]/255;
% symbol_color_1{4} = [0 0 0];
% symbol_color_1{5} = [255 192 0]/255;
% symbol_color_1{6} = [47 255 255]/255;
% symbol_color_1{7} = [112 255 10]/255;
% symbol_color_1{8} = [1 0 1];
% 
% 
% for st = 1:3
%     for ct = 1:length(conflict_time_list)
%         if ~isempty(Thresh_psy_raw{ct,st}([index_not_nan{ct,st}]))          
%             
%             if bar_count == 2
%                 bar_count = bar_count + 1;
%                 Thresh_psy_mean_bar_2 = [Thresh_psy_mean_bar_2 Thresh_pred_mean];
%                 Thresh_psy_ste_bar_2 = [Thresh_psy_ste_bar_2 Thresh_pred_ste];
%                 Thresh_psy_raw_bar_2{bar_count} = Thresh_pred_raw;
%                 bar(bar_count*gain-offset,Thresh_pred_mean,'facecolor','none','edgecolor',symbol_color_1{bar_count},'linewidth',line_width_bar);hold on;
%             end
%             
%             bar_count = bar_count + 1;
%             %(1) threshold                       
%             bar(bar_count*gain-offset,Thresh_psy_mean(ct,st),'facecolor','none','edgecolor',symbol_color_1{bar_count},'linewidth',line_width_bar);hold on;
% %             plot(ones(1,length(Thresh_psy_raw{ct,st}))*bar_count,Thresh_psy_raw{ct,st},'linestyle','none','marker','.','markersize',15,'color',symbol_color{ct,st})
%             Thresh_psy_mean_bar_2 = [Thresh_psy_mean_bar_2 Thresh_psy_mean(ct,st)];
%             Thresh_psy_ste_bar_2 = [Thresh_psy_ste_bar_2 Thresh_psy_ste(ct,st)];
%             Thresh_psy_raw_bar_2{bar_count} = [Thresh_psy_raw{ct,st}];
%         end
%     end
% end
% 
% bar_total = bar_count ;
% errorbar((1:bar_total)*gain-offset,Thresh_psy_mean_bar_2,Thresh_psy_ste_bar_2 ,'color','k','linestyle','none','linewidth',line_width_errorbar);hold on;
% 
% % set(gca,'ylim',y_range);
% % for i = 3:bar_total-1
% %     bar_pairs(i-2,1:2) = [bar_total-i+2,bar_total];% bar pair
% % end
% % bar_pairs = [1 2;3 4;4 5;3 5;bar_pairs ]
% % % bar_pairs = [1 2;2 3;1 3;3 4]
% % star_bar(gca,Thresh_psy_mean_bar,Thresh_psy_ste_bar,Thresh_p,bar_pairs);
% 
% malimalihong('','Threshold(бу)')
% set(gca,'xtick',(1:bar_count)*gain-offset)
% label_string = {'Ves' 'Vis','Pred' '-750' '-500' '-250' '0' '+250'}
% set(gca,'xticklabel',label_string)
% % ylim([4 11])
end