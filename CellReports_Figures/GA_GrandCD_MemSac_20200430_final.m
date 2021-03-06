% ------------Keo 20171008--------------
clc;clear;
% example:
% =========== ?????????? ===========

%% Memory saccade
%% 1.1 Memory saccade - Messi
FolderPath_1 = 'Z:\Data\Tempo\Batch\BatchFile_20200430_Messi_FEF_MemSac\';

% 0.1 load data in the fold
Files = dir(fullfile(FolderPath_1,'*.mat'));
File_work_count = 0;
list_01 = load([FolderPath_1 Files(1).name]);
for ct = 1:length(Files)
    if ct==1 & isfield(list_01,'errorFiles')
        continue;
    else
        File_work_count = File_work_count+1;
        File_normal_name_1{File_work_count} = Files(ct).name;
        list_1{File_work_count} = load([FolderPath_1 Files(ct).name]);
    end
end
File_normal_total_1 = File_work_count;
display(['1.1 Memory Saccade (Messi), N = ',num2str(File_normal_total_1)])

for f = 1:File_normal_total_1
    DDI(f) = list_1{1,f}.result.DDI;
end

%% 1.2 Memory saccade - Fara
FolderPath_12 = 'Z:\Data\Tempo\Batch\BatchFile_20200510_Fara_FEF_MemSac\';

% 0.1 load data in the fold
Files = dir(fullfile(FolderPath_12,'*.mat'));
File_work_count = 0;
list_01 = load([FolderPath_12 Files(1).name]);
for ct = 1:length(Files)
    if ct==1 & isfield(list_01,'errorFiles')
        continue;
    else
        File_work_count = File_work_count+1;
        File_normal_name_1{File_work_count} = Files(ct).name;
        list_1{File_work_count+File_normal_total_1} = load([FolderPath_12 Files(ct).name]);
    end
end
File_normal_total_12 = File_work_count;
display(['1.2 Memory Saccade (Fara), N = ',num2str(File_normal_total_12)])

for f = 1:File_normal_total_12
    DDI(f+File_normal_total_1) = list_1{1,f+File_normal_total_1}.result.DDI;
end

%% Grand CD(choice divergence)
%% 2.1 Grand CD - Messi
FolderPath_2 = 'Z:\Data\Tempo\Batch\BatchFile_20200506_Messi_FEF_GrandCD\';

% 0.1 load data in the fold
Files = dir(fullfile(FolderPath_2,'*.mat'));
File_work_count = 0;
list_02 = load([FolderPath_2 Files(1).name]);
for ct = 1:length(Files)
    if ct==1 & isfield(list_02,'errorFiles')
        continue;
    else
        File_work_count = File_work_count+1;
        File_normal_name_2{File_work_count} = Files(ct).name;
        list_2{File_work_count} = load([FolderPath_2 Files(ct).name]);
    end
end
File_normal_total_2 = File_work_count;
display(['2.1 Grand Choice Divergence (Messi), N = ',num2str(File_normal_total_2)])

for f = 1:File_normal_total_2
    CD_grand(f,:) = [list_2{1,f}.result.CD_grand];
    CD_p(f,:) = [list_2{1,f}.result.CD_p];
end

stim_type_name = list_2{1,f}.result.stim_type_names;

%% 2.2 Grand CD - Fara
FolderPath_22 = 'Z:\Data\Tempo\Batch\BatchFile_20200510_Fara_GrandCD\';

% 0.1 load data in the fold
Files = dir(fullfile(FolderPath_22,'*.mat'));
File_work_count = 0;
list_02 = load([FolderPath_22 Files(1).name]);
for ct = 1:length(Files)
    if ct==1 & isfield(list_02,'errorFiles')
        continue;
    else
        File_work_count = File_work_count+1;
        File_normal_name_2{File_work_count} = Files(ct).name;
        list_2{File_work_count + File_normal_total_2} = load([FolderPath_22 Files(ct).name]);
    end
end
File_normal_total_22 = File_work_count;
display(['2.2 Grand Choice Divergence (Fara), N = ',num2str(File_normal_total_22)])

for f = 1:File_normal_total_22
    CD_grand(f + File_normal_total_2,:) = [list_2{1,f}.result.CD_grand];
    CD_p(f + File_normal_total_2,:) = [list_2{1,f}.result.CD_p];
end

% ATTENTION! take absolute value of grand choice divergence.
CD_grand = abs(CD_grand);
disp('ATTENTION! take absolute value of grand choice divergence.')

%% memory saccade DDI v.s. grand choice divergence
stim_type_colors = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1];
ScrSize = get(0, 'screensize');
set(figure(666),'pos',[ScrSize(3)*0.01, ScrSize(4)*0.35, ScrSize(3)*0.98, ScrSize(4)*0.30],'name','Raster plot & Grand choice divergence (Grand CD)','color','white'); clf;

for st = 1:5 % stim type
    % === raw data
    subplot('position',[0.08+0.18*(st-1), 0.15, 0.13, 0.65]);hold on;
%     title(stim_type_name(st))
    
    % Messi
    find_p_small = find(CD_p(1:File_normal_total_2,st)<0.05);
    find_p_big = find(CD_p(1:File_normal_total_2,st)>=0.05);    
    plot(DDI(find_p_small), CD_grand(find_p_small,st),'color',stim_type_colors(st,:),'marker','o','markersize',5,'linestyle','none','markerfacecolor',stim_type_colors(st,:))
    plot(DDI(find_p_big), CD_grand(find_p_big,st),'color',stim_type_colors(st,:),'marker','o','markersize',5,'linestyle','none','markerfacecolor','none')
       
    % Fara
    find_p_small = find(CD_p(1+File_normal_total_2:end,st)<0.05);
    find_p_big = find(CD_p(1+File_normal_total_2:end,st)>=0.05);    
    plot(DDI(find_p_small+File_normal_total_2), CD_grand(find_p_small+File_normal_total_2,st),'color',stim_type_colors(st,:),'marker','v','markersize',5,'linestyle','none','markerfacecolor',stim_type_colors(st,:))
    plot(DDI(find_p_big+File_normal_total_2), CD_grand(find_p_big+File_normal_total_2,st),'color',stim_type_colors(st,:),'marker','v','markersize',5,'linestyle','none','markerfacecolor','none')
    
    example_index = 41
    plot(DDI(example_index), CD_grand(example_index,st),'color','k','marker','o','markersize',5,'linestyle','none','markerfacecolor',stim_type_colors(st,:))
        
    box off;
    
%     if st == 1
%         ylabel('Grand CD')
%     elseif st == 3
%         xlabel('Memory Saccade DDI')
%     end
    xlabel('Memory Saccade DDI')
    ylabel('Grand CD')
    
    % linear regression
    x = DDI;
    y = CD_grand(:,st)';  
    n = length(x);
    X = [ones(n,1),x'];
    [b, bint, r, rint, s] = regress(y',X,0.05);
    y_reg = b(1) + b(2) * x;
    plot(x, y_reg,'color','k','linestyle', '-')
    
    R = s(1); p = s(3);
    ylim([0 1])
    x_lim = get(gca, 'xlim');
    y_lim = get(gca, 'ylim');
    
    text(x_lim(1) + diff(x_lim) * 0.05, y_lim(1) + diff(y_lim) * 1.22, ['y = ' num2str(b(1)) ' + ' num2str(b(2)) ' * x'])
    text(x_lim(1) + diff(x_lim) * 0.45, y_lim(1) + diff(y_lim) * 1.16, [stim_type_name(st)])
    text(x_lim(1) + diff(x_lim) * 0.05, y_lim(1) + diff(y_lim) * 1.10, ['p = ' num2str(p) ])
    text(x_lim(1) + diff(x_lim) * 0.05, y_lim(1) + diff(y_lim) * 1.04, ['R = ' num2str(R)  ])
    
    
end