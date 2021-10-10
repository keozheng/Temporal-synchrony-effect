clear;clc;
dt = [-750,-500,-250,0,+250];
for k = 1:length(dt)
    data = load(['Z:\Data\Acceleration\System_deltaT_test\20200424_nsig3.5dur2250amp0.12_' num2str(dt(k)) '.mat']);
    
    figure(k);clf;
    set(gcf,'color','white')
    set(gcf,'name', ['Delta T Measurement (dt=' num2str(dt(k)) ')']);

    a1 = subplot(2,1,1);
    a2 = subplot(2,1,2);

    markers = data.V20200424_nsig3_Ch32;
    start_idxs = find(markers.codes == 03); % 03 = IN_FIX_WIN_CD
    Fixation_to_StimulusOnset = 0.1; % Fixation->Visual stimulus, unit: second
    start_times = markers.times(start_idxs);

    interval = 0.01; % unit: second(s)
    time_range = [-0.4, 3.1]; % unit: second(s)
    times = ((time_range(1):interval:time_range(2)) - Fixation_to_StimulusOnset)*1000; % set 'Stim On' as '0', unit: ms

    COLOR{1} = [0 0 1]; % blue
    COLOR{2} = [0.7 0.7 1]; % light blue
    COLOR{3} = [1 0 0]; % red
    COLOR{4} = [1 0.7 0.7]; % light red
    COLOR{5} = [0.7 0.7 0.7]; % gray
    line_width_0 = 1;
    line_width_1 = 5;

    acc_raw = data.V20200424_nsig3_Ch4.values * 9.8; % from 'g' to m/s^2
    vel_raw = data.V20200424_nsig3_Ch3.values;

    span_acc = 100;
    span_vel = 1;
    acc_raw = smooth(acc_raw, span_acc);
    vel_raw = smooth(vel_raw, span_vel);

    first_peak_idxs = [1200 1100 980 800 700]/10 + 50;
    count_left = 0;
    count_right = 0;
    for i = 1:length(start_idxs)
        start_idx = ceil(start_times(i)/interval);    
    %     display(['trial ' num2str(i) ': start = ' num2str(num2str(start_idx + time_range(1)/interval)) ', end = ' num2str(start_idx + time_range(2)/interval)])
        acc_trial(i,:) = acc_raw(start_idx + time_range(1)/interval : start_idx + time_range(2)/interval);
        vel_trial(i,:) = vel_raw(start_idx + time_range(1)/interval : start_idx + time_range(2)/interval);

        if acc_trial(i,first_peak_idxs(k)) < 0 % whether left or right, 'floor(100+(-dt)/10)' ¡Ö first_peak_time
            count_left = count_left + 1;

            acc_trial_left(count_left,:) = acc_trial(i,:);
            axes(a1); plot(times, acc_trial(i,:),'color',COLOR{2},'linestyle','--','linewidth',line_width_0); hold on;

            vel_trial_left(count_left,:) = vel_trial(i,:);
            axes(a2); plot(times, vel_trial(i,:),'color',COLOR{4},'linestyle','--','linewidth',line_width_0); hold on;
        else
            count_right = count_right + 1;
            acc_trial_right(count_right,:) = acc_trial(i,:);
            axes(a1); plot(times, acc_trial(i,:),'color',COLOR{2},'linestyle','-','linewidth',line_width_0); hold on;

            vel_trial_right(count_right,:) = vel_trial(i,:);
            axes(a2); plot(times, vel_trial(i,:),'color',COLOR{4},'linestyle','-','linewidth',line_width_0); hold on;
        end
        
    end
    display(['dt = ', num2str(dt(k)) ', count_left = ', num2str(count_left) ', count_right = ', num2str(count_left)])
    
    % resample data for bootstrap
    for b = 1:1000
        resampled_num = size(acc_trial_left,1);
        resampled_idx = randi([1,resampled_num],resampled_num,1); % resample 'resampled_num' integers from [1, resampled_num]
        acc_trial_left_resampled = acc_trial_left(resampled_idx,:);
        
        resampled_num = size(vel_trial_left,1);
        resampled_idx = randi([1,resampled_num],resampled_num,1); % resample 'resampled_num' integers from [1, resampled_num]
        vel_trial_left_resampled = vel_trial_left(resampled_idx,:);

        resampled_num = size(acc_trial_right,1);
        resampled_idx = randi([1,resampled_num],resampled_num,1); % resample 'resampled_num' integers from [1, resampled_num]
        acc_trial_right_resampled = acc_trial_right(resampled_idx,:);

        resampled_num = size(vel_trial_right,1);
        resampled_idx = randi([1,resampled_num],resampled_num,1); % resample 'resampled_num' integers from [1, resampled_num]
        vel_trial_right_resampled = vel_trial_right(resampled_idx,:);

        % find the peak time of acceleration
        peak_time_acc(b) = (times(find(mean(acc_trial_left_resampled)==max(mean(acc_trial_left_resampled)))) + times(find(mean(acc_trial_right_resampled)==max(mean(acc_trial_right_resampled)))))/2;
%         axes(a1); xlabel('Time to stimulus onset (ms)'); ylabel('Vestibular acceleration (m/s^2)'); box off;
%         plot(times,mean(acc_trial_left_resampled),'color',COLOR{1},'linestyle',':','linewidth',line_width_1);% average curve, left
%         plot(times,mean(acc_trial_right_resampled),'color',COLOR{1},'linestyle','-','linewidth',line_width_1);% average curve, right
%         plot(times,zeros(length(times)),'k'); % horizontal line
%         y_range = get(gca, 'ylim');
%         plot([peak_time_acc peak_time_acc], [y_range(1) y_range(2)],'color',COLOR{1},'linewidth',line_width_1)
%         
        % find the peak time of velocity
        peak_time_vel(b) = (times(find(mean(vel_trial_left_resampled)==max(mean(vel_trial_left_resampled)))) + times(find(mean(vel_trial_right_resampled)==max(mean(vel_trial_right_resampled)))))/2;
%         axes(a2); xlabel('Time to stimulus onset (ms)'); ylabel('Visual voltage (V)'); box off;
%         plot(times,mean(vel_trial_left_resampled),'color',COLOR{3},'linestyle',':','linewidth',line_width_1);% average curve, left
%         plot(times,mean(vel_trial_right_resampled),'color',COLOR{3},'linestyle','-','linewidth',line_width_1);% average curve, right
%         plot(times,zeros(length(times)),'k'); % horizontal line
%         y_range = get(gca, 'ylim');
%         plot([peak_time_vel peak_time_vel], [y_range(1) y_range(2)],'color',COLOR{3},'linewidth',line_width_1)

        dt_measured(k, b) = peak_time_vel(b) - peak_time_acc(b);
    end
    
    dt_measured_mean(k) = mean(dt_measured(k,:));
    dt_measured_ste(k) = std(dt_measured(k,:));
    
end

set(figure(666),'color','white');clf;
plot(dt,dt_measured_mean,'k.','markersize',10);hold on;
errorbar(dt,dt_measured_mean,dt_measured_ste,'Marker','none','LineStyle','none','color','black')
set(gca,'XTick',dt); 
set(gca,'YTick',dt); 
x_range = get(gca,'xlim')
y_range = get(gca,'ylim')
plot(x_range,y_range,'k--')
box off;
malimalihong('¦¤T_e_x_p_e_c_t_e_d (ms)','¦¤T_m_e_a_s_u_r_e_d (ms)')

set(figure(999),'color','white');clf;
bar(dt,dt-dt_measured_mean,'facecolor','none');hold on;
errorbar(dt,dt-dt_measured_mean,dt_measured_ste,'Marker','none','LineStyle','none','color','black')
set(gca,'XTick',dt); 
x_range = get(gca,'xlim');
y_range = get(gca,'ylim');
box off;
malimalihong('¦¤T_e_x_p_e_c_t_e_d (ms) ','¦¤T_e_x_p_e_c_t_e_d - ¦¤T_m_e_a_s_u_r_e_d (ms)')
for k = 1:length(dt)
    [~, p(k)] = ttest(dt(k)-dt_measured(k,:))
end