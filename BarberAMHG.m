clear all
close all
clc

plotit=0; % For AMHG plots
plots = 0; % For data anlysis plots

%folder
masterfolder=dir('D:\box sync\Published papers\Caitline\Barber_data\USGS_Final_Rivers\');

for x=4:length(masterfolder)
    
    rivername=masterfolder(x).name;
     clear coeffs allw alld allv allQ
    %select data folder
    root=dir(strcat('D:\box sync\Published papers\Caitline\Barber_data\USGS_Final_Rivers\',rivername,'\*.txt'));
    
    %loop through USGS station text files and pull data 2004 and onward
    for k=1:length(root)
        % count how many text files exist for each river
        txt(x,1)=length(root);
        %
        filename=root(k).name;
        fid = fopen(strcat('D:\box sync\Published papers\Caitline\Barber_data\USGS_Final_Rivers\',rivername,'\',filename));
        format='%*s %*s %*s %s %*s %*s %*s %*s %f %f %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %f %*s %f %*s %*s %*s %*s %*s %*s %*s %*s ';
        data=textscan(fid,format,'HeaderLines',16,'delimiter','\t');
        fclose(fid);
        
        date=data{1};
        depth=data{2};
        discharge=data{3};
        width=data{4};
        velocity=data{5};
        
        numberdates=datenum(date);
        
        
        %pull good values------------------------------
        counter=0;
        for i=1:length(numberdates)
            
            if isnan(depth(i))==0 && isnan(width(i))==0 && isnan(velocity(i))==0 && ...
                    isnan(discharge(i))==0 && depth(i) >0 && width(i) >0 && velocity(i) >0 && discharge(i) >0
                counter=counter+1;
                finalwidth(counter,1)=width(i);
                finaldepth(counter,1)=depth(i);
                finalvelocity(counter,1)=velocity(i);
                finaldischarge(counter,1)=discharge(i);
            end
            
        end
        
        if exist('finalwidth','var')==0
            clear finalwidth finaldepth finalvelocity finaldischarge
            continue
        end
        
         if length(finalwidth)<20
            clear finalwidth finaldepth finalvelocity finaldischarge
            continue
        end
        
        %convert to metric--------------------
        finaldischarge=finaldischarge.*(.3048^3);
        finalwidth=finalwidth.*.3048;
        finaldepth=finaldepth.*.3048;
        finalvelocity=finalvelocity.*.3048;
        %--------------------------------------
        
        
        
        if isempty(finalwidth)
            clear finalwidth finaldepth finalvelocity finaldischarge
            continue
        end
        
        %convert to hydraulic mean depth (via continuity)
        finaldepth=finaldischarge./(finalwidth.*finalvelocity);
        
        
        allw{k}=finalwidth;
        alld{k}=finaldepth;
        allv{k}=finalvelocity;
        allQ{k}=finaldischarge;
        
        n{k}=finaldischarge;
        o{k}=finalwidth;
        
        logn=log10(n{k});
        logo=log10(o{k});
        coeffs(k,1:2)=polyfit(logn,logo,1);
        
        
        p{k}=finaldischarge;
        q{k}=finaldepth;
        
        logp=log10(p{k});
        logq=log10(q{k});
        coeffs(k,3:4)=polyfit(logp,logq,1);
        
        r{k}=finaldischarge;
        s{k}=finalvelocity;
        
        logr=log10(r{k});
        logs=log10(s{k});
        coeffs(k,5:6)=polyfit(logr,logs,1);
        
        holdname{k}=filename;
        
        
        %     plot(logn,logo,'bx')
        %     hold on
        %     plot(logp,logq,'rx')
        %     plot(logr,logs,'kx')
        %
        %
        %
        %     plot(finaldepth)
        %     plot(finalvelocity)
        %     plot(finalwidth)
        %     plot(finaldepth,finalvelocity,'x')
        
        % Calculate average values for each river
        avg_discharge(x,1) = mean(finaldischarge);
        avg_width(x,1) = mean(finalwidth);
        avg_depth(x,1) = mean(finaldepth);
        avg_velocity(x,1) = mean(finalvelocity);
        
        % Calculate min and max discharge
        min_discharge(x,1) = min(finaldischarge);
        max_discharge(x,1) = max(finaldischarge);
        
        % Create histogram of all f exponents for all stations on all rivers
        
        
        
       
        
        
        clear finalwidth finaldepth finalvelocity finaldischarge
        %close all
    end
    
    
    %filter extreme low b (<0.02)
    coeffs(coeffs(:,1)<0.02,:) =0;
    coeffs(coeffs(:,1)<0.00,:)=0;
    coeffs(coeffs(:,1)>1,:)=0;
    coeffs(coeffs(:,3)<0.00,:)=0;
    coeffs(coeffs(:,5)<0.00,:)=0;
    coeffs(coeffs(:,3)>1,:)=0;
    coeffs(coeffs(:,5)>1,:)=0;
    %---------------------------------------
    
    %if there were any stations removed, clean up data
    alld=alld(coeffs(:,1)~=0);
    allv=allv(coeffs(:,1)~=0);
    allw=allw(coeffs(:,1)~=0);
    allQ=allQ(coeffs(:,1)~=0);
    holdname=holdname(coeffs(:,1)~=0);
    coeffs=coeffs(coeffs(:,1)~=0,:);
    %---------------------------------------
    
    %calculate p_int-----------------------
     p_int(x,1)= calc_p_int(allQ,allw);
    
    %---------------------------------------
    
    
    %AMHG----------------------------------------------
    %make the free amhg
    widthAMHG_free=polyfit(coeffs(:,2),coeffs(:,1),1);
    logwc_free(x,1)=-widthAMHG_free(2)/widthAMHG_free(1);
    logQc_wfree(x,1)=-1/widthAMHG_free(1);
    
    
    depthAMHG_free=polyfit(coeffs(:,4),coeffs(:,3),1);
    logdc_free(x,1)=-depthAMHG_free(2)/depthAMHG_free(1);
    logQc_dfree(x,1)=-1/depthAMHG_free(1);
    
    
    velocityAMHG_free=polyfit(coeffs(:,6),coeffs(:,5),1);
    logvc_free(x,1)=-velocityAMHG_free(2)/velocityAMHG_free(1);
    logQc_vfree(x,1)=-1/velocityAMHG_free(1);
    
    %from the depth Qc, calc AMHG. This means it has a fixed slope
    %using definition of linear least squares fit
    
    %using single wc to prove congruence
    velocityAMHG_c(2)=mean(coeffs(:,5))-depthAMHG_free(1)*mean(coeffs(:,6));
    velocityAMHG_c(1)=depthAMHG_free(1);
    logvc_c(x,1)=-velocityAMHG_c(2)/velocityAMHG_c(1);
    
    %use this dc and wc and qc to define velocityAMHG
    %by continuity qc=dc*vc*wc log qc = log dc + log vc + log wc
    logwc_c(x,1) = logQc_dfree(x,1) - logvc_c(x,1) -logdc_free(x,1);
    widthAMHG_c=[depthAMHG_free(1),logwc_c(x,1)/logQc_dfree(x,1)];
    
    Qc(x,1)=10^logQc_dfree(x,1);
    
    vc_free(x,1)=10^logvc_free(x,1);
    vc_c(x,1)=10^logvc_c(x,1);
    
    dc_free(x,1)=10^logdc_free(x,1);
    dc_c(x,1)=10^logdc_free(x,1);
    
    wc_free(x,1)=10^logwc_free(x,1);
    wc_c(x,1)=10^logwc_c(x,1);
   
        
  
    %---------------------------------------------------
       
    if plotit==1
        
        figure;
        avect=-2:0.1:3;
        subplot(1,2,1)
        plot(avect,(avect*widthAMHG_free(1))+widthAMHG_free(2),'b')
        hold on
        depthAMHG_freeLINE = (avect*depthAMHG_free(1))+depthAMHG_free(2);
        plot(avect,depthAMHG_freeLINE,'r')
        plot(avect,(avect*velocityAMHG_free(1))+velocityAMHG_free(2),'k')
        plot(coeffs(:,2),coeffs(:,1),'bx')
        plot(coeffs(:,4),coeffs(:,3),'rx')
        plot(coeffs(:,6),coeffs(:,5),'kx')
        title('Gleason/Smith AMHG')
        xlabel('log a, c, k')
        ylabel('b, f, m')
        legend('width','depth','velocity')
        set(gca,'FontSize',14)
        ylim([0,0.8])
        
        
        subplot(1,2,2)
        plot(avect,(avect*widthAMHG_c(1))+widthAMHG_c(2),'b')
        hold on
        plot(avect,(avect*depthAMHG_free(1))+depthAMHG_free(2),'r')
        plot(avect,(avect*velocityAMHG_c(1))+velocityAMHG_c(2),'k')
        plot(coeffs(:,2),coeffs(:,1),'bx')
        plot(coeffs(:,4),coeffs(:,3),'rx')
        plot(coeffs(:,6),coeffs(:,5),'kx')
        title('Congruent AMHG')
        xlabel('log a, c, k')
        ylabel('b, f, m')
        legend('width','depth','velocity')
        set(gca,'FontSize',14)
        ylim([0,0.8])
        
    end  
    
    % Calculate r squared
    % for depthAMHG_free
     yfit_d = polyval(depthAMHG_free, coeffs(:,4));
     yresid_d = coeffs(:,3) - yfit_d;
     SSresid_d = sum(yresid_d.^2);
     SStotal_d = (length(coeffs(:,3))-1)*var(coeffs(:,3));
     rsq_depthAMHG_free(x,1) = 1 - SSresid_d/SStotal_d;
    
     % for widthAMHG_free 
     yfit_w = polyval(widthAMHG_free, coeffs(:,2));
     yresid_w = coeffs(:,1) - yfit_w;
     SSresid_w = sum(yresid_w.^2);
     SStotal_w = (length(coeffs(:,1))-1)*var(coeffs(:,1));
     rsq_widthAMHG_free(x,1) = 1 - SSresid_w/SStotal_w;
     
     % for velocityhAMHG_free 
     yfit_v = polyval(velocityAMHG_free, coeffs(:,6));
     yresid_v = coeffs(:,5) - yfit_v;
     SSresid_v = sum(yresid_v.^2);
     SStotal_v = (length(coeffs(:,5))-1)*var(coeffs(:,5));
     rsq_velocityAMHG_free(x,1) = 1 - SSresid_v/SStotal_v;
    
     
end

plot(p_int(5:end),rsq_depthAMHG_free,'bo')


logmax_discharge=log10(max_discharge);
logmax_discharge=logmax_discharge(5:end);
logavg_discharge=log10(avg_discharge);
logavg_discharge=logavg_discharge(5:end);

% %
% figure;
% subplot(1,2,1)
% histogram(logavg_discharge)
% %xlabel('log $\bar{Q}$', 'interpreter', 'Latex', 'FontSize', 18)
% xlim([-2 6])
% %ylim([0 60])
% set(gca, 'FontSize', 15)
% hold on
% subplot(1,2,2)
% histogram(logmax_discharge, 'BinWidth', 0.3)
% %xlabel('log $Q_{peak}$', 'interpreter', 'Latex','FontSize', 18)
% xlim([0 5])
% %ylim([0 60])
% set(gca, 'FontSize', 15)


% Delete unecessary rows
% delete first 4 rows of rsq values
Qc=Qc(5:end);
wc_free=wc_free(5:end);
dc_free=dc_free(5:end);
vc_free=vc_free(5:end);

wc_c=wc_c(5:end);
dc_c=dc_c(5:end);
vc_c=vc_c(5:end);

avg_discharge = avg_discharge(5:end);

rsq_depthAMHG_free=rsq_depthAMHG_free(5:end);
rsq_widthAMHG_free=rsq_widthAMHG_free(5:end);
rsq_velocityAMHG_free=rsq_velocityAMHG_free(5:end);
avg_depth = avg_depth(5:end);
avg_width = avg_width(5:end);
avg_velocity = avg_velocity(5:end);

min_discharge = min_discharge(5:end);
max_discharge = max_discharge(5:end);

% Calculate Froude number (Fr_c)
Fr_c_free = vc_free./((9.81.*dc_free).^0.5);

%plot(Qc(Qc<10^4),'bo')

%now that we're here, we need to plot AMHG vs .....


% save x y locations
% read in lat/long from excel file
% FinalMetaData_rsq = xlsread('/Users/aleanderkamb/Documents/AMHG/FinalMetaData_rsq.xlsx');
% FinalMetaData_rsq = horzcat(FinalMetaData_rsq, rsq_depthAMHG_free, rsq_widthAMHG_free, rsq_velocityAMHG_free);   %Combine AMHG values with river metadata;
% FinalMetaData_GIS = [FinalMetaData_rsq Qc]; % Combine with Qc

% Without including lat/long info
Final_data = horzcat(rsq_depthAMHG_free, rsq_widthAMHG_free, rsq_velocityAMHG_free);

% Add Qc
Final_data = [Final_data Qc];
% Add froude number to dataframe file
Final_data = [Final_data Fr_c_free];
% Add average discharge to dataframe file
Final_data = [Final_data avg_discharge];
% Add average width, depth, velocity
Final_data = [Final_data avg_width];
Final_data = [Final_data avg_depth];
Final_data = [Final_data avg_velocity];
% Add congruent width, depth, velocity;
Final_data = [Final_data wc_free];
Final_data = [Final_data dc_free];
Final_data = [Final_data vc_free];
% Add min and max discharge
Final_data = [Final_data min_discharge];
Final_data = [Final_data max_discharge];

% Create matrix of numbers 1:191
river_number=reshape(1:190, [190 1]);
% Add to Final_data
Final_data = [Final_data river_number];

% Create matrix with only rivers rsq < 0.6
lowrsq = Final_data;
LR = lowrsq(:,1)>0.6;
lowrsq(LR,:)=[];

rsq_depthAMHG_free_LOW = lowrsq(:,1);
Qc_LOW = lowrsq(:,4);
Fr_c_free_LOW = lowrsq(:,5);
avg_discharge_LOW = lowrsq(:,6);
avg_width_LOW = lowrsq(:,7);
avg_depth_LOW = lowrsq(:,8);
avg_velocity_LOW = lowrsq(:,9);
wc_free_LOW = lowrsq(:,10);
dc_free_LOW = lowrsq(:,11);
vc_free_LOW = lowrsq(:,12);
min_discharge_r_LOW = lowrsq(:,13);
max_discharge_r_LOW = lowrsq(:,14);
river_number_LOW = lowrsq(:,15);

% Create matrix with only rivers with rsq > 0.6 (column 1 in Final_data)
TF = Final_data(:,1)<0.6;
Final_data(TF,:)=[];

rsq_depthAMHG_free_r = Final_data(:,1);
Qc_r = Final_data(:,4);
Fr_c_free_r = Final_data(:,5);
avg_discharge_r = Final_data(:,6);
avg_width_r = Final_data(:,7);
avg_depth_r = Final_data(:,8);
avg_velocity_r = Final_data(:,9);
wc_free_r = Final_data(:,10);
dc_free_r = Final_data(:,11);
vc_free_r = Final_data(:,12);
min_discharge_r = Final_data(:,13);
max_discharge_r = Final_data(:,14);
river_number_r = Final_data(:,15);

% Create excel file to import into GIS
%xlswrite('GIS.xlsx' , GIS_data)

% Calculate log(Qc) and log(avg_discharge)
logQc = log10(Qc);

logavg_depth = log10(avg_depth);
logavg_width = log10(avg_width);
logavg_velocity = log10(avg_velocity);

logFr_c = log10(Fr_c_free);
logwc = log10(wc_free);
logvc = log10(vc_free);
logdc = log10(dc_free);


% Determine if Qc is between min and max discharge for rsq > 0.6 
% 1 = Qc is within the range of Q
for i=1:length(Final_data);
    Final_data(i,16) = Final_data(i,4)>=Final_data(i,13) && Final_data(i,4)<=Final_data(i,14);
end

% Determine if Qc is between min and max discharge for rsq < 0.6
for i=1:length(lowrsq);
    lowrsq(i,16) = lowrsq(i,4)>=lowrsq(i,13) && lowrsq(i,4)<=lowrsq(i,14);
end

% Combine Final_data and lowrsq 
All_Rivers_Qrange_1 = vertcat(Final_data, lowrsq);
All_Rivers_Qrange_2 = vertcat(Final_data, lowrsq);

QIR = All_Rivers_Qrange_1(:,16)<0.5;
All_Rivers_Qrange_1(QIR,:)=[];

QOR = All_Rivers_Qrange_2(:,16)>0.5;
All_Rivers_Qrange_2(QOR,:)=[];

% % Create Plot of Qc within range
% xline_r = [0 200];
% yline_r = [0.6 0.6];
% figure;
% set(gca, 'FontSize', 15)
% line(xline_r, yline_r, 'LineWidth', 1, 'Color', 'k')
% hold on
% scatter(All_Rivers_Qrange_1(:,15),All_Rivers_Qrange_1(:,1), 75);
% hold on
% scatter(All_Rivers_Qrange_2(:,15),All_Rivers_Qrange_2(:,1), 75, 'x');
% set(gca, 'XTickLabel', {});
% set(gca, 'XTick', [])
% %legend('Qmin < Qc < Qmax', ' Qc > Qmax or Qc < Qmin');

%---------------------------------------------------

% Calculate vc*dc*wc
% vcdcwc = vc_free.*dc_free.*wc_free;
% logvcdcwc = log10(vcdcwc);
% 
% logavg_discharge_r=log10(avg_discharge_r);
% logQc_r=log10(Qc_r);
% logavg_depth_r=log10(avg_depth_r);
% logavg_velocity_r=log10(avg_velocity_r);
% logavg_width_r=log10(avg_width_r);
% logdc_free_r=log10(dc_free_r);
% logvc_free_r=log10(vc_free_r);
% logwc_free_r=log10(wc_free_r);

% Plot vcdcwc vs. Qc
% Create 1:1 line
% xline = [-100 100];
% yline = [-100 100];
% 
% figure;
% scatter(logvcdcwc(:,1),logQc,75);
% xlim([-25 50])
% ylim([-25 50])
% set(gca, 'FontSize', 15)
% %xlabel('log $v_{c}d_{c}w_{c}$', 'interpreter', 'Latex','FontSize', 18)
% %ylabel('log $Q_{c}$', 'interpreter', 'Latex','FontSize', 18)
% hold on
% line(xline, yline, 'LineWidth', 1, 'Color', 'k');

% Boxplot of congruent parameters for strong AMHG rivers
% figure;
% subplot(2,2,1);
% boxplot([logavg_discharge_r logQc_r]);
% set(gca, 'FontSize', 15);
% set(gca, 'XTickLabel', {});
% %set(gca, 'XTick', [])
% hold on
% 
% 
% subplot(2,2,2);
% boxplot([logavg_depth_r logdc_free_r]);
% set(gca, 'FontSize', 15);
% set(gca, 'XTickLabel', {});
% %set(gca, 'XTick', [])

% 
% subplot(2,2,3)
% boxplot([logavg_velocity_r logvc_free_r]);
% set(gca, 'FontSize', 15);
% set(gca, 'XTickLabel', {});
% ylim([-9 10])
% 
% subplot(2,2,4)
% boxplot([logavg_width_r logwc_free_r]);
% set(gca, 'FontSize', 15)
% set(gca, 'XTickLabel', {});
% ylim([-10 8])
% 
% 

% figure;
% subplot(1,2,1);
% boxplot(logavg_width_r);
% set(gca, 'FontSize', 15)
% set(gca, 'XTickLabel', {});
% hold on
% 
% subplot(1,2,2);
% boxplot(logwc_free_r);
% ylim([-9 8]);
% set(gca, 'FontSize', 15)
% set(gca, 'XTickLabel', {});






% 
% 
% %---------------------------------------------------
% % Analyze the data
% 
% % % Plot rsq vs. log(Qc)
% % figure;
% % title('R Squared vs. Qc');
% % xlabel('R Squared');
% % ylabel('Log Qc');
% % axis([0.6 1 0 5]);
% % hold on
% % scatter(rsq_depthAMHG_free, logQc);
% % 


% %Plot Qc vs average discharge
% figure;
% subplot(1,2,1);
% scatter(logavg_discharge_r, logQc_r, 75);
% hold on
% set(gca, 'FontSize', 15);
% %xlabel('log $\bar{Q}$', 'interpreter', 'Latex', 'FontSize', 18);
% %ylabel('log $Q_c$', 'interpreter', 'Latex', 'FontSize', 18);
% axis([0 5 0 5]);
% 
% % dc vs average depth
% subplot(1,2,2);
% set(gca, 'FontSize', 15)
% %xlabel('log $\bar{d}$', 'interpreter', 'Latex', 'FontSize', 18);
% %ylabel('log $d_c$', 'interpreter', 'Latex', 'FontSize', 18);
% %axis([0 5 0 5]);
% hold on
% scatter(logavg_depth_r, logdc_free_r, 75);
% 
% % wc vs average width
% subplot(2,2,3)
% xlabel('log $\bar{w}$', 'interpreter', 'Latex', 'FontSize', 18);
% ylabel('log $w_c$', 'interpreter', 'Latex', 'FontSize', 18);
% %axis([0 5 0 5]);
% hold on
% scatter(logavg_width, logwc);
% 
% % vc vs average velocity
% subplot(2,2,4)
% xlabel('log $\bar{v}$', 'interpreter', 'Latex', 'FontSize', 18);
% ylabel('log $v_c$', 'interpreter', 'Latex', 'FontSize', 18);
% %axis([0 5 0 5]);
% hold on
% scatter(logavg_velocity, logvc);

% % 
% % % Plot rsq vs. average discharge
% % figure;
% % title('R Squared vs Average Discharge');
% % xlabel('R squared');
% % ylabel('Log Average Discharge');
% % axis([0.6 1 0 5]);
% % hold on
% % scatter(rsq_depthAMHG_free, logavg_discharge);
% % 
% % % Plot rsq vs. Froude number
% % figure;
% % title('Rsq vs Froude Number');
% % xlabel('R squared');
% % ylabel('Congruent Froude Number');
% % axis([0.6 1 0 3]);
% % hold on
% % scatter(rsq_depthAMHG_free, Fr_c_free);
% % 
% % % Plot average discharge vs. Froude number
% % figure;
% % title('Average Discharge vs Froude Number');
% % xlabel('Log Average Discharge');
% % ylabel('Congruent Froude Number');
% % axis([0 5 0 2]);
% % hold on
% % scatter(logavg_discharge, Fr_c_free);
% % 
% % 
% 
% if plots == 1
% 
%     % Plots vs Qc
%     figure;
%     subplot(3,3,1)
%     scatter(logavg_discharge,logQc,'b')
%     hold on
%     title('Qc vs. Q avg')
%     xlabel('log Q avg')
%     ylabel('log Qc')
%     
%     subplot(3,3,2)
%     scatter(logFr_c,logQc,'b')
%     hold on
%     title('Qc vs. log Fr_c')
%     xlabel('Fr_c')
%     ylabel('log Qc')
%     
%     subplot(3,3,3)
%     scatter(logwc, logQc,'b')
%     hold on
%     title('Qc vs. log wc')
%     xlabel('wc')
%     ylabel('log Qc')
%     
%     subplot(3,3,4)
%     scatter(dc_free, logQc,'b')
%     hold on
%     title('Log Qc vs. dc')
%     xlabel('dc')
%     ylabel('log Qc')
%     
%     subplot(3,3,5)
%     scatter(logvc, logQc,'b')
%     hold on
%     title('Qc vs. log vc')
%     xlabel('vc')
%     ylabel('log Qc')
%     
%     subplot(3,3,6)
%     scatter(avg_width, logQc,'b')
%     hold on
%     title('Qc vs. Average width')
%     xlabel('Average width')
%     ylabel('log Qc')
%     
%     subplot(3,3,7)
%     scatter(avg_depth, logQc,'b')
%     hold on
%     title('Qc vs. Average depth')
%     xlabel('Average depth')
%     ylabel('log Qc')
%     
%     subplot(3,3,8)
%     scatter(avg_velocity, logQc,'b')
%     hold on
%     title('Qc vs. Average Velocity')
%     xlabel('Average velocity')
%     ylabel('log Qc')
%     
%     subplot(3,3,9)
%     scatter(w.d, logQc,'b')
%     hold on
%     title('Qc vs. w/d')
%     xlabel('Average width / Average depth')
%     ylabel('log Qc')
%     
%     
%     % Plots vs. R squared
%     figure;
%     subplot(3,3,1)
%     scatter(logavg_discharge,rsq_depthAMHG_free,'b')
%     hold on
%     title('R squared vs. Q avg')
%     xlabel('log Q avg')
%     ylabel('R squared')
%     
%     subplot(3,3,2)
%     scatter(logFr_c,rsq_depthAMHG_free,'b')
%     hold on
%     title('R squared vs. Fr_c')
%     xlabel('Log Fr_c')
%     ylabel('R squared')
%     
%     subplot(3,3,3)
%     scatter(logwc, rsq_depthAMHG_free,'b')
%     hold on
%     title('R squared vs. wc')
%     xlabel('Log wc')
%     ylabel('R squared')
%     
%     subplot(3,3,4)
%     scatter(dc_free, rsq_depthAMHG_free,'b')
%     hold on
%     title('R squared vs. dc')
%     xlabel('dc')
%     ylabel('R squared')
%     
%     subplot(3,3,5)
%     scatter(logvc, rsq_depthAMHG_free,'b')
%     hold on
%     title('R squared vs. vc')
%     xlabel('Log vc')
%     ylabel('R squared')
%     
%     subplot(3,3,6)
%     scatter(avg_width, rsq_depthAMHG_free,'b')
%     hold on
%     title('R squared vs. Average width')
%     xlabel('Average width')
%     ylabel('R squared')
%     
%     subplot(3,3,7)
%     scatter(avg_depth, rsq_depthAMHG_free,'b')
%     hold on
%     title('R squared vs. Average depth')
%     xlabel('Average depth')
%     ylabel('R squared')
%     
%     subplot(3,3,8)
%     scatter(avg_velocity, rsq_depthAMHG_free,'b')
%     hold on
%     title('R squared vs. Average Velocity')
%     xlabel('Average velocity')
%     ylabel('R squared')
%     
%     subplot(3,3,9)
%     scatter(w.d, rsq_depthAMHG_free,'b')
%     hold on
%     title('R squared vs. w/d')
%     xlabel('Average width / Average depth')
%     ylabel('R squared')
%     
%     
%     % Plots vs. dc
%     figure;
%     subplot(3,3,1)
%     scatter(logavg_discharge,logdc,'b')
%     hold on
%     title('log dc vs. Q avg')
%     xlabel('log Q avg')
%     ylabel('log dc')
%     
%     subplot(3,3,2)
%     scatter(logFr_c,logdc,'b')
%     hold on
%     title('log dc vs. Fr_c')
%     xlabel('Log Fr_c')
%     ylabel('log dc')
%     
%     subplot(3,3,3)
%     scatter(logwc, logdc,'b')
%     hold on
%     title('log dc vs. wc')
%     xlabel('Log wc')
%     ylabel('log dc')
%     
%     subplot(3,3,4)
%     scatter(dc_free, logdc,'b')
%     hold on
%     title('log dc vs. dc')
%     xlabel('dc')
%     ylabel('log dc')
%     
%     subplot(3,3,5)
%     scatter(logvc, logdc,'b')
%     hold on
%     title('log dc vs. vc')
%     xlabel('Log vc')
%     ylabel('log dc')
%     
%     subplot(3,3,6)
%     scatter(avg_width, logdc,'b')
%     hold on
%     title('log dc vs. Average width')
%     xlabel('Average width')
%     ylabel('log dc')
%     
%     subplot(3,3,7)
%     scatter(avg_depth, logdc,'b')
%     hold on
%     title('log dc vs. Average depth')
%     xlabel('Average depth')
%     ylabel('log dc')
%     
%     subplot(3,3,8)
%     scatter(avg_velocity, logdc,'b')
%     hold on
%     title('log dc vs. Average Velocity')
%     xlabel('Average velocity')
%     ylabel('log dc')
%     
%     subplot(3,3,9)
%     scatter(w.d, logdc,'b')
%     hold on
%     title('log dc vs. w/d')
%     xlabel('Average width / Average depth')
%     ylabel('log dc')
%     
%     
%     % Log dc / Log Qc
%     figure;
%     dc.Qc = logdc./logQc;
%     plot(dc.Qc, 'bo');
%     title('Log dc / Log Qc')
% 
% end
