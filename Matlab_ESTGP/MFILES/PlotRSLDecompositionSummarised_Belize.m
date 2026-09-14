function[]=PlotRSLDecompositionSummarised(filepath,toplot,group,f2s,sd2s,V2s,testlocs,unique_index,imageres)

    % DIR = '/Users/fangyi/Desktop/Codes/Erica STEHM/Belize_EST_GP/RSLDecompositionPlots';\
    DIR=filepath;
    if ~exist(DIR,'dir')
        mkdir(DIR);
    end
    cd(DIR);

    for k=1:length(toplot)
    
        index = toplot(k);
    
        %%%%%% Summary plots for chosen model %%%%%%
        
        % For the chosen model, overlay the RSL predictions (total, global, local)
        % from various sites
        
        colourmap = brewermap(length(testlocs.sites)-1, 'Dark2');  % Colormap for plotting
        stepsize=length(testlocs.X)/length(testlocs.names); % no. of rows per site to extract each site's data
        allsites_index = length(testlocs.X)-stepsize+1:length(testlocs.X); % last set of data are for all sites combined
        legend_all=[]; % initiate handles to store legends
        % initiate handles to store rates and set x axis limits for rate plots
        df2s_all=[];
        dsd2s_all=[];
        % define y-axis limits based on total RSL signal
        ymin=min(f2s{unique_index(index)}{1}/1000-2*sd2s{unique_index(index)}{1}/1000)-2;
        ymax=max(f2s{unique_index(index)}{1}/1000+2*sd2s{unique_index(index)}{1}/1000)+2;
        
        figure
        % plot global RSL
        subplot(2,4,2)
        t=1950-testlocs.X(allsites_index,3);
        minage=min(t);
        maxage=max(t);
        RSLmean=f2s{unique_index(index)}{2}(1:stepsize)/1000; % f2s{m}{1}: the {1} extracts the full posterior output without white
        % noise; use {2} for global component only {3} for regional non-linear component only and {4} for local component only
        RSLsd = sd2s{unique_index(index)}{2}(1:stepsize)/1000;
        h=plot(t,RSLmean,'Color',[0.5 0.5 0.5]); % RSL mean line
        legend_all = [legend_all; h];
        hold on
        fill([t; flipud(t)], [RSLmean+RSLsd; flipud(RSLmean-RSLsd)],[0.5 0.5 0.5],'FaceAlpha',0.1,'EdgeColor',[0.5 0.5 0.5],'EdgeAlpha',0.5); % RSL 1sd
        fill([t; flipud(t)], [RSLmean+2*RSLsd; flipud(RSLmean-2*RSLsd)],[0.5 0.5 0.5],'FaceAlpha',0.1,'EdgeColor',[0.5 0.5 0.5],'EdgeAlpha',0.5); % RSL 2sd
        xlim([minage maxage]); % set x axis limits 
        ylim([ymin ymax]); % set y axis limits
        ylabel('RSL (m)');
        xlabel('Age (BP)')
        set(gca, 'XDir', 'reverse');% Reverse x axes
        title('Global RSL')
        hold on
        
        % Calculate RSL rates
        difftimestep=testlocs.X(2,3)-testlocs.X(1,3); % Temporal spacing (in years) between time steps for computing and plotting rate estimates.
        % Find time-paired rows separated by exactly 'difftimestep' units
        Mdiff = bsxfun(@eq,testlocs.X(allsites_index,3),testlocs.X(allsites_index,3)')-bsxfun(@eq,testlocs.X(allsites_index,3),testlocs.X(allsites_index,3)'+ difftimestep);
        % Remove rows without any matching time-difference pairs (i.e., only
        % want rows with difftimestep
        sub=find(sum(Mdiff,2)==0);  
        Mdiff=Mdiff(sub,:);
        % Compute average time difference (denominator for rate calculation)
        difftimes=1950-bsxfun(@rdivide,abs(Mdiff)*testlocs.X(allsites_index,3),sum(abs(Mdiff),2));
        % Normalize each row to compute finite differences (df/dt)
        Mdiff=bsxfun(@rdivide,Mdiff,Mdiff*testlocs.X(allsites_index,3));
        
        % Plot global RSL rate
        df2s=Mdiff*f2s{unique_index(index)}{2}(allsites_index);
        dV2s=Mdiff*V2s{unique_index(index)}{2}(allsites_index,allsites_index)*Mdiff';
        dsd2s=sqrt(diag(dV2s(:,:)));
        subplot(2,4,6);
        plot(difftimes,df2s,'Color',[0.5 0.5 0.5]); % RSL rate mean
        hold on
        fill([difftimes;flipud(difftimes)],[df2s+dsd2s; flipud(df2s-dsd2s)],[0.5 0.5 0.5],'FaceAlpha',0.1,'EdgeColor',[0.5 0.5 0.5],'EdgeAlpha',0.5); % RSL rate 1sd
        fill([difftimes;flipud(difftimes)],[df2s+2*dsd2s; flipud(df2s-2*dsd2s)],[0.5 0.5 0.5],'FaceAlpha',0.1,'EdgeColor',[0.5 0.5 0.5],'EdgeAlpha',0.5) % RSL rate 2sd
        xlim([minage maxage]); % set x axis limits but let y axis float
        ylabel({'RSL Rate', '(mm/yr = m/ka)'});
        xlabel('Age (BP)')
        set(gca, 'XDir', 'reverse');% Reverse x axes
        title('Global RSL rate');
        
        % Append rates to calculate limits
        df2s_all = [df2s_all; df2s];
        dsd2s_all = [dsd2s_all; dsd2s];
        
        for n=1:length(testlocs.sites)-1 % for each site (exclude last row as that is average position for whole Twin Cays)
            
            siteindex = stepsize*(n-1)+1:stepsize*n; % rows for that site
            col=colourmap(n,:);
            t=1950-testlocs.X(siteindex,3); % time for predictions in yrs BP
            minage=min(t);
            maxage=max(t);
        
            % Plot site-specific RSL
            subplot(2,4,1)
            RSLmean=f2s{unique_index(index)}{1}(siteindex)/1000; % f2s{m}{1}: the {1} extracts the full posterior output without white
            % noise; use {2} for global component only {3} for regional non-linear component only and {4} for local component only
            RSLsd = sd2s{unique_index(index)}{1}(siteindex)/1000;
            h=plot(t,RSLmean,'Color',col,'DisplayName',sprintf('Site %d',n)); % RSL mean line
            legend_all = [legend_all; h];
            hold on
            fill([t; flipud(t)], [RSLmean+RSLsd; flipud(RSLmean-RSLsd)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 1sd
            fill([t; flipud(t)], [RSLmean+2*RSLsd; flipud(RSLmean-2*RSLsd)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 2sd
            xlim([minage maxage]); % set x axis limits 
            ylim([ymin ymax]); % set y axis limits
            ylabel('RSL (m)');
            xlabel('Age (BP)')
            set(gca, 'XDir', 'reverse');% Reverse x axes
            title('RSL')
            hold on
        
            % Plot regional non-linear RSL component for that site 
            subplot(2,4,3)
            localmean=f2s{unique_index(index)}{3}(siteindex)/1000; % f2s{m}{1}: the {1} extracts the full posterior output without white
            % noise; use {2} for global component only {3} for regional non-linear component only and {4} for local component only
            localsd = sd2s{unique_index(index)}{3}(siteindex)/1000;
            plot(t,localmean,'Color',col); % RSL mean line
            hold on
            fill([t; flipud(t)], [localmean+localsd; flipud(localmean-localsd)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 1sd
            fill([t; flipud(t)], [localmean+2*localsd; flipud(localmean-2*localsd)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 2sd
            xlim([minage maxage]); % set x axis limits 
            ylim([ymin ymax]); % set y axis limits
            ylabel('RSL (m)');
            xlabel('Age (BP)')
            set(gca, 'XDir', 'reverse');% Reverse x axes
            title('Regional non-linear RSL')
            hold on
        
            % Plot local RSL component for that site 
            subplot(2,4,4)
            localmean=f2s{unique_index(index)}{4}(siteindex)/1000; % f2s{m}{1}: the {1} extracts the full posterior output without white
            % noise; use {2} for global component only {3} for regional non-linear component only and {4} for local component only
            localsd = sd2s{unique_index(index)}{4}(siteindex)/1000;
            plot(t,localmean,'Color',col); % RSL mean line
            hold on
            fill([t; flipud(t)], [localmean+localsd; flipud(localmean-localsd)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 1sd
            fill([t; flipud(t)], [localmean+2*localsd; flipud(localmean-2*localsd)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 2sd
            xlim([minage maxage]); % set x axis limits 
            ylim([ymin ymax]); % set y axis limits
            ylabel('RSL (m)');
            xlabel('Age (BP)')
            set(gca, 'XDir', 'reverse');% Reverse x axes
            title('Local RSL')
            hold on
        
            % Plot site-specific RSL rate
            df2s=Mdiff*f2s{unique_index(index)}{1}(siteindex);
            dV2s=Mdiff*V2s{unique_index(index)}{1}(siteindex,siteindex)*Mdiff';
            dsd2s=sqrt(diag(dV2s(:,:)));
            subplot(2,4,5);
            plot(difftimes,df2s,'Color',col); % RSL rate mean
            hold on
            fill([difftimes;flipud(difftimes)],[df2s+dsd2s; flipud(df2s-dsd2s)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL rate 1sd
            fill([difftimes;flipud(difftimes)],[df2s+2*dsd2s; flipud(df2s-2*dsd2s)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5) % RSL rate 2sd
            xlim([minage maxage]); % set x axis limits but let y axis float
            ylabel({'RSL Rate', '(mm/yr = m/ka)'});
            xlabel('Age (BP)')
            set(gca, 'XDir', 'reverse');% Reverse x axes
            title('RSL rate');
        
             % Append rates to calculate limits
            df2s_all = [df2s_all; df2s];
            dsd2s_all = [dsd2s_all; dsd2s];
         
            % Plot regional non-linear rates
            df2s=Mdiff*f2s{unique_index(index)}{3}(siteindex);
            dV2s=Mdiff*V2s{unique_index(index)}{3}(siteindex,siteindex)*Mdiff';
            dsd2s=sqrt(diag(dV2s(:,:)));
            subplot(2,4,7);
            plot(difftimes,df2s,'Color',col); % RSL rate mean
            hold on
            fill([difftimes;flipud(difftimes)],[df2s+dsd2s; flipud(df2s-dsd2s)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL rate 1sd
            fill([difftimes;flipud(difftimes)],[df2s+2*dsd2s; flipud(df2s-2*dsd2s)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5) % RSL rate 2sd
            xlim([minage maxage]); % set x axis limits but let y axis float
            ylabel({'RSL Rate', '(mm/yr = m/ka)'});
            xlabel('Age (BP)')
            set(gca, 'XDir', 'reverse');% Reverse x axes
            title('Regional non-linear RSL rate');
        
             % Append rates to calculate limits
            df2s_all = [df2s_all; df2s];
            dsd2s_all = [dsd2s_all; dsd2s];
        
            % Plot local rates
            df2s=Mdiff*f2s{unique_index(index)}{4}(siteindex);
            dV2s=Mdiff*V2s{unique_index(index)}{4}(siteindex,siteindex)*Mdiff';
            dsd2s=sqrt(diag(dV2s(:,:)));
            subplot(2,4,8);
            plot(difftimes,df2s,'Color',col); % RSL rate mean
            hold on
            fill([difftimes;flipud(difftimes)],[df2s+dsd2s; flipud(df2s-dsd2s)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL rate 1sd
            fill([difftimes;flipud(difftimes)],[df2s+2*dsd2s; flipud(df2s-2*dsd2s)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5) % RSL rate 2sd
            xlim([minage maxage]); % set x axis limits but let y axis float
            ylabel({'RSL Rate', '(mm/yr = m/ka)'});
            xlabel('Age (BP)')
            set(gca, 'XDir', 'reverse');% Reverse x axes
            title('Local RSL rate');
        
             % Append rates to calculate limits
            df2s_all = [df2s_all; df2s];
            dsd2s_all = [dsd2s_all; dsd2s];
        
        end
        
        % set y axis limits for rate plots
        minrate = min(df2s_all-2*dsd2s_all);
        maxrate = max(df2s_all+2*dsd2s_all);
        buffer=(maxrate-minrate)*0.1;
        for sp = 5:8
            subplot(2, 4, sp);
            ylim([minrate-buffer maxrate+buffer]);
        end
        
        % Now add legend inside figure, below subplots
        lgd = legend(legend_all, 'Orientation', 'horizontal', 'NumColumns', ceil(length(legend_all)/2));
        
        % Adjust legend position (tweak these numbers as needed)
        lgd.Position = [0.3, 0, 0.4, 0.06]; % [left bottom width height]
            
        % save figure as pdf 
        label=sprintf('RSL Decomposition model %d (%s).pdf',index,group);
        fig = gcf;
        fig.Units = 'inches';         % Set units to inches
        % fig.Position = [1, 1, 10, 8]; % [left, bottom, width, height] — adjust width (e.g., 10 inches)
        % exportgraphics(gcf, label, 'ContentType', 'vector');

        if strcmpi(imageres, 'low')
        fig.Position = [1, 1, 10, 8]; % [left, bottom, width, height] — adjust width (e.g., 10 inches)
        exportgraphics(gcf, label, 'ContentType', 'image', 'Resolution', 300);
        else 
        fig.Position = [1, 1, 10, 6]; % [left, bottom, width, height] — adjust width (e.g., 10 inches)
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 7 6]; % [left, bottom, width, height] — adjust width (e.g., 10 inches)
        fig.PaperSize = [7 6];
        print(gcf, label, '-dpdf', '-painters', '-r300');
        end

    
    
    end

end