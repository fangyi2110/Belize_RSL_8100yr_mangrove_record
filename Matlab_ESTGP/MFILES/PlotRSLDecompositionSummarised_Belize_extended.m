function[]=PlotRSLDecompositionSummarised(filepath,datasets,toplot,group,f2s,sd2s,V2s,testlocs,unique_index,imageres,RSLmaxage)

    % DIR = '/Users/fangyi/Desktop/Codes/Erica STEHM/Belize_EST_GP/RSLDecompositionPlots';\
    DIR=filepath;
    if ~exist(DIR,'dir')
        mkdir(DIR);
    end
    cd(DIR);

    % extract data for plotting SLIPs and limiting data
    data=datasets{1,1};

    for k=1:length(toplot)
    
        index = toplot(k);
    
        %%%%%% Summary plots for chosen model %%%%%%
        
        % For the chosen model, overlay the RSL predictions (total, global, local)
        % from various sites
        
        colourmap = brewermap(length(testlocs.sites)-1, 'Dark2');  % Colormap for plotting

        allsites_index = find(testlocs.reg==99999); % find index for the average across all sites in Belize
        legend_all=[]; % initiate handles to store legends
        % initiate handles to store rates and set x axis limits for rate plots
        df2s_all=[];
        dsd2s_all=[];

        figure
        % plot global RSL
        subplot(2,4,2)
        t=1950-testlocs.X(allsites_index,3);
        RSLmean=f2s{unique_index(index)}{2}(allsites_index)/1000; % f2s{m}{1}: the {1} extracts the full posterior output without white
        % noise; use {2} for global component only {3} for regional non-linear component only and {4} for local component only
        RSLsd = sd2s{unique_index(index)}{2}(allsites_index)/1000;
        ageidx = find(t<=RSLmaxage); % find index where prediction age is younger than max age of Belize SLIPs
        % filter predictions to time span of Belize
        t_cut = t(ageidx); 
        RSLmean_cut = RSLmean(ageidx);
        RSLsd_cut = RSLsd(ageidx);
        h=plot(t_cut,RSLmean_cut,'Color',[0.5 0.5 0.5]); % RSL mean line
        legend_all = [legend_all; h];
        hold on
        fill([t_cut; flipud(t_cut)], [RSLmean_cut+RSLsd_cut; flipud(RSLmean_cut-RSLsd_cut)],[0.5 0.5 0.5],'FaceAlpha',0.1,'EdgeColor',[0.5 0.5 0.5],'EdgeAlpha',0.5); % RSL 1sd
        fill([t_cut; flipud(t_cut)], [RSLmean_cut+2*RSLsd_cut; flipud(RSLmean_cut-2*RSLsd_cut)],[0.5 0.5 0.5],'FaceAlpha',0.1,'EdgeColor',[0.5 0.5 0.5],'EdgeAlpha',0.5); % RSL 2sd
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
        ageidx = find(difftimes<=RSLmaxage); % find index where prediction age is younger than max age of Belize SLIPs
        % filter predictions to time span of Belize
        difftimes_cut = difftimes(ageidx); 
        df2s_cut = df2s(ageidx);
        dsd2s_cut = dsd2s(ageidx);
        subplot(2,4,6);
        plot(difftimes_cut,df2s_cut,'Color',[0.5 0.5 0.5]); % RSL rate mean
        hold on
        fill([difftimes_cut;flipud(difftimes_cut)],[df2s_cut+dsd2s_cut; flipud(df2s_cut-dsd2s_cut)],[0.5 0.5 0.5],'FaceAlpha',0.1,'EdgeColor',[0.5 0.5 0.5],'EdgeAlpha',0.5); % RSL rate 1sd
        fill([difftimes_cut;flipud(difftimes_cut)],[df2s_cut+2*dsd2s_cut; flipud(df2s_cut-2*dsd2s_cut)],[0.5 0.5 0.5],'FaceAlpha',0.1,'EdgeColor',[0.5 0.5 0.5],'EdgeAlpha',0.5) % RSL rate 2sd
        ylabel({'RSL Rate', '(mm/yr = m/ka)'});
        xlabel('Age (BP)')
        set(gca, 'XDir', 'reverse');% Reverse x axes
        title('Global RSL rate');
        
        % Append rates to calculate limits
        df2s_all = [df2s_all; df2s];
        dsd2s_all = [dsd2s_all; dsd2s];
        
        for n=1:length(testlocs.sites)-2 % loop for all local sites (remove last two rows as that is for Barbados and mean across Twin Cays)
            
            siteindex = find(testlocs.reg==testlocs.sites(n,1)); % rows for that site
            col=colourmap(n,:);
            t=1950-testlocs.X(siteindex,3); % time for predictions in yrs BP
        
            % Plot site-specific RSL
            subplot(2,4,1)
            RSLmean=f2s{unique_index(index)}{1}(siteindex)/1000; % f2s{m}{1}: the {1} extracts the full posterior output without white
            % noise; use {2} for global component only {3} for regional non-linear component only and {4} for local component only
            RSLsd = sd2s{unique_index(index)}{1}(siteindex)/1000;
            ageidx = find(t<=RSLmaxage); % find index where prediction age is younger than max age of Belize SLIPs
            % filter predictions to time span of Belize
            t_cut = t(ageidx); 
            RSLmean_cut = RSLmean(ageidx);
            RSLsd_cut = RSLsd(ageidx);
            % plot RSL data 
            for a=1:length(data.datid)
                hold on
                if data.limiting(a)==0 % plot SLIPs
                    xmin=min(1950-data.time1(a),1950-data.time2(a));
                    ymin=(data.Y(a)-2*data.dY(a))/1000;
                    width=data.time2(a)-data.time1(a);
                    height=(4*data.dY(a))/1000;
                    rectangle('Position',[xmin,ymin,width,height]);
                elseif data.limiting(a)==-1 % plot marine limiting (navy)
                    plot([1 1]*(1950-data.meantime(a)),1e-3*[data.Y(a)-2*data.dY(a)-500,data.Y(a)-2*data.dY(a)],'Color',[.1,.35,.7]); 
                    plot([1950-data.time1(a) 1950-data.time2(a)],[1 1]*1e-3*(data.Y(a)-2*data.dY(a)),'Color',[.1,.35,.7]);
                elseif data.limiting(a)==1 % plot terrestrial limiting (sea green)
                    plot([1 1]*(1950-data.meantime(a)),1e-3*[data.Y(a)+2*data.dY(a),data.Y(a)+2*data.dY(a)+500],'Color',[.1,.35,.7]); 
                    plot([1950-data.time1(a) 1950-data.time2(a)],[1 1]*1e-3*(data.Y(a)+2*data.dY(a)),'Color',[.1,.35,.7]);
                end
            end
            hold on
            h=plot(t_cut,RSLmean_cut,'Color',col,'DisplayName',sprintf('Site %d',n)); % RSL mean line
            legend_all = [legend_all; h];
            hold on
            fill([t_cut; flipud(t_cut)], [RSLmean_cut+RSLsd_cut; flipud(RSLmean_cut-RSLsd_cut)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 1sd
            fill([t_cut; flipud(t_cut)], [RSLmean_cut+2*RSLsd_cut; flipud(RSLmean_cut-2*RSLsd_cut)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 2sd
            ylabel('RSL (m)');
            xlabel('Age (BP)')
            % define x-axis limits to maximum age of Belize data, including
            % limiting data
            BelizeID = find(data.datid~=1029000);
            xmin = 1950-max(data.time2(BelizeID));
            xmax = 1950-min(data.time1(BelizeID));
            xlim([xmin xmax]);
            % set y axis limits for all RSL subplots based on full RSL
            minRSL = min(RSLmean_cut-2*RSLsd_cut);
            maxRSL = max(RSLmean_cut+2*RSLsd_cut);
            buffer=(maxRSL-minRSL)*0.05;
            ylim([minRSL-buffer maxRSL+buffer]);
            set(gca, 'XDir', 'reverse');% Reverse x axes
            title('RSL')
            hold on
            
            % standardise x axis limits for all subsequent RSL subplots to be standa
            minage = min(t_cut);
            maxage= max(t_cut);
        
            % Plot regional non-linear RSL component for that site 
            subplot(2,4,3)
            localmean=f2s{unique_index(index)}{3}(siteindex)/1000; % f2s{m}{1}: the {1} extracts the full posterior output without white
            % noise; use {2} for global component only {3} for regional non-linear component only and {4} for local component only
            localsd = sd2s{unique_index(index)}{3}(siteindex)/1000;
            ageidx = find(t<=RSLmaxage); % find index where prediction age is younger than max age of Belize SLIPs
            % filter predictions to time span of Belize
            t_cut = t(ageidx); 
            localmean_cut = localmean(ageidx);
            localsd_cut = localsd(ageidx);
            plot(t,localmean,'Color',col); % RSL mean line
            hold on
            fill([t_cut; flipud(t_cut)], [localmean_cut+localsd_cut; flipud(localmean_cut-localsd_cut)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 1sd
            fill([t_cut; flipud(t_cut)], [localmean_cut+2*localsd_cut; flipud(localmean_cut-2*localsd_cut)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 2sd
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
            ageidx = find(t<=RSLmaxage); % find index where prediction age is younger than max age of Belize SLIPs
            % filter predictions to time span of Belize
            t_cut = t(ageidx); 
            localmean_cut = localmean(ageidx);
            localsd_cut = localsd(ageidx);
            plot(t,localmean,'Color',col); % RSL mean line
            hold on
            fill([t_cut; flipud(t_cut)], [localmean_cut+localsd_cut; flipud(localmean_cut-localsd_cut)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 1sd
            fill([t_cut; flipud(t_cut)], [localmean_cut+2*localsd_cut; flipud(localmean_cut-2*localsd_cut)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 2sd
            ylabel('RSL (m)');
            xlabel('Age (BP)')
            set(gca, 'XDir', 'reverse');% Reverse x axes
            title('Local RSL')
            hold on
        
            % Plot site-specific RSL rate
            % Calculate RSL rates
            difftimestep=testlocs.X(2,3)-testlocs.X(1,3); % Temporal spacing (in years) between time steps for computing and plotting rate estimates.
            % Find time-paired rows separated by exactly 'difftimestep' units
            Mdiff = bsxfun(@eq,testlocs.X(siteindex,3),testlocs.X(siteindex,3)')-bsxfun(@eq,testlocs.X(siteindex,3),testlocs.X(siteindex,3)'+ difftimestep);
            % Remove rows without any matching time-difference pairs (i.e., only
            % want rows with difftimestep
            sub=find(sum(Mdiff,2)==0);  
            Mdiff=Mdiff(sub,:);
            % Compute average time difference (denominator for rate calculation)
            difftimes=1950-bsxfun(@rdivide,abs(Mdiff)*testlocs.X(siteindex,3),sum(abs(Mdiff),2));
            % Normalize each row to compute finite differences (df/dt)
            Mdiff=bsxfun(@rdivide,Mdiff,Mdiff*testlocs.X(siteindex,3));
            df2s=Mdiff*f2s{unique_index(index)}{1}(siteindex);
            dV2s=Mdiff*V2s{unique_index(index)}{1}(siteindex,siteindex)*Mdiff';
            dsd2s=sqrt(diag(dV2s(:,:)));
            subplot(2,4,5);
            plot(difftimes,df2s,'Color',col); % RSL rate mean
            hold on
            fill([difftimes;flipud(difftimes)],[df2s+dsd2s; flipud(df2s-dsd2s)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL rate 1sd
            fill([difftimes;flipud(difftimes)],[df2s+2*dsd2s; flipud(df2s-2*dsd2s)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5) % RSL rate 2sd
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
            ylabel({'RSL Rate', '(mm/yr = m/ka)'});
            xlabel('Age (BP)')
            set(gca, 'XDir', 'reverse');% Reverse x axes
            title('Local RSL rate');
        
             % Append rates to calculate limits
            df2s_all = [df2s_all; df2s];
            dsd2s_all = [dsd2s_all; dsd2s];
        
        end

         % set y and x axis limits for RSL plots
        for sp = 2:4
            subplot(2, 4, sp);
            ylim([minRSL-buffer maxRSL+buffer]);
            xlim([minage maxage]);
        end
        
        % set y and x axis limits for rate plots
        minrate = min(df2s_all-2*dsd2s_all);
        maxrate = max(df2s_all+2*dsd2s_all);
        buffer=(maxrate-minrate)*0.1;
        for sp = 5:8
            subplot(2, 4, sp);
            ylim([minrate-buffer maxrate+buffer]);
            xlim([minage maxage]);
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