function[]=PlotRSLDecompositionSummarised(runs_to_plot,RSLmaxage, datasets,f2s,sd2s,V2s,testlocs,legendrows,modelnumbers,subplotcolno,filename,imageres)

% select model to plot
n_plot = length(runs_to_plot); % number of models being plotted
colourmap = brewermap(n_plot, 'Dark2');  % Colormap for plotting

% initiate handles to store legends
legend_all=[];

% extract data for plotting SLIPs and limiting data
data=datasets{1,1};

% determine panels for subplots
subplotrowno = ceil((2*(length(testlocs.names)-1)+1)/subplotcolno);

%% Plot RSL 

figure

for k = 1:n_plot

    m=runs_to_plot(k);
    modelnumber = modelnumbers(k);

    if isempty(f2s{m}) || isempty(sd2s{m}) || isempty(V2s{m})
        continue; % Skip if model data is missing
    end

    % color map for plotting
    col = colourmap(k, :);

    % row index for all sites ('Holo-TwinCays_all' in testlocs.names)
    index = find(testlocs.reg==99999);
    
    %% Plot RSL

    % Plot Full RSL across all sites 
    subplot(subplotrowno,subplotcolno,1);
    t=1950-testlocs.X(index,3);
    RSLmean=f2s{m}{1}(index)/1000; % f2s{m}{1}: the {1} extracts the full posterior output without white
    % noise; use {2} for global component only {3} for regional non-linear component only and {4} for local component only
    RSLsd = sd2s{m}{1}(index)/1000;
    for n=1:length(data.datid)
            hold on
            if data.limiting(n)==0 % plot SLIPs
                xmin=min(1950-data.time1(n),1950-data.time2(n));
                ymin=(data.Y(n)-2*data.dY(n))/1000;
                width=data.time2(n)-data.time1(n);
                height=(4*data.dY(n))/1000;
                rectangle('Position',[xmin,ymin,width,height]);
            elseif data.limiting(n)==-1 % plot marine limiting (navy)
                plot([1 1]*(1950-data.meantime(n)),1e-3*[data.Y(n)-2*data.dY(n)-500,data.Y(n)-2*data.dY(n)],'Color',[.1,.35,.7]); 
                plot([1950-data.time1(n) 1950-data.time2(n)],[1 1]*1e-3*(data.Y(n)-2*data.dY(n)),'Color',[.1,.35,.7]);
            elseif data.limiting(n)==1 % plot terrestrial limiting (sea green)
                plot([1 1]*(1950-data.meantime(n)),1e-3*[data.Y(n)+2*data.dY(n),data.Y(n)+2*data.dY(n)+500],'Color',[.1,.35,.7]); 
                plot([1950-data.time1(n) 1950-data.time2(n)],[1 1]*1e-3*(data.Y(n)+2*data.dY(n)),'Color',[.1,.35,.7]);
            end
    end
    hold on
    h=plot(t,RSLmean,'Color',col,'DisplayName',sprintf('Model %d',modelnumber)); % RSL mean line
    legend_all = [legend_all; h];
    fill([t; flipud(t)], [RSLmean+RSLsd; flipud(RSLmean-RSLsd)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 1sd
    fill([t; flipud(t)], [RSLmean+2*RSLsd; flipud(RSLmean-2*RSLsd)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 2sd
    ylabel('RSL (m)','Color','k');
    xlabel('Age (BP)')
    set(gca, 'XDir', 'reverse');% Reverse x axes
    title('Full RSL')

    % Plot global RSL component
    subplot(subplotrowno,subplotcolno,2);
    t=1950-testlocs.X(index,3);
    RSLmean=f2s{m}{2}(index)/1000; % f2s{m}{1}: the {1} extracts the full posterior output without white
    % noise; use {2} for global component only {3} for regional non-linear component only and {4} for local component only
    RSLsd = sd2s{m}{2}(index)/1000;
    ageidx = find(t<=RSLmaxage); % find index where prediction age is younger than max age of Belize SLIPs
    % filter predictions to time span of Belize
    t_cut = t(ageidx); 
    RSLmean_cut = RSLmean(ageidx);
    RSLsd_cut = RSLsd(ageidx);
    % plot
    plot(t_cut,RSLmean_cut,'Color',col); % RSL mean line
    hold on
    fill([t_cut; flipud(t_cut)], [RSLmean_cut+RSLsd_cut; flipud(RSLmean_cut-RSLsd_cut)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 1sd
    fill([t_cut; flipud(t_cut)], [RSLmean_cut+2*RSLsd_cut; flipud(RSLmean_cut-2*RSLsd_cut)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 2sd
    ylabel('RSL (m)','Color','k');
    xlabel('Age (BP)')
    set(gca, 'XDir', 'reverse');% Reverse x axes
    title('Global')

    % set y axis limits for global, local and regional subplots based on global component
    minRSL = min(RSLmean_cut-2*RSLsd_cut);
    maxRSL = max(RSLmean_cut+2*RSLsd_cut);
    buffer=(maxRSL-minRSL)*0.2;
    minage = min(t_cut);
    maxage= max(t_cut);
    
     % Plot regional components for each site
    for i=1:length(testlocs.names)-2 % loop for all local sites (remove last two rows as that is for Barbados and mean across Twin Cays)
        
        siteindex = find(testlocs.reg==testlocs.sites(i,1));
        
        subplot(subplotrowno,subplotcolno,2+i); % skip first two panels 
        t=1950-testlocs.X(siteindex,3);
        RSLmean=f2s{m}{3}(siteindex)/1000; % f2s{m}{1}: the {1} extracts the full posterior output without white
        % noise; use {2} for global component only {3} for regional non-linear component only and {4} for local component only
        RSLsd = sd2s{m}{3}(siteindex)/1000;
        ageidx = find(t<=RSLmaxage); % find index where prediction age is younger than max age of Belize SLIPs
        % filter predictions to time span of Belize
        t_cut = t(ageidx); 
        RSLmean_cut = RSLmean(ageidx);
        RSLsd_cut = RSLsd(ageidx);
        % plot
        plot(t_cut,RSLmean_cut,'Color',col); % RSL mean line
        hold on
        fill([t_cut; flipud(t_cut)], [RSLmean_cut+RSLsd_cut; flipud(RSLmean_cut-RSLsd_cut)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 1sd
        fill([t_cut; flipud(t_cut)], [RSLmean_cut+2*RSLsd_cut; flipud(RSLmean_cut-2*RSLsd_cut)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 2sd
        ylabel('RSL (m)','Color','k');
        xlabel('Age (BP)')
        set(gca, 'XDir', 'reverse');% Reverse x axes
        title(sprintf('Regional (Site %d)',i));
    end
    
    % Plot local components for each site
    for i=1:length(testlocs.names)-2 % loop for all local sites (remove last two rows as that is for Barbados and mean across Twin Cays)
        
        siteindex = find(testlocs.reg==testlocs.sites(i,1));
        
        subplot(subplotrowno,subplotcolno,2+length(testlocs.names)-2+i); % skip first two panels and regional panels
        t=1950-testlocs.X(siteindex,3);
        RSLmean=f2s{m}{4}(siteindex)/1000; % f2s{m}{1}: the {1} extracts the full posterior output without white
        % noise; use {2} for global component only {3} for regional non-linear component only and {4} for local component only
        RSLsd = sd2s{m}{4}(siteindex)/1000;
        ageidx = find(t<=RSLmaxage); % find index where prediction age is younger than max age of Belize SLIPs
        % filter predictions to time span of Belize
        t_cut = t(ageidx); 
        RSLmean_cut = RSLmean(ageidx);
        RSLsd_cut = RSLsd(ageidx);
        % plot
        plot(t_cut,RSLmean_cut,'Color',col); % RSL mean line
        hold on
        fill([t_cut; flipud(t_cut)], [RSLmean_cut+RSLsd_cut; flipud(RSLmean_cut-RSLsd_cut)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 1sd
        fill([t_cut; flipud(t_cut)], [RSLmean_cut+2*RSLsd_cut; flipud(RSLmean_cut-2*RSLsd_cut)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 2sd
        ylabel('RSL (m)','Color','k');
        xlabel('Age (BP)')
        set(gca, 'XDir', 'reverse');% Reverse x axes
        title(sprintf('Local (Site %d)',i));
    end
end

% Loop through all subplots for the global, local and regional component and apply the same y-limits and x-limits
numplots = (subplotrowno * subplotcolno); % Total subplot slots
for sp = 2:numplots
    subplot(subplotrowno, subplotcolno, sp);
    ylim([minRSL-buffer maxRSL+buffer]);
    xlim([minage maxage]);
end

% % Now add legend inside figure, below subplots
% lgd = legend(legend_all, 'Orientation', 'horizontal', 'NumColumns', ceil(length(legend_all)/legendrows));
% 
% % Adjust legend position (tweak these numbers as needed)
% lgd.Position = [0.3, 0.0, 0.4, 0.06]; % [left bottom width height]

% save figure as pdf 
fig = gcf;
fig.Units = 'inches';         % Set units to inches
if strcmpi(imageres, 'low')
    fig.Position = [1, 1, 15, 8]; % [left, bottom, width, height] — adjust width (e.g., 10 inches)
    exportgraphics(gcf, strcat(filename,'_RSL.pdf'), 'ContentType', 'image', 'Resolution', 300);
else 
    fig.Position = [1, 1, 15, 6]; % [left, bottom, width, height] — adjust width (e.g., 10 inches)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 10 6]; % [left, bottom, width, height] — adjust width (e.g., 10 inches)
    fig.PaperSize = [10 6];
    print(gcf, strcat(filename,'_RSL.pdf'), '-dpdf', '-painters', '-r300');
end

%% Plot RSL rates

% initiate handles to store legends
legend_all=[];

% initiate handles to store rates
df2s_all=[];
dsd2s_all=[];

figure

for k = 1:n_plot

    m=runs_to_plot(k);
    modelnumber = modelnumbers(k);

    if isempty(f2s{m}) || isempty(sd2s{m}) || isempty(V2s{m})
        continue; % Skip if model data is missing
    end

    % color map for plotting
    col = colourmap(k, :);

    % row index for all sites ('Holo-TwinCays_all' in testlocs.names)
    index = find(testlocs.reg==99999);

    % Calculate average RSL rates across all sites as derivatives of RSL
    difftimestep=100; % Temporal spacing (in years) between time steps for computing and plotting rate estimates.
    % Find time-paired rows separated by exactly 'difftimestep' units
    Mdiff = bsxfun(@eq,testlocs.X(index,3),testlocs.X(index,3)')-bsxfun(@eq,testlocs.X(index,3),testlocs.X(index,3)'+ difftimestep);
    % Remove rows without any matching time-difference pairs (i.e., only
    % want rows with difftimestep
    sub=find(sum(Mdiff,2)==0);  
    Mdiff=Mdiff(sub,:);
    % Compute average time difference (denominator for rate calculation)
    difftimes=1950-bsxfun(@rdivide,abs(Mdiff)*testlocs.X(index,3),sum(abs(Mdiff),2));
    % Normalize each row to compute finite differences (df/dt)
    Mdiff=bsxfun(@rdivide,Mdiff,Mdiff*testlocs.X(index,3));
    
    % Calculate rates
    df2s=Mdiff*f2s{m}{1}(index);
    dV2s=Mdiff*V2s{m}{1}(index,index)*Mdiff';
    dsd2s=sqrt(diag(dV2s(:,:)));

    % Plot RSL rates
    subplot(subplotrowno,subplotcolno,1);
    h=plot(difftimes,df2s,'Color',col,'Display',sprintf('Model %d',modelnumber)); % RSL rate mean
    legend_all=[legend_all; h];
    hold on
    fill([difftimes;flipud(difftimes)],[df2s+dsd2s; flipud(df2s-dsd2s)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL rate 1sd
    fill([difftimes;flipud(difftimes)],[df2s+2*dsd2s; flipud(df2s-2*dsd2s)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5) % RSL rate 2sd
    ylabel({'RSL Rate', '(mm/yr = m/ka)'}, 'Color', 'k');
    xlabel('Age (BP)')
    set(gca, 'XDir', 'reverse');% Reverse x axes 
    title(sprintf('Full RSL rate'))
    ratebuffer = (max(df2s+2*dsd2s)-min(df2s-2*dsd2s))*0.1;
    ylim([min(df2s-2*dsd2s)-ratebuffer max(df2s+2*dsd2s)+ratebuffer])

    % Calculate global rates
    df2s=Mdiff*f2s{m}{2}(index);
    dV2s=Mdiff*V2s{m}{2}(index,index)*Mdiff';
    dsd2s=sqrt(diag(dV2s(:,:)));

    % Append rates to calculate limits
    df2s_all = [df2s_all; df2s];
    dsd2s_all = [dsd2s_all; dsd2s];


    % Plot global RSL rates
    ageidx = find(difftimes<=RSLmaxage); % find index where prediction age is younger than max age of Belize SLIPs
    % filter predictions to time span of Belize
    difftimes_cut = difftimes(ageidx); 
    df2s_cut = df2s(ageidx);
    dsd2s_cut = dsd2s(ageidx);
    subplot(subplotrowno,subplotcolno,2);
    plot(difftimes_cut,df2s_cut,'Color',col); % RSL rate mean
    hold on
    fill([difftimes_cut;flipud(difftimes_cut)],[df2s_cut+dsd2s_cut; flipud(df2s_cut-dsd2s_cut)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL rate 1sd
    fill([difftimes_cut;flipud(difftimes_cut)],[df2s_cut+2*dsd2s_cut; flipud(df2s_cut-2*dsd2s_cut)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5) % RSL rate 2sd
    ylabel({'RSL Rate (mm/yr)'}, 'Color', 'k');
    xlabel('Age (BP)')
    set(gca, 'XDir', 'reverse');% Reverse x axes
    title(sprintf('Global'));

    % Set x-axis limits for global, regional and local subplots based on
    % age range of global component
    minage = min(difftimes_cut);
    maxage= max(difftimes_cut);

    % Plot regional components for each site
    for i=1:length(testlocs.names)-2 % loop for all local sites (remove last two rows as that is for Barbados and mean across Twin Cays)

        siteindex = find(testlocs.reg==testlocs.sites(i,1));

            % Calculate regional rates
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
        
            df2s=Mdiff*f2s{m}{3}(siteindex);
            dV2s=Mdiff*V2s{m}{3}(siteindex,siteindex)*Mdiff';
            dsd2s=sqrt(diag(dV2s(:,:)));

            % Append rates to calculate limits
            df2s_all = [df2s_all; df2s];
            dsd2s_all = [dsd2s_all; dsd2s];
        
            % Plot regional RSL rates
            ageidx = find(difftimes<=RSLmaxage); % find index where prediction age is younger than max age of Belize SLIPs
            % filter predictions to time span of Belize
            difftimes_cut = difftimes(ageidx); 
            df2s_cut = df2s(ageidx);
            dsd2s_cut = dsd2s(ageidx);
            subplot(subplotrowno,subplotcolno,2+i); % skip first two panels 
            plot(difftimes_cut,df2s_cut,'Color',col); % RSL rate mean
            hold on
            fill([difftimes_cut;flipud(difftimes_cut)],[df2s_cut+dsd2s_cut; flipud(df2s_cut-dsd2s_cut)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL rate 1sd
            fill([difftimes_cut;flipud(difftimes_cut)],[df2s_cut+2*dsd2s_cut; flipud(df2s_cut-2*dsd2s_cut)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5) % RSL rate 2sd
            ylabel({'RSL Rate (mm/yr)'}, 'Color', 'k');
            xlabel('Age (BP)')
            set(gca, 'XDir', 'reverse');% Reverse x axes
            title(sprintf('Regional (Site %d)',i));
    end

    % Plot local components for each site
    for i=1:length(testlocs.names)-2 % loop for all local sites (remove last two rows as that is for Barbados and mean across Twin Cays)

        siteindex = find(testlocs.reg==testlocs.sites(i,1));

            % Calculate local rates
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
        
            df2s=Mdiff*f2s{m}{4}(siteindex);
            dV2s=Mdiff*V2s{m}{4}(siteindex,siteindex)*Mdiff';
            dsd2s=sqrt(diag(dV2s(:,:)));

            % Append rates to calculate limits
            df2s_all = [df2s_all; df2s];
            dsd2s_all = [dsd2s_all; dsd2s];
        
            % Plot local RSL rates
            ageidx = find(difftimes<=RSLmaxage); % find index where prediction age is younger than max age of Belize SLIPs
            % filter predictions to time span of Belize
            difftimes_cut = difftimes(ageidx); 
            df2s_cut = df2s(ageidx);
            dsd2s_cut = dsd2s(ageidx);
            subplot(subplotrowno,subplotcolno,2+length(testlocs.names)-2+i); % skip first two panels and regional panels
            plot(difftimes_cut,df2s_cut,'Color',col); % RSL rate mean
            hold on
            fill([difftimes_cut;flipud(difftimes_cut)],[df2s_cut+dsd2s_cut; flipud(df2s_cut-dsd2s_cut)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL rate 1sd
            fill([difftimes_cut;flipud(difftimes_cut)],[df2s_cut+2*dsd2s_cut; flipud(df2s_cut-2*dsd2s_cut)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5) % RSL rate 2sd
            ylabel({'RSL Rate (mm/yr)'}, 'Color', 'k');
            xlabel('Age (BP)')
            set(gca, 'XDir', 'reverse');% Reverse x axes
            title(sprintf('Local (Site %d)',i));
    end

end

% set y axis limits
minrate = min(df2s_all-2*dsd2s_all);
maxrate = max(df2s_all+2*dsd2s_all);
buffer=(maxrate-minrate)*0.1;

% Loop through all subplots and apply the same y-limits
numplots = (subplotrowno * subplotcolno); % Total subplot slots
for sp = 2:numplots
    subplot(subplotrowno, subplotcolno, sp);
    ylim([minrate-buffer maxrate+buffer]);
    xlim([minage maxage])
end

% % Now add legend inside figure, below subplots
% lgd = legend(legend_all, 'Orientation', 'horizontal', 'NumColumns', ceil(length(legend_all)/legendrows));
% 
% % Adjust legend position (tweak these numbers as needed)
% lgd.Position = [0.3, 0.0, 0.4, 0.06]; % [left bottom width height]
    
% save figure as pdf 
fig = gcf;
fig.Units = 'inches';         % Set units to inches
if strcmpi(imageres, 'low') 
    fig.Position = [1, 1, 15, 8]; % [left, bottom, width, height] — adjust width (e.g., 10 inches)
    exportgraphics(gcf, strcat(filename,'_rates.pdf'), 'ContentType', 'image', 'Resolution', 300);
else 
    fig.Position = [1, 1, 15, 6]; % [left, bottom, width, height] — adjust width (e.g., 10 inches)
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 10 6]; % [left, bottom, width, height] — adjust width (e.g., 10 inches)
    fig.PaperSize = [10 6];
    print(gcf, strcat(filename,'_rates.pdf'), '-dpdf', '-painters', '-r300');
end

