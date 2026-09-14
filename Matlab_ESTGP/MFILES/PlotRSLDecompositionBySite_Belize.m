function[]=PlotRSLDecompositionSummarised(runs_to_plot,maxage,minage,datasets,f2s,sd2s,V2s,testlocs,legendrows,modelnumbers,subplotcolno,filename,imageres)

% convert max and min age to years BP
    maxageBP=1950-maxage;
    minageBP=1950-minage;

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

    % set y axis limits to be standardised for global and local components
    minRSL = min(min(cell2mat(cellfun(@(c) [c{3} c{2}], f2s(runs_to_plot), 'UniformOutput', false))-cell2mat(cellfun(@(c) [c{3} c{2}], sd2s(runs_to_plot), 'UniformOutput', false))))/1000;
    maxRSL = max(max(cell2mat(cellfun(@(c) [c{3} c{2}], f2s(runs_to_plot), 'UniformOutput', false))+cell2mat(cellfun(@(c) [c{3} c{2}], sd2s(runs_to_plot), 'UniformOutput', false))))/1000;
    buffer=(maxRSL-minRSL)*0.2;

    % row index for all sites ('Holo-TwinCays_all' in testlocs.names)
    index = which(testlocs.reg==99999);
    % stepsize=length(testlocs.X)/length(testlocs.names);
    % index = length(testlocs.X)-stepsize+1:length(testlocs.X); % last set of data are for all sites combined
    
    %% Plot RSL

    % Plot average RSL across all sites 
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
    xlim([minageBP maxageBP]); % set x axis limits 
    ylabel('RSL (m)','Color','k');
    xlabel('Age (BP)')
    set(gca, 'XDir', 'reverse');% Reverse x axes
    title('Mean RSL for Belize')

    % Plot global RSL component
    subplot(subplotrowno,subplotcolno,2);
    t=1950-testlocs.X(index,3);
    RSLmean=f2s{m}{2}(index)/1000; % f2s{m}{1}: the {1} extracts the full posterior output without white
    % noise; use {2} for global component only {3} for regional non-linear component only and {4} for local component only
    RSLsd = sd2s{m}{2}(index)/1000;
    plot(t,RSLmean,'Color',col); % RSL mean line
    hold on
    fill([t; flipud(t)], [RSLmean+RSLsd; flipud(RSLmean-RSLsd)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 1sd
    fill([t; flipud(t)], [RSLmean+2*RSLsd; flipud(RSLmean-2*RSLsd)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 2sd
    xlim([minageBP maxageBP]); % set x axis limits
    ylim([minRSL-buffer maxRSL+buffer])
    ylabel('RSL (m)','Color','k');
    xlabel('Age (BP)')
    set(gca, 'XDir', 'reverse');% Reverse x axes
    title('Global')

     % Plot regional components for each site
    for i=1:length(testlocs.names)-1 % loop for all local sites (remove last row as that is for mean across Twin Cays)
        
        stepsize=length(testlocs.X)/length(testlocs.names); % find no. of rows per site
        siteindex = stepsize*(i-1)+1:stepsize*i; % rows for site
        
        subplot(subplotrowno,subplotcolno,2+i); % skip first two panels 
        t=1950-testlocs.X(siteindex,3);
        RSLmean=f2s{m}{3}(siteindex)/1000; % f2s{m}{1}: the {1} extracts the full posterior output without white
        % noise; use {2} for global component only {3} for regional non-linear component only and {4} for local component only
        RSLsd = sd2s{m}{3}(siteindex)/1000;
        plot(t,RSLmean,'Color',col); % RSL mean line
        hold on
        fill([t; flipud(t)], [RSLmean+RSLsd; flipud(RSLmean-RSLsd)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 1sd
        fill([t; flipud(t)], [RSLmean+2*RSLsd; flipud(RSLmean-2*RSLsd)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 2sd
        xlim([minageBP maxageBP]); % set x axis limits 
        ylim([minRSL-buffer maxRSL+buffer])
        ylabel('RSL (m)','Color','k');
        xlabel('Age (BP)')
        set(gca, 'XDir', 'reverse');% Reverse x axes
        title(sprintf('Regional (Site %d)',i));
    end
    
    % Plot local components for each site
    for i=1:length(testlocs.names)-1 % loop for all local sites (remove last row as that is for mean across Twin Cays)
        
        stepsize=length(testlocs.X)/length(testlocs.names); % find no. of rows per site
        siteindex = stepsize*(i-1)+1:stepsize*i; % rows for site
        
        subplot(subplotrowno,subplotcolno,2+length(testlocs.names)-1+i); % skip first two panels and regional panels
        t=1950-testlocs.X(siteindex,3);
        RSLmean=f2s{m}{4}(siteindex)/1000; % f2s{m}{1}: the {1} extracts the full posterior output without white
        % noise; use {2} for global component only {3} for regional non-linear component only and {4} for local component only
        RSLsd = sd2s{m}{4}(siteindex)/1000;
        plot(t,RSLmean,'Color',col); % RSL mean line
        hold on
        fill([t; flipud(t)], [RSLmean+RSLsd; flipud(RSLmean-RSLsd)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 1sd
        fill([t; flipud(t)], [RSLmean+2*RSLsd; flipud(RSLmean-2*RSLsd)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL 2sd
        xlim([minageBP maxageBP]); % set x axis limits 
        ylim([minRSL-buffer maxRSL+buffer])
        ylabel('RSL (m)','Color','k');
        xlabel('Age (BP)')
        set(gca, 'XDir', 'reverse');% Reverse x axes
        title(sprintf('Local (Site %d)',i));
    end
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
    stepsize=length(testlocs.X)/length(testlocs.names);
    index = length(testlocs.X)-stepsize+1:length(testlocs.X); % last set of data are for all sites combined

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

    % Append rates to calculate limits
    df2s_all = [df2s_all; df2s];
    dsd2s_all = [dsd2s_all; dsd2s];

    % Plot RSL rates
    subplot(subplotrowno,subplotcolno,1);
    h=plot(difftimes,df2s,'Color',col,'Display',sprintf('Model %d',modelnumber)); % RSL rate mean
    legend_all=[legend_all; h];
    hold on
    fill([difftimes;flipud(difftimes)],[df2s+dsd2s; flipud(df2s-dsd2s)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL rate 1sd
    fill([difftimes;flipud(difftimes)],[df2s+2*dsd2s; flipud(df2s-2*dsd2s)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5) % RSL rate 2sd
    xlim([minageBP maxageBP]); % set x axis limits but let y axis float
    ylabel({'RSL Rate', '(mm/yr = m/ka)'}, 'Color', 'k');
    xlabel('Age (BP)')
    set(gca, 'XDir', 'reverse');% Reverse x axes 
    title(sprintf('Mean RSL rate for Belize'))

    % Calculate global rates
    df2s=Mdiff*f2s{m}{2}(index);
    dV2s=Mdiff*V2s{m}{2}(index,index)*Mdiff';
    dsd2s=sqrt(diag(dV2s(:,:)));

    % Append rates to calculate limits
    df2s_all = [df2s_all; df2s];
    dsd2s_all = [dsd2s_all; dsd2s];

    % Plot global RSL rates
    subplot(subplotrowno,subplotcolno,2);
    plot(difftimes,df2s,'Color',col); % RSL rate mean
    hold on
    fill([difftimes;flipud(difftimes)],[df2s+dsd2s; flipud(df2s-dsd2s)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL rate 1sd
    fill([difftimes;flipud(difftimes)],[df2s+2*dsd2s; flipud(df2s-2*dsd2s)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5) % RSL rate 2sd
    xlim([minageBP maxageBP]); % set x axis limits but let y axis float
    ylabel({'RSL Rate (mm/yr)'}, 'Color', 'k');
    xlabel('Age (BP)')
    set(gca, 'XDir', 'reverse');% Reverse x axes
    title(sprintf('Global'));

    % Plot regional components for each site
    for i=1:length(testlocs.names)-1 % loop for all local sites (remove last row as that is for mean across Twin Cays)

        stepsize=length(testlocs.X)/length(testlocs.names); % find no. of rows per site
        siteindex = stepsize*(i-1)+1:stepsize*i; % rows for site

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
            subplot(subplotrowno,subplotcolno,2+i); % skip first two panels 
            plot(difftimes,df2s,'Color',col); % RSL rate mean
            hold on
            fill([difftimes;flipud(difftimes)],[df2s+dsd2s; flipud(df2s-dsd2s)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL rate 1sd
            fill([difftimes;flipud(difftimes)],[df2s+2*dsd2s; flipud(df2s-2*dsd2s)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5) % RSL rate 2sd
            xlim([minageBP maxageBP]); % set x axis limits but let y axis float
            ylabel({'RSL Rate (mm/yr)'}, 'Color', 'k');
            xlabel('Age (BP)')
            set(gca, 'XDir', 'reverse');% Reverse x axes
            title(sprintf('Regional (Site %d)',i));
    end

    % Plot local components for each site
    for i=1:length(testlocs.names)-1 % loop for all local sites (remove last row as that is for mean across Twin Cays)

        stepsize=length(testlocs.X)/length(testlocs.names); % find no. of rows per site
        siteindex = stepsize*(i-1)+1:stepsize*i; % rows for site

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
            subplot(subplotrowno,subplotcolno,2+length(testlocs.names)-1+i); % skip first two panels and regional panels
            plot(difftimes,df2s,'Color',col); % RSL rate mean
            hold on
            fill([difftimes;flipud(difftimes)],[df2s+dsd2s; flipud(df2s-dsd2s)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5); % RSL rate 1sd
            fill([difftimes;flipud(difftimes)],[df2s+2*dsd2s; flipud(df2s-2*dsd2s)],col,'FaceAlpha',0.1,'EdgeColor',col,'EdgeAlpha',0.5) % RSL rate 2sd
            xlim([minageBP maxageBP]); % set x axis limits but let y axis float
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
for sp = 1:numplots
    subplot(subplotrowno, subplotcolno, sp);
    ylim([minrate-buffer maxrate+buffer]);
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

