% Last updated by Erica Ashe, Tuesday, 10-4-2022

    minlat=-36;
    maxlat=-13;
    minlong=112;
    maxlong=128;
    [minlat maxlat minlong maxlong]
    lat_incr = .45;
    long_incr = .45;
    Flat=minlat:lat_incr:maxlat;
    Flong=minlong:long_incr:maxlong;
    scattersize = 200;

nulldataset=SubsetDataStructure(wdataset,1,1);
nulldataset.meantime=0; nulldataset.dt=0; nulldataset.dY=200e3; nulldataset.limiting=0;

trainsub = find((wdataset.limiting==0));
firstyears=-10050:500:1450; %[-8550:1000:%4200bp to 4400bp
lastyears=-9550:500:1950;%[-7550:1000:
for ii=1:length(firstyears)
    nsitepts{ii}=[];
    td = lastyears(ii)-firstyears(ii);
    addt=(3000-td)/2;
    fy=firstyears(ii)-addt;
    ly=lastyears(ii)+addt;
    subt = intersect(find(wdataset.time1<=ly),find(wdataset.time2>=fy));
    %find the number of sites with data influencing the model close to this time period
    wmodelspec = modelspec;

    [fslopeF,sdslopeF,fsF,sdsF,~,~,~,~,passderivs,invcv,fmeanF,fsdF] = RegressRateField_ea(wdataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),Flat,Flong,firstyears(ii),lastyears(ii),trainsub);
    [priorslope,sdpriorslope,priorfsF,priorsdsF] = RegressRateField_ea(nulldataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),mean(Flat),mean(Flong),firstyears(ii),lastyears(ii));

    fprobOverZero = normcdf(-meanSL,fmeanF,fsdF,'upper');
% Start mapping %%%

    clf;
    ax = worldmap([minlat maxlat],[minlong maxlong]);
    setm(ax,'meridianlabel','off','parallellabel','off','flinewidth',3);
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(ax, land, 'FaceColor', [0.7 0.7 0.7]);
    hold on;

    Flat1=min(Flat):max(Flat);
    Flong1=min(Flong):max(Flong);

    [FLONG,FLAT]=meshgrid(Flong,Flat);
    [FLONG1,FLAT1]=meshgrid(Flong1,Flat1);
    mapped = griddata(FLONG(:),FLAT(:),fprobOverZero(:),Flong,Flat(:),'linear'); %+fslopeGSL;
    sdmapped = griddata(FLONG(:),FLAT(:),sdslopeF(:),Flong,Flat(:),'linear');

    hs1=scatterm(FLAT(:),FLONG(:),scattersize,mapped(:),'filled','marker','s');
    hold on;

    sublong=find((mod(FLONG1(:),5)==0).*(mod(FLAT1(:),5)==0));

    geoshow(ax, land, 'FaceColor', [0.9 0.9 0.9]);
    hold on;
    
    subt = intersect(find(wdataset.time1<=ly),find(wdataset.time2>=fy));
    ud=unique(wdataset.datid(find((wdataset.time2>=firstyears(ii)).*(wdataset.time1<=lastyears(ii)))));
    sub1=find(ismember(wdataset.siteid,ud));

    % map the regions where there are data and how much
    for j=1:length(wdataset.siteid)
      nsitepts{ii}(j)=length(intersect(intersect(intersect(find(wdataset.limiting==0),find(wdataset.time1<=ly)),find(wdataset.time2>=fy)),find(wdataset.datid==wdataset.siteid(j))));
    end 

    subsite = find(nsitepts{ii}>0);
    ndat = nsitepts{ii}(subsite)';
    lat_site = wdataset.sitecoords(subsite,1);
    long_site = wdataset.sitecoords(subsite,2);

    uns = unique(ndat);
    ss_1_3=intersect(find(ndat>=1),find(ndat<=3));
    ss_4_10=intersect(find(ndat>=4),find(ndat<=10));
    ss_11_20=intersect(find(ndat>=11),find(ndat<=20));
    ss_21=find(ndat>=21);

    hgh5=scatterm(lat_site(ss_21),long_site(ss_21),5*ndat(ss_21),'MarkerEdgeColor','k','MarkerFaceColor','w');
    hgh3=scatterm(lat_site(ss_11_20),long_site(ss_11_20),5*ndat(ss_11_20),'MarkerEdgeColor','k','MarkerFaceColor','w');
    hgh2=scatterm(lat_site(ss_4_10),long_site(ss_4_10),5*ndat(ss_4_10),'MarkerEdgeColor','k','MarkerFaceColor','w');
    hgh1=scatterm(lat_site(ss_1_3),long_site(ss_1_3),5*ndat(ss_1_3),'MarkerEdgeColor','k','MarkerFaceColor','w');
    hgh4=scatterm(lat_site(ss_1_3),long_site(ss_1_3),5*ndat(ss_1_3),'MarkerEdgeColor','k','MarkerFaceColor','w');

    colormap(jet);
    axis tight;
    hcb=colorbar;

    box on;
    caxis([0 1]);

    datehs = (1950-firstyears(ii)+1950-lastyears(ii))/2;
    datefirst = 1950-firstyears(ii);
    datelast = 1950-lastyears(ii);

    title({['Probability of Sea Level > 0, from ' num2str(datefirst) ' to ' num2str(datelast) ' BP']});
    pdfwrite(['ProbHighstand_' labl '_' num2str(datelast)]);
end
