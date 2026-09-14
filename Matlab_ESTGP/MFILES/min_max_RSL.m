    firstyears = [-5050:100:-3050];
    lastyears = [-4950:100:2-950];

    fid=fopen(['HeightMinMax_small.tsv'],'w');
    fprintf(fid,['Max or Min\tLat\tLong\tHeight\tHeight SD\tfirstyear\tlastyear']);
    fprintf(fid,'\n');

for i=1:length(firstyears)
%    Flat=5:.3:27;
%    Flong=-95:.3:-50;
    sub=find(testsites(:,2)<=360);
    ulat=unique(round(testsites(sub,2)));
    ulong=unique(round(testsites(sub,3)));
    ulat=unique(bsxfun(@plus,ulat,[-4:2:4]));
    ulong=unique(bsxfun(@plus,ulong,[-4:2:4]));
%    sub=find(abs(ulat)<90); ulat=ulat(sub); ulat=ulat(:)'; ulong=ulong(:)';
    Flat=ulat;
    Flong=ulong;
    %union(mod(Flong,360),mod(ulong,360));
    wmodelspec = modelspec;
    [fslopeF,sdslopeF,fsF,sdsF,~,~,~,~,passderivs,invcv,fmeanF,fsdF] = RegressRateField_ea(wdataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),Flat,Flong,firstyears(i),lastyears(i),trainsub);
    [FLONG,FLAT]=meshgrid(Flong,Flat);
    submax=find(max(fmeanF)==fmeanF);
    submin=find(min(fmeanF)==fmeanF);
    maxH=fmeanF(submax)+meanSL;
    maxSD=fsdF(submax);
    maxLat=FLAT(submax);
    maxLong=FLONG(submax);
    minH=fmeanF(submin)+meanSL;
    minSD=fsdF(submin);
    minLat=FLAT(submin);
    minLong=FLONG(submin);
        fprintf(fid,['Max']);
        fprintf(fid,'\t%0.2f',[maxLat]);
        fprintf(fid,'\t%0.2f',[maxLong]);
        fprintf(fid,'\t%0.2f',[maxH]);
        fprintf(fid,'\t%0.2f',[maxSD]);
        fprintf(fid,'\t%0.2f',[1950-firstyears(i)]);
        fprintf(fid,'\t%0.2f',[1950-lastyears(i)]);
        fprintf(fid,'\n');
        fprintf(fid,['Min']);
        fprintf(fid,'\t%0.2f',[minLat]);
        fprintf(fid,'\t%0.2f',[minLong]);
        fprintf(fid,'\t%0.2f',[minH]);
        fprintf(fid,'\t%0.2f',[minSD]);
        fprintf(fid,'\t%0.2f',[firstyears(i)]);
        fprintf(fid,'\t%0.2f',[lastyears(i)]);
        fprintf(fid,'\n');
end;
fclose(fid);


    
    
    %[fslopeGSL,sdslopeGSL]=SLRateCompare(f2s{iii}(GSLdatsub,1),V2s{iii}(GSLdatsub,GSLdatsub,1),testsites(GSLsitesub),testreg(GSLdatsub),testX(GSLdatsub,3),firstyears,lastyears);
%    [priorslope,sdpriorslope,priorfsF,priorsdsF] = RegressRateField(nulldataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),-80,20,firstyears,lastyears);

    Flat1=min(Flat):max(Flat);
    Flong1=min(Flong):max(Flong);

    [FLONG,FLAT]=meshgrid(Flong,Flat);
    [FLONG1,FLAT1]=meshgrid(Flong1,Flat1);
    mapped = griddata(FLONG(:),FLAT(:),fmeanF(:)+meanSL,Flong,Flat(:),'linear'); %+fslopeGSL;
    sdmapped = griddata(FLONG(:),FLAT(:),fsdF(:),Flong,Flat(:),'linear');

