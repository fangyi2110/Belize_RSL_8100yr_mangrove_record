    firstyears = [-6050:500:-2550];
    lastyears = [-5550:500:-2050];
    trainsub = find((wdataset.limiting==0));
    wmodelspec=modelspec;
%    Flat=[5.6];
%    Flong=[-55.75];
%     Flat=[6.82];
%     Flong=[-58.1];
%     Flat=5.35;
%     Flong=-55.4;
    results = [];
    fid=fopen(['HeightProb.tsv'],'w');
    fprintf(fid,['Height\tHeight SD\tProb\tfirstyear\tlastyear']);
    fprintf(fid,'\n');

    for i=1:length(firstyears)
        [fslopeF,sdslopeF,fsF,sdsF,~,~,~,~,passderivs,invcv,fmeanF,fsdF] = RegressRateField_ea(wdataset,wmodelspec,thetTGG{jj},noiseMasks(1,:),Flat,Flong,firstyears(i),lastyears(i),trainsub);
        fprobOverZero = normcdf(-meanSL,fmeanF,fsdF,'upper');
        fprintf(fid,'%0.2f',[fmeanF+meanSL]);
        fprintf(fid,'\t%0.2f',[fsdF]);
        fprintf(fid,'\t%0.4f',[fprobOverZero]);
        fprintf(fid,'\t%0.2f',[1950-firstyears(i)]);
        fprintf(fid,'\t%0.2f',[1950-lastyears(i)]);
        fprintf(fid,'\n');
    end;
    fclose(fid);
%     res = [fmeanF fsdF fprobOverZero 1950-firstyears(i) 1950-lastyears(i)];
%     results = [results res'];
% 
%     fid=fopen(['HeightProb.tsv'],'w');
%     fprintf(fid,['Height\tHeight SD\tProb\tfirstyear\tlastyear']);
%     fprintf(fid,'\n');
%         for pp=1:length(firstyears)
%             fprintf(fid,'\t%0.2f',[fslopeavg(kk,pp) sdslopeavg(kk,pp)]);
% %             fprintf(fid,'\t%0.2f',[fslopeavg(kk,pp) sdslopeavg(kk,pp)]);
% %             fprintf(fid,'\t%0.2f',[fslopeavg(kk,pp) sdslopeavg(kk,pp)]);           
%         end
%     %    for pp=1:length(diffplus)
%     %        fprintf(fid,'\t%0.2f',[fslopeavgdiff(kk,pp) 2*sdslopeavgdiff(kk,pp)]);
%     %    end
%         fprintf(fid,'\n');
%     end
%     fclose(fid);
