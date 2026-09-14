% EA (09/11/20): Added different time periods to average rates over
% Generate table of rates at each site, and difference sin rate at each site
% if exist('ka4')
%     fy = -2500;
% elseif min(testt)<0
%     fy=-1000;
% else
%     fy=0;
% end
% fslopeavg={};
% sdslopeavg={};
% iii=1;
% for kkk = 1:4
%     if kkk == 1
%         xx = '100yrs';    
%         firstyears = fy:100:1900;
%         lastyears= fy+100:100:2000;
%     elseif kkk == 2
%         xx = '40yrs'; 
%         firstyears = fy:40:1960;
%         lastyears= fy+40:40:2000;
%     elseif kkk == 3
%         xx = '20yrs';
%         firstyears = fy:20:1980;
%         lastyears= fy+20:20:2000;
%     elseif kkk == 4
%         xx = '60yrs';    
%         firstyears = fy+mod(1940-fy,60):60:1940;
%         lastyears= fy+60+mod(1940-fy,60):60:2000;
%     end
firstyears = 950:1000:11950;
lastyears = -50:1000:10950;

    fslopeavg={};
    sdslopeavg={};
    for jjj = 1:length(f2s{iii}(1,:))
    %for jjj = 1:length(5)
        % rate calculation
        % put this in the same format as f2s
        [fslopeavg{iii}(:,jjj),sdslopeavg{iii}(:,jjj)]=SLRateCalc(f2s{iii}(:,jjj),V2s{iii}(:,:,jjj),testsites,testreg,testX(:,3),firstyears,lastyears);
        %[fslopeavg{iii}(:,jjj),sdslopeavg{iii}(:,jjj),tavg,treg]=SLRateCalc(f2s{iii}(:,jjj),V2s{iii}(:,:,jjj),testsites,testreg,testX(:,3),firstyears,lastyears);
        %[fslopeavg{iii}(:,jjj),sdslopeavg{iii}(:,jjj),tavg,treg]=SLRateCalc(f2s{iii}(:,jjj),V2s{iii}(:,:,jjj),testsites,testreg,testX(:,3),firstyears,lastyears);
    end    
    fid=fopen(['sldecomp_rates_' xx '.tsv'],'w');
    fprintf(fid,'Site Name\tLat\tLong\tAvg Year');
    for ssss=1:length(noiseMasklabels)
        fprintf(fid,['\t' noiseMasklabels{ssss} ' (mm/yr)\t1s']);
    end
    fprintf(fid,'\n');

    for iii=1:length(fslopeavg)
        for rrrr=1:size(fslopeavg{iii},1) % times rate predicted
            subname=find(testsites(:,1)==treg(rrrr));
            %subname=treg(rrrr);
            fprintf(fid,testnames2{subname(1)});
%            fprintf(fid,num2str(testnames2(subname)));
            fprintf(fid,'\t%0.2f',testsites(subname(1),2:3));
            fprintf(fid,'\t%0.0f',tavg(rrrr)); 
            for ssss=1:size(fslopeavg{iii},2) % decomposition type
                fprintf(fid,'\t%0.2f',fslopeavg{iii}(rrrr,ssss));
                fprintf(fid,'\t%0.2f',sdslopeavg{iii}(rrrr,ssss));        
        %         fprintf(fid,'\t%0.2f',f2s{iii}(rrrr,ssss));
        %         fprintf(fid,'\t%0.2f',sd2s{iii}(rrrr,ssss));        
            end
            fprintf(fid,'\n');
        end
    end
    fclose(fid);
%end

% fid=fopen(['linrates' labl '.tsv'],'w');
% fprintf(fid,['Rates (mm/y), ' labl '\n']);
% fprintf(fid,'Site\tSiteID\tLat\tLong\tOldest\tYoungest\tICE5G VM2-90\tRate (linear)\t2s\tOffset (mm)\t2s');
% for pp=1:length(firstyears)
%     fprintf(fid,'\tRate (avg, %0.0f-%0.0f)\t2s\tP>0',[firstyears(pp) lastyears(pp)]);
% end
% for pp=1:length(diffplus)
%     fprintf(fid,['\tRate Diff. (avg, %0.0f-%0.0f minus ' ...
%                     '%0.0f-%0.0f)\t2s\tP>0'],[firstyears(diffplus(pp)) ...
%                         lastyears(diffplus(pp)) firstyears(diffless(pp)) lastyears(diffless(pp))]);
% end
% 
% fprintf(fid,'\n');
% fprintf(fid,noiseMasklabels{1});
% fprintf(fid,'\n');
% for kk=1:size(testsites,1)
%     fprintf(fid,testnames2{kk});
%     fprintf(fid,'\t%0.2f',testsitedef.sites(kk,:));
%     fprintf(fid,'\t%0.0f',testsitedef.oldest(kk));
%     fprintf(fid,'\t%0.0f',testsitedef.youngest(kk));       
%     fprintf(fid,'\t%0.2f',testsitedef.GIA(kk));
%     fprintf(fid,'\t%0.2f',[fslopelin(kk) 2*sdslopelin(kk)]);
%     fprintf(fid,'\t%0.2f',[foffset(kk) 2*sdoffset(kk)]);
%     for pp=1:length(firstyears)
%         fprintf(fid,'\t%0.2f',[fslopeavg(kk,pp) 2*sdslopeavg(kk,pp)]);
%         fprintf(fid,'\t%0.3f',[normcdf(fslopeavg(kk,pp)/sdslopeavg(kk,pp))]);
%     end
%     for pp=1:length(diffplus)
%         fprintf(fid,'\t%0.2f',[fslopeavgdiff(kk,pp) 2*sdslopeavgdiff(kk,pp)]);
%         fprintf(fid,'\t%0.3f',[normcdf(fslopeavgdiff(kk,pp)/sdslopeavgdiff(kk,pp))]);
%     end
%     fprintf(fid,'\n');
% end
% 
% fclose(fid);


