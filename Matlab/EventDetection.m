function Events = EventDetection(Trialname,FTable,ForceRatio,MTable,Threshold)
if contains(Trialname,"IsoM")
    [indx,c]=find(abs(FTable(:,10))>ForceRatio.*max(abs(FTable(:,10))));
    Stime=FTable(indx([1;find(diff(indx)>10)+1]),1);
    Etime=FTable(indx([(find(diff(indx)>10));end]),1);
elseif contains(Trialname,"IsoK") 
    if contains(Trialname,"Fl")
        [indx,c]=find(MTable(:,6)>Threshold(1) & MTable(:,6)<Threshold(2));
        % & [0;diff(MTable(:,6))>0]
        Sindx=indx([1;find(diff(indx)>10)+1]);
        Eindx=indx([(find(diff(indx)>10));end]);
        count=0;
        Stime=[];
        Etime=[];
        for ww=1:length(Sindx)
            if (MTable(Sindx(ww),6)-MTable(Sindx(ww)-1,6))>0
                count=count+1;
                Stime=[Stime;FTable(Sindx(ww),1)];
                Etime=[Etime;FTable(Eindx(ww),1)];
            end
        end
    else
        [indx,c]=find(MTable(:,6)>Threshold(1) & MTable(:,6)<Threshold(2));
        Sindx=indx([1;find(diff(indx)>10)+1]);
        Eindx=indx([(find(diff(indx)>10));end]);
        count=0;
        Stime=[];
        Etime=[];
        for ww=1:length(Sindx)
            if (MTable(Sindx(ww)+1,6)-MTable(Sindx(ww),6))<0
                count=count+1;
                Stime=[Stime;FTable(Sindx(ww),1)];
                Etime=[Etime;FTable(Eindx(ww),1)];
            end
        end
    end
else
     fprintf('\nERROR: %s unknown trail ...\n\n', Trialname);
end
Events=[Stime,Etime];
end