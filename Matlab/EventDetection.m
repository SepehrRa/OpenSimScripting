function Events = EventDetection(Trialname,FTable,ForceRatio,MTable,Threshold)
if contains(Trialname,"IsoM")
    [indx,c]=find(abs(FTable.data(:,10))>ForceRatio.*max(abs(FTable.data(:,10))));
    Stime=FTable.data(indx([1;find(diff(indx)>10)+1]),1);
    Etime=FTable.data(indx([(find(diff(indx)>10));end]),1);
elseif contains(Trialname,"IsoK") 
    if contains(Trialname,"Fl")
        [indx,c]=find(MTable.data(:,6)>Threshold(1) & MTable.data(:,6)<Threshold(2));
        % & [0;diff(MTable.data(:,6))>0]
        Sindx=indx([1;find(diff(indx)>10)+1]);
        Eindx=indx([(find(diff(indx)>10));end]);
        count=0;
        Stime=[];
        Etime=[];
        for ww=1:length(Sindx)
            if (MTable.data(Sindx(ww),6)-MTable.data(Sindx(ww)-1,6))>0
                count=count+1;
                Stime=[Stime;FTable.data(Sindx(ww),1)];
                Etime=[Etime;FTable.data(Eindx(ww),1)];
            end
        end
    else
        [indx,c]=find(MTable.data(:,6)>Threshold(1) & MTable.data(:,6)<Threshold(2));
        Sindx=indx([1;find(diff(indx)>10)+1]);
        Eindx=indx([(find(diff(indx)>10));end]);
        count=0;
        Stime=[];
        Etime=[];
        for ww=1:length(Sindx)
            if (MTable.data(Sindx(ww)+1,6)-MTable.data(Sindx(ww),6))<0
                count=count+1;
                Stime=[Stime;FTable.data(Sindx(ww),1)];
                Etime=[Etime;FTable.data(Eindx(ww),1)];
            end
        end
    end
else
     fprintf('\nERROR: %s unknown trail ...\n\n', Trialname);
end
Events=[Stime,Etime];
end