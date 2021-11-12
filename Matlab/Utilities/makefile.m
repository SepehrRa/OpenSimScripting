function makefile (folder,F_fnames,Dataheader,Datadata,Resulotion,delimiterIn)
% Resulotion number of digit after zero
Datah=[];
            [r,c]=size(Datadata);
            for t=1:length(Dataheader)
            Datah=[Datah char(Dataheader(t)),'\t']
            end
            Title='\nversion=1\nnRows=%d\nnColumns=%d\nInDegrees=no\nendheader\n';
            fid=fopen(append(folder,'\',F_fnames), "w");
            if fid < 0
                fprintf('\nERROR: %s could not be opened for writing...\n\n', fname);
            return
            end
            [r,c]=size(Datadata);
            ft=[];
            for f=1:c
                ft=[ft '%.' num2str(Resulotion) 'f' delimiterIn];
            end
            fprintf(fid,[char(F_fnames) Title],[r,c]);
            fprintf(fid,[Datah '\n']);
                for i = 1:r
                    fprintf(fid,ft,Datadata(i,:));
                    fprintf(fid, '\n');
                end
            fclose(fid);
            fprintf('Saved %s\n', append(folder,'\',F_fnames))
end