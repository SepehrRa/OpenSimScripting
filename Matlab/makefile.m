function makefile (folder,F_fnames,Title,Titledata,Dataheader,Datadata)
           fid=fopen([folder '\' char(F_fnames)], "w");
            if fid < 0
                fprintf('\nERROR: %s could not be opened for writing...\n\n', fname);
            return
            end
            [r,c]=size(Datadata);
            ft=[];
            for f=1:c
                ft=[ft '%.5f\t'];
            end
            fprintf(fid,[char(F_fnames) Title],Titledata);
            fprintf(fid,[Dataheader '\n']);
                for i = 1:r
                    fprintf(fid,ft,Datadata(i,:));
                    fprintf(fid, '\n');
                end
            fclose(fid);
            fprintf('Saved %s\n', [folder char(F_fnames)])
end