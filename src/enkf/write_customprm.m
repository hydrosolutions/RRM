function [E] = write_customprm(paths,prm,Eout)
% Writing back updated custom parameters after successful assimilation.
%
% tobias siegfried, 26/11/2014
%% Custom prm file
cd(strcat(paths.main,'prm/'))
fid = fopen('prm-Custom.txt', 'w');
% print back file values
fprintf(fid, 'customprm.dt=%d \n',prm.customprm.dt);
%fprintf(fid, 'customprm.T=%d \n',prm.customprm.T);
%fprintf(fid, 'customprm.use_compiled=%d \n',prm.customprm.use_compiled);
fprintf(fid, 'customprm.limCutoff=%e \n',prm.customprm.limCutoff);
fprintf(fid, 'customprm.nStorage=%d \n',prm.customprm.nStorage);
fprintf(fid, 'customprm.simStartY=%d \n',prm.customprm.simStartY);
fprintf(fid, 'customprm.simStartM=%d \n',prm.customprm.simStartM);
fprintf(fid, 'customprm.simStartD=%d \n',prm.customprm.simStartD);
fprintf(fid, 'customprm.TmeanAvail=%d\n',prm.customprm.TmeanAvail);
fprintf(fid, 'customprm.obsVarQ=%f\n',prm.customprm.obsVarQ);
fprintf(fid, 'customprm.obsVarS=%f\n',prm.customprm.obsVarS);
fprintf(fid, 'customprm.obsVarG=%f\n',prm.customprm.obsVarG);
fprintf(fid, 'customprm.obsVarSN=%f\n',prm.customprm.obsVarSN);
fprintf(fid, 'customprm.obsVarET=%f\n',prm.customprm.obsVarET);
fprintf(fid, 'customprm.snow=%d\n',prm.customprm.snow);
fprintf(fid, 'customprm.calculateET=%d\n',prm.customprm.calculateET);
fprintf(fid, 'customprm.fclength=%f\n',prm.customprm.fclength);
fprintf(fid, 'customprm.Fmin=%f\n',prm.customprm.Fmin);
fprintf(fid, 'customprm.Fmax=%f\n',prm.customprm.Fmax);
%eval(['fprintf(fid, ''customprm.integrator=%s \n'', ''''''' prm.customprm.integrator ''''''');']);
eval(['fprintf(fid, ''customprm.fname_samples=%s \n'', ''''''' prm.customprm.fname_samples ''''''');']);
eval(['fprintf(fid, ''customprm.loc_connMat=%s \n'', ''''''' prm.customprm.loc_connMat ''''''');']);
eval(['fprintf(fid, ''customprm.locS0G0=%s \n'', ''''''' prm.customprm.locS0G0 ''''''');']);
eval(['fprintf(fid, ''customprm.runType=%s\n'', ''''''' prm.customprm.runType ''''''');']);
eval(['fprintf(fid, ''customprm.geometry=%s\n'', ''''''' prm.customprm.geometry ''''''');']);
fclose(fid);
%% WRITE BACK UPDATED E and SG MATRICES
E = Eout; save(prm.path.E,'E');
end
