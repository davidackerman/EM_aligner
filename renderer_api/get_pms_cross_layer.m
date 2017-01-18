function jj = get_pms_cross_layer(pm, sID1, sID2, wopts)
%% get point-matches between two groups (layers) sID1 and sID2
%%%%%%%%



urlChar = sprintf('%s/owner/%s/matchCollection/%s/group/%s/matchesWith/%s', ...
    pm(1).server, pm(1).owner, pm(1).match_collection, sID1, sID2);

U = matlab.net.URI(urlChar);

if numel(pm)>1
    data_options_str = '?';
    for ix = 2:numel(pm)
        if ix ==2
            data_options_str = [data_options_str 'mergeCollection=' pm(ix).match_collection];
        else
            data_options_str = [data_options_str '&mergeCollection=' pm(ix).match_collection];
        end
    end
    QPs = matlab.net.QueryParameter(data_options_str);
    U.Query = QPs;
end

try
    jj = webread(char(U), wopts);
catch err_fetch_pm
    kk_disp_err(err_fetch_pm)
    pause(1);
    disp('trying again');
    jj = webread(char(U),options); % try again
end



% % %disp([sID{isix} ' ---- ' sID{jsix}]);
% % %%jj = webread(urlChar, wopts);
% % try
% %     jj = webread(urlChar, wopts);
% % catch err_fetch_pm
% %     kk_disp_err(err_fetch_pm)
% %     pause(1);
% %     jj = webread(urlChar,wopts); % try again
% % end

%%%%%%%%%%%%%%%
