%% Test cases definition, index is vector of scalar test case number
% Each test case is defined by creating objects for the basic RAT settings,
% enviroment to use and the channel and interference settings together with
% simulator settings and what input and output to be used. There is one
% template / generic testcase that is default. Each testcase defines an
% overlay of parameters. If parameters not are defined in a test case the
% generic ones will be used.

function tc = tc(index)

    tc = cell(length(index), 1);
    
    %% Default test case definition - later test cases overrides this one
    t_g = struct;                % define a structure array
    t_g.name = 'tc_default';     % 
    t_g.description = 'Default parameters';
    
    % numerology (signal)
    t_g.nr = nr(2);
    t_g.slotNumber = randi(20,1)-1;
    
    % TRPs setting
    t_g.env_fix.trpselect='Auto';
    t_g.env_fix.num_trp=6;                                                 % number of trps used
    t_g.env_fix.trp_pos = [...
        10, 30, 50, 70, 90, 110, 10, 30, 50, 70, 90, 110;...
        15, 15, 15, 15, 15,  15, 35, 35, 35, 35, 35,  35;...
        3,  3,  3,  3,  3,  3,   3,  3,  3,  3,  3,   3];
    t_g.env_fix.trp_bearing = [30; 150; -90]-90;
    t_g.env_fix.max_dim = [120;50;3];                                      % The size of the office
    
    % UE setting
    t_g.env_fix.ue_pos = 'random';                                          % UEs position mode: 'random' or 'fixed' 
    t_g.env_fix.ue_nbr = 10;                                                % number of UEs to run
    
    % Channel setting
    t_g.env_fix.ch = 'cdl-d';                                                % channel mode
    
    % decoder setting
    t_g.ue_conf.combiner = 1;                                              %Zero forcing
    t_g.ue_conf.eta = 0.8; % 10^(-0.15);    % First to peak ratio to detect the first path % to determine the threshold, a small test can be done where running the simulation only for AWGN and then for ETU then take a value which will b ealmost inbetween based on CDF plots
    t_g.ue_conf.dr = 1;
    t_g.ue_conf.searchGrid = 100;           % this parameter is determinded by the range of the scnario!!! 
    t_g.ue_conf.doMPS = true;
    t_g.ue_conf.corr_method = 'time';
    t_g.ue_conf.interp_fac=1;
    % local server setting 
    t_g.ls_conf.method = 'one';
    % data 
    t_g.sink_conf.save_file = 'false';
    t_g.sink_conf.debug_store = 'false';
    t_g.sink_conf.plot_cdl = 'true';
    t_g.sink_conf.plot_geom = 'true';    
    
    for i = 1:length(index)
        
        t = t_g;
        switch index(i)
            case 1
                
            %% Test case 1 - Stationary indoor FR1 NR 
            % Define NR with 15-kHz sub-carrier spacing numerology
            t.name = 'tc1';
            t.description = 'Indoor open office, DL OTDOA, FR1 - 3 GHz, 15kHz SCS';

            otherwise
                error('Expected test case not defined - abort')
        end
        
        % Add to list of test cases
         tc{i,1} = t;
         
    end
    
end
