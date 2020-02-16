%% Test cases definition, index is vector of scalar test case number
function tc = tc(index)
    %% Default test case definition - later test cases overrides this one
    t = struct;                % define a structure array
    t.name = 'tc1';
    t.description = 'Indoor open office, DL OTDOA, FR1 - 3 GHz, 15kHz SCS';
    t.nr = nr(1);
    % Define the environment
    trp_pos = [...
        10, 30, 50, 70, 90, 110, 10, 30, 50, 70, 90, 110;...
        15, 15, 15, 15, 15,  15, 35, 35, 35, 35, 35,  35;...
        3,  3,  3,  3,  3,  3,   3,  3,  3,  3,  3,   3];
    trp_bearing = [30; 150; -90] - 90;
    ue_nbr = 10;
    ue_pos = 'random';
    max_dim = [120;50;3];
    ch = 'los'; %%%  
    t.env = environment(trp_pos, trp_bearing, max_dim, ch, ue_nbr, ue_pos);
    t.slotNumber = randi(20,1)-1;

    t.ue_conf.combiner = 1;                                              %Zero forcing
    t.ue_conf.eta = 0.8; % 10^(-0.15);    % First to peak ratio to detect the first path % to determine the threshold, a small test can be done where running the simulation only for AWGN and then for ETU then take a value which will b ealmost inbetween based on CDF plots
    t.ue_conf.dr = 1;
    t.ue_conf.searchGrid = 100;           % this parameter is determinded by the range of the scnario!!! 
    t.ue_conf.doMPS = true;
    t.ue_conf.corr_method = 'time';
    
    t.ls_conf.method = 'one';
    
    t.sink_conf.save_file = 'false';
    t.sink_conf.debug_store = 'false';
    t.sink_conf.plot_cdl = 'true';
    t.sink_conf.plot_geom = 'true';
    tc = t;

    
end
