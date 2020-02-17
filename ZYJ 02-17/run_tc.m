%% Main function to run a set of testcases defined in tc.m
% tc  is a integer array of test cases to be run from the tc.m file
function result = run_tc(tc_idx)
tic
    t = tc(tc_idx);
    if nargin == 1
        result = cell(length(t), 1);
        for i = 1:length(t)
            result{i} = run_simulator(t{i});
        end
        % save tc configuration and result:
        
        
              
        
    else
        error('No test case added or test case not indicated as integer vector');
    end
toc    
end

%% Run one testcase sequentially, for each UE signals from all TRP
function result = run_simulator(tc)
    
    % Config Sink (debug and result output) TODO
    %sink_o = sink(tc);

    % Config gNBs
    gnb_o = gnb_tx(tc.env.trp, tc.nr);
    
    % Config UEs
    ue_o = ue(tc.ue_conf,tc.nr);
    
    % Generate signals and beams from all TRPs for 1 slot
    sig_gnb = gnb_o.gen_sf(tc.slotNumber);
    
    % Setup the location calculation entity / location server (ls)
    ls_o = loc_est(tc.ls_conf, tc.env.trp, tc.env.max_dim, tc.nr);
    
    % Loop over UE, TRPs(& sectors)/cells
    est_coord = nan(tc.env.ue_nbr, 3);
    ue_pos = nan(tc.env.ue_nbr, 3);
    for ue_idx = 1:tc.env.ue_nbr
        fprintf('%d\n',ue_idx); %dddddddd
                            
        ue_pos(ue_idx, :) = tc.env.dropUE();
        
        % Create channel at UE location
        sig_ch = [];
        for TRPIndex = 1:length(tc.env.trp)
            sig_ch = [sig_ch tc.env.apply_ch(sig_gnb(:, :, TRPIndex), TRPIndex)];
        end
        
        % TRP: Change to select best beam per TRP or better, only generate 
        % the best beam. 
        % ToDo: 1 here is just test - must be updated
        sig_gnb_one = squeeze(sig_gnb(:,1,:));         

        % Let UE decode prs and measure TOA
        otoa = ue_o.decode_prs(sig_ch, sig_gnb_one);

        % Calculate the geo-coordinates - Location server 
        est_coord(ue_idx,:) = ls_o.localize(otoa);
        
    end
    
    % Analyze result TODO
    result = est_coord-ue_pos;
    %result = sink_o.analyze_result(est_coord);
    
    % Plot result graphs TODO
    %sink_o.plot();
    
    % Store results to file TODO
    %sink_o.save_tc_file();
    
end
    



