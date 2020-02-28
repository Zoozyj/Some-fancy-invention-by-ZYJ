%% Main function to run a set of testcases defined in tc.m
% tc  is a integer array of test cases to be run from the tc.m file

%result is a 3D-cell structure which is storing the test result sort by date

%% clear
clear
% Set Timer
Time=datetime('now');
clock1=tic;
Time_apply_channel=0;

% Initialization
Conf1= tc(1);          % t,t_conf stores all the unchanged configuration
Con=Conf1{1};
Pos=struct;           % Pos stores all the configuration changed in every iteration
Res = cell(1,7);      % Time, configuration, OTOA, corr_test, (x,y,z)matrix, Error

% Iteration
% iteration for each Pos
% Buffer for ue_pos and est_pos
ue_pos=nan(Con.env.ue_nbr,3);
est_coord=nan(Con.env.ue_nbr,3);
OTOA=nan(Con.env.ue_nbr,Con.env.num_trp);
Pos=cell(Con.env.ue_nbr,1);
% constructor for Environment and Ue
env_o = environment(Con.env.trp_pos , Con.env.trp_bearing,...
    Con.env.max_dim , Con.env.ch , Con.env.ue_pos,...
    Con.env.trpselect , Con.env.num_trp , Con.nr ,...
    Con.env.Sig_int_end-Con.env.Sig_int_star+1);
ue_o = ue(Con.ue_conf,Con.nr);


for ue_idx = 1:Con.env.ue_nbr
    %ue_idx = 1;
    fprintf('%d \n',ue_idx); %Atlease something can show up during the boring waiting time
    
    Pos_tmp=env_o.setenvironment;
    Pos_tmp.Num=num2str(ue_idx);
    
    % constructor for gnb ue ls object
    gnb_o = gnb_tx(Pos_tmp.trp, Con.nr);
    ls_o = ls(Con.ls_conf, Pos_tmp.trp, Con.env.max_dim, Con.nr);
    
    % signal generator
    [sig_gnb,~]= gnb_o.gen_sf(Con.slotNumber);    %size 30720X32X6
    
    %extract the ue position
    ue_pos(ue_idx, :) = Pos_tmp.UTs.ArrayPosition';   %assign ue_pos by Pos.UTs.ArrayPosition
    
    % Create channel at UE location
    clock2=tic;
    sig_ch = zeros(length(sig_gnb(:,1,1)),length(Pos_tmp.trp));     %sig_ch buffer
    for TRPIndex = 1:length(Pos_tmp.trp)
        sig_ch(Con.env.Sig_int_star:Con.env.Sig_int_end,TRPIndex) ...
        = env_o.apply_ch(sig_gnb(Con.env.Sig_int_star:Con.env.Sig_int_end, :, TRPIndex) , Pos_tmp.trp(TRPIndex) , Pos_tmp.UTs);  %only consider the sample range that is interested
    end
    Time_apply_channel=Time_apply_channel+toc(clock2);           % time elapse to generate the received signal
    
    % only the signal from the first array element is used in this case
    sig_gnb_one = squeeze(sig_gnb(:,1,:));
    
    % Let UE decode prs and measure TOA
    [OTOA(ue_idx,:),corr_test] = ue_o.decode_prs(sig_ch, sig_gnb_one);
    
    OTOA_tmp=ls_o.OTOASelect(OTOA(ue_idx,:));
    
    % Calculate the geo-coordinates - Location server
    est_coord(ue_idx,:) = ls_o.localize(OTOA_tmp);
    
    % save the case info
    Pos{ue_idx}=Pos_tmp;
end
fprintf('\n Average execution time for apply_ch: %2.2f s \n',Time_apply_channel/length(Pos_tmp.trp)/Con.env.ue_nbr);

% save tc configuration and result
Res{1} = Time;
Res{2} = Con;
Res{3} = Pos;
Res{4} = OTOA;
Res{5} = sqrt((est_coord(:,1)-ue_pos(:,1)).^2+(est_coord(:,2)-ue_pos(:,2)).^2);
Res{6} = corr_test;
Res{7} = [];       % leave some comment in here
Time_total=toc(clock1);

%% result analyse
%Constructor
Conf1= tc(1);          % t,t_conf stores all the unchanged configuration
Con=Conf1{1};
sink_o = sink(Con.sink_conf.save_file , Con.sink_conf.plot_cdl ,...
              Con.sink_conf.plot_geom , Con.sink_conf.plot_cox);

Res=sink_o.analyze_result(Res);

















