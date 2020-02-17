%% Main function to run a set of testcases defined in tc.m
% tc  is a integer array of test cases to be run from the tc.m file

%result is a 3D-cell structure which is storing the test result sort by date 

%% clear
clear
%% Set Timer
Time=datetime('now');
clock1=tic;
Time_apply_channel=0;

%% Initialization
Conf1= tc(1);          % t,t_conf stores all the unchanged configuration
Con=Conf1{1};
Pos=struct;           % Pos stores all the configuration changed in every iteration
Res = cell(1,6);      % Time, configuration, OTOA, corr_test, (x,y,z)matrix, Error

%% Iteration 
%iteration for each Pos
% Buffer for ue_pos and est_pos
ue_pos=nan(Con.env_fix.ue_nbr,3);
est_coord=nan(Con.env_fix.ue_nbr,3);
for ue_idx = 1:Con.env_fix.ue_nbr
    %ue_idx = 1;
    fprintf('%d\n',ue_idx); %Atlease something can show up during the boring waiting time 
    
    % constructor for Pos
    Pos = environment(Con.env_fix.trp_pos , Con.env_fix.trp_bearing,...
                      Con.env_fix.max_dim , Con.env_fix.ch, Con.env_fix.ue_pos,...
                      Con.env_fix.trpselect , Con.env_fix.num_trp);
   
    % constructor for gnb ue ls object
    gnb_o = gnb_tx(Pos.trp, Con.nr);
    ue_o = ue(Con.ue_conf,Con.nr);      
    ls_o = loc_est(Con.ls_conf, Pos.trp, Con.env_fix.max_dim, Con.nr);
   
    % signal generator
    sig_gnb = gnb_o.gen_sf(Con.slotNumber);    %size 30720X32X6
                   
    ue_pos(ue_idx, :) = Pos.UTs.ArrayPosition';   %assign ue_pos by Pos.UTs.ArrayPosition
        
    % Create channel at UE location
    clock2=tic;
    sig_ch = nan(length(sig_gnb(:,1,1)),length(Pos.trp));     %sig_ch buffer
    for TRPIndex = 1:length(Pos.trp)         
        sig_ch(:,TRPIndex) = Pos.apply_ch(sig_gnb(:, :, TRPIndex), TRPIndex);      
    end
    Time_apply_channel=Time_apply_channel+toc(clock2);           % time elapse to generate the received signal
    
    % apply noise
    
    sig_gnb_one = squeeze(sig_gnb(:,1,:));      % only the signal from the first array element is used in this case       

    % Let UE decode prs and measure TOA
    [OTOA,corr_test] = ue_o.decode_prs(sig_ch, sig_gnb_one);

    % Calculate the geo-coordinates - Location server 
    est_coord(ue_idx,:) = ls_o.localize(OTOA);    
    
    end 
        
Res{1} = Time; 
Res{2} = Con;
Res{3} = OTOA;
%Res{4} = corr_test;   %The correlation from the last measurement  
%Res{5} = est_coord-ue_pos;
Res{6} = sqrt((est_coord(:,1)-ue_pos(:,1)).^2+(est_coord(:,2)-ue_pos(:,2)).^2);
               
% save tc configuration and result:      
Time_total=toc(clock1);    

%% The ratio between OTOA and the FirstPathDelay (descibe how many precentage off from measurement to the expected)
FirstPathDelay=nan(1,length(OTOA));
for ii=1:length(OTOA)
    FirstPathDelay(ii)=Pos.trp(ii).Channel.FirstPathDelay;
end
Off = OTOA./FirstPathDelay
Res{4}=Pos;
Res{5}=Off;


%% check if the result is resonable or not.
%% The correlation plot
clf
%figure(1)
%stem(abs(Res{4}(:,1)));

% The cdf
figure(1)
cdfplot(Res{6})
% Plot the environment
Pos.plotRoom(Pos.trp,Pos.UTs,1.2,2);
%}

%{              
% Check for the position generator!!! (This part will be removed)
  Pos_check=Pos.UTs.ArrayPosition;
  for ii=1:6              
      Pos_check=[Pos_check Pos.trp(ii).ArrayPosition];
  end
  Pos_check
 %}

%% add the result into the historic data file (sort by date)
%{
SaveData=input('Do you want to save these measurement? \n');
if SaveData == 1   
Displayplot=0;  %would you like to see the plot?

Res_tmp=Res;
emp=[];
if ~exist('ResultZYJ.mat')      % if the mat file doesn't exit, create a new empyty file 'ResultZYJ.mat'
    save('ResultZYJ.mat','emp');
end
clear('emp')

Lastday=date;
Lastday([3 7])='_';
Res_name=['Res_' Lastday];

load('ResultZYJ.mat',Res_name);   %if the variable is not exist in the file, then a new variable cooresponding with the date will be created.
Res_tot={};
if ~exist(Res_name)
    assignin('base',Res_name,Res_tmp);
    save('ResultZYJ.mat',Res_name,'-append');
else
    assignin('base','Res_tot',eval(Res_name)); 
    if Res_tot{end,1}~= Res_tmp{1}  %if the time from the last measurement in the file is not the same with the current result, then the current result would be added into the file
       Res_tot=[Res_tot; Res_tmp];
       assignin('base',Res_name,Res_tot);
       save('ResultZYJ.mat',Res_name,'-append');
    end
end

whos('-file','ResultZYJ.mat')

if Displayplot == 1
% plot
clf
% The cdf
figure(1)
cdfplot(Res_tot{end,6})

Trp_info=Res_tot{end,4}.trp;
UTs_info=Res_tot{end,4}.UTs;
Res_tot{end,4}.plotRoom(Trp_info,UTs_info(end),1.2,2);

OTOA=Res_tot{end,3};
FirstPathDelay=nan(1,length(OTOA));
for ii=1:length(OTOA)
    FirstPathDelay(ii)=Pos.trp(ii).Channel.FirstPathDelay;
end
Off = OTOA./FirstPathDelay
end

end
%}

%% Historic data reader








