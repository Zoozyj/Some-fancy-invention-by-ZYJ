% Usage manual:   
%
% 1. Store the measurement into one file call 'ResultZYJ.mat' 
%    >>in tc, set sink_conf.save_file = 'true' or 'false'
%    >>Display 'Do you want to save blablabla...?',then type '1':save; '0':not save 
%
% 2. plot the cdf of the measurement 
%    >>in tc, set plot_cdl = 'true'  or 'false'
%
% 3. plot the geometry of the room with UEs and TRPs 
%    >>in tc, sink_conf.plot_geom = 'default'  or 'Manual' or 'false'
%      if =='default'
%         >> plot the last UE's geometry
%      if =='Manual'
%         >> Display 'Which one....', then type a number n
%            >>plot the n UE's geometry
%
% 4. Check The measured OTOA off from the real OTOA
%  
%
% 5. plot the cross-corelation (optional)
%    >>in tc, set sink_conf.plot_cox = 'true'  or 'false'
%
% 6. Read the historic measurement data from the file'ResultZYJ.mat' 
%    >>call the function CaseReader externally
%
% 7. Plot the 3D bar for a resource gird
%    >>call the function CaseReader externally
%
%
%
%
%
classdef sink < handle
    
    properties
        save_file
        plot_cdl
        plot_geom
        plot_cox
    end

    methods
      %% Constructor
      function this = sink(save_file,plot_cdl,plot_geom,plot_cox)
          if nargin == 4
              this.save_file = save_file;
              this.plot_cdl = plot_cdl;
              this.plot_geom = plot_geom;
              this.plot_cox = plot_cox;
          else
              error('Wrong number of input argument for class ue constructor')
          end
      end   
      %% All in one
      function Res=analyze_result(obj,Res)
          % 1.save file
          switch obj.save_file
              case 'true'
                  Res=obj.saveRes(Res);
                  
              case 'false'
              otherwise
                  warning('save_file >>>> Unexpected input arguement.')
          end
          % 2. plot cdl
          switch obj.plot_cdl
              case 'true'
                  obj.plotaCDF(Res{5},1);
              case 'false'
              otherwise
                  warning('plot_cdl >>>> Unexpected input arguement.')
          end
          % 3. plot room and check Off
          switch obj.plot_geom
              case 'default'
                  Pos=Res{3};Pos=Res{3};
                  obj.plotRoom(Pos{end}.trp ,Pos{end}.UTs , 1.2 , 2);
                  title(['The geometry of the indoor office, the last case, ', Res{2}.env_fix.ch]);
                  FirstPathDelay=obj.Off(Res{3}{end},Res{4}(end,:));
              case 'manual'
                  sim_num_range=length(Res{3});
                  Sim_num=input(['Choose one case: (1~' , num2str(sim_num_range) , ') \n ']);
                  obj.plotRoom(Res{3}{Sim_num}.trp , Res{3}{Sim_num}.UTs , 1.2 , 2);
                  title(['The geometry of the indoor office, case' ,num2str(sim_num_range),', ', Res{2}.env_fix.ch]);
                  FirstPathDelay=obj.Off(Res{3}{Sim_num},Res{4}(Sim_num,:));
              case 'false'   
              otherwise
                  warning('plot_geom >>>> Unexpected input arguement.')
          end
          
          % 4. check coxx (not implemented yet)
             switch obj.plot_cox
              case 'true'
         
                  Marker1=Res{4}(end,1:6) * Res{2}.ue_conf.interp_fac * Res{2}.nr.Fsamp + floor((Res{2}.ue_conf.searchGrid *2)*Res{2}.ue_conf.interp_fac/2);
                  Marker2=FirstPathDelay(1:6) * Res{2}.ue_conf.interp_fac * Res{2}.nr.Fsamp + floor((Res{2}.ue_conf.searchGrid *2)*Res{2}.ue_conf.interp_fac/2);
                  Marker2=round(Marker2);
                  Coxx(obj,Res{6},3,Marker1,Marker2);
                  
              case 'false'
              otherwise
                  warning('plot_cox >>>> Unexpected input arguement.')
          end
          
       
      end
      %% Store the measurement into one file call 'ResultZYJ.mat'
      function Res=saveRes(obj,Res)
          
          SaveData=input('Do you want to save the measurement? (Yes:1, No:0)\n');
          if SaveData == 1
              % leave some comment in the Res{7}
              Res{7} = input('leave some comment to this measurement: \n','s');
              % if you want to check that the data you just saved are correct, then turn it to 1
              Displayplot=0;
              % if the 'ResultZYJ.mat' file doesn't exit, create a new empyty file 'ResultZYJ.mat'
              Res_tot=cell(1,7);
              if ~exist('ResultZYJ.mat')
                  Res_tot{1}=datetime('now');
                  save('ResultZYJ.mat','Res_tot');
              else
                  load('ResultZYJ.mat','Res_tot');
              end
              
              % if the variable is not exist in the file, then a new variable cooresponding with the date will be created.
              if Res_tot{1,1} ~= Res{1}
                  Res_tot=[Res; Res_tot];
                  save('ResultZYJ.mat','Res_tot');
              else
                  Res_tot{1,7} = Res{7};
                  save('ResultZYJ.mat','Res_tot');
              end
              Res_tot
              
              if Displayplot == 1
                  clf
                  figure(1)
                  cdfplot(Res_tot{1,5});
                  Trp_info=Res_tot{1,3}{end}.trp;
                  UTs_info=Res_tot{1,3}{end}.UTs;
                  plotRoom(Trp_info,UTs_info,1.2,2);
                  OTOA=Res_tot{1,4};
                  FirstPathDelay=nan(1,length(OTOA(end,:)));
                  for ii=1:length(OTOA(end,:))
                      FirstPathDelay(ii)=Trp_info(ii).Channel.FirstPathDelay;
                  end
                  Off = OTOA./FirstPathDelay
              end
          end  
      end
      
      %% plot the cdf of the measurement
      function plotaCDF(obj,error,fig_num)
          if (nargin < 3)
              figure;
              clf
          else
              figure(fig_num);
              clf(fig_num);
          end
          hold on;
          grid on;
          cdfplot(error);
          xlabel('Horizontal Error/m');
          ylabel('CDF');
          title('CDF measurement');
          
      end
      
      %% plot the geometry of the room with UEs and TRPs
      function plotRoom(obj,TRPs,UTs,lineWidth,fig_num)
          if (nargin < 4)
              lineWidth = 1.2;
          end
          if (nargin < 5)
              figure;
              clf;
              hold on;
          else
              figure(fig_num);
              clf(fig_num);
              hold on;
          end
          R = 3;
          hold on;
          line([0 120 120 0 0],[0 0 50 50 0]);                           %sketch the wall of the office
          for TRPIndex = 1:length(TRPs)
              plot(TRPs(TRPIndex).ArrayPosition(1),TRPs(TRPIndex).ArrayPosition(2),'bo','LineWidth',lineWidth);
              h = line(TRPs(TRPIndex).ArrayPosition(1)+1.5*R*[0,cos(deg2rad(TRPs(TRPIndex).ArrayBearing+90))],TRPs(TRPIndex).ArrayPosition(2)+1.5*R*[0,sin(deg2rad(TRPs(TRPIndex).ArrayBearing+90 ))]);
              set(h,'Color','r');
              set(h,'LineWidth',lineWidth);
              phiGCS = 37.5 + 15*(7:-1:0);
              if (lineWidth>1)
                  for AngleIdx = 1:length(phiGCS)
                      hold on;
                      h = line(TRPs(TRPIndex).ArrayPosition(1) + R*[0,cos(deg2rad(TRPs(TRPIndex).ArrayBearing + phiGCS(AngleIdx)))],TRPs(TRPIndex).ArrayPosition(2)+R*[0,sin(deg2rad(TRPs(TRPIndex).ArrayBearing + phiGCS(AngleIdx)))]);
                      set(h,'Color','k');
                      set(h,'LineWidth',lineWidth/2);
                  end
              end
              if (nargin>3)
                  text(TRPs(TRPIndex).ArrayPosition(1)+2.5*R*cos(deg2rad(TRPs(TRPIndex).ArrayBearing+90)),...
                      TRPs(TRPIndex).ArrayPosition(2)+2.5*R*sin(deg2rad(TRPs(TRPIndex).ArrayBearing+90 )),...
                      num2str(TRPIndex));    %test the number of the trps alongside them
                  %text(TRPs(TRPIndex).ArrayPosition(1)+3,TRPs(TRPIndex).ArrayPosition(2),num2str(TRPIndex));
              end
          end
          plot(UTs(1).ArrayPosition(1,:),UTs.ArrayPosition(2,:),'gx');
          h = line(UTs.ArrayPosition(1)+R*[0,cos(deg2rad(UTs.UTDirectionOfTravel(1)))],UTs.ArrayPosition(2)+R*[0,sin(deg2rad(UTs.UTDirectionOfTravel(1)))]);
          set(h,'Color','k');
          set(h,'LineWidth',lineWidth);
          axis equal;
          grid on;
          axis([0 120 0 50]);
          xlabel('x/m');
          ylabel('y/m');
      end
      
      %% plot the cross-corelation (not implement completely)
      function Coxx(obj,corr,fig_num,Marker1,Marker2)
          if (nargin<3)
              figure;
              clf;
          else
              figure(fig_num);
              clf(fig_num);
          end
          if (nargin<4)
              Marker1=0;
          end
          if (nargin<5)
              Marker2=0;
          end
          [len_x,num_subfig]=size(corr);
          if num_subfig>6
              corr=corr(:,1:6); % if the input is larger than the size n X 6, then just take the first 6 of them.
          end
          legend()
          for subfig_idx=1:6
              subplot(3,2,subfig_idx);
              plot(1:len_x , corr(:,subfig_idx) , 'b');
              axis tight
              hold on;
              if Marker1(1)~=0
                  if (Marker1(subfig_idx)<=len_x && Marker1(subfig_idx)>=1)
                      p1=plot(Marker1(subfig_idx),corr(Marker1(subfig_idx),subfig_idx),'^','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r');
                      M1=1;
                  else
                      M1=0;
                  end
              end
              if Marker2(1)~=0
                  if (Marker2(subfig_idx)<=len_x && Marker2(subfig_idx)>=1)
                      p2=plot(Marker2(subfig_idx),corr(Marker2(subfig_idx),subfig_idx),'o','MarkerSize',5,'MarkerEdgeColor','g','MarkerFaceColor','g');
                              M2=1;
                  else
                      M2=0;
                  end
              end
              if (subfig_idx==1 && M1==1 && M2==1 )
                  legend([p1 p2],{'Measure','Real'});
                 legend('boxoff');
              end
          
              grid on;
              xlabel('Sample');
              axis([1 len_x 0 max(corr(:,subfig_idx)*1.1)]);
              title(['Corr of Trp' num2str(subfig_idx)]);
          end
        
      end
  
      %% The ratio between OTOA and the FirstPathDelay (descibe how many precentage off from measurement to the expected)
      function  FirstPathDelay=Off(obj,Pos,OTOA)
          FirstPathDelay=nan(1,length(OTOA(end,:)));
          for ii=1:length(OTOA(end,:)) 
              FirstPathDelay(ii)=Pos.trp(ii).Channel.FirstPathDelay;
          end
          disp(['Off  = ' , num2str(OTOA(end,:)./FirstPathDelay) ]);
      end
      
      %% case Reader
      function Res=CaseReader(obj)
          Res=cell(1,7);
          if exist('ResultZYJ.mat')
              load('ResultZYJ.mat');
              [Case_num_range,~]=size(Res_tot);
              Res_tot
              Case_num=input(['Which case want to check? (1~', num2str(Case_num_range), ') \n']);
              for ii=1:7
                Res{ii}=Res_tot{Case_num,ii};
              end
          end          
      end
      
      %% Draw the resource grid
      function RBPlot(obj,sig_barplot)
          t=1:14;
          colorplate=['r', 'g', 'b', 'c', 'm', 'y', 'k','w'];
          figure(4)
          clf(4)
          for ii=1:6
              bar3(4:11,abs(squeeze(sig_barplot(:,4:11,ii)))',colorplate(ii));
              axis equal
              hold on
          end
          h=bar3([1:3 12:14],abs(squeeze(sig_barplot(:,[1:3 12:14],ii)))','w');
          title('Resource grid');
          xlabel('Subcarriers');
          ylabel('Symbol Number');
      end
      
    end
end
