%% Read the historic measurement
%
function readcase
reader_o = sink(0,0,0,0);
Res_rr=cell(1,7);
if exist('ResultZYJ.mat')
    load('ResultZYJ.mat');
    [Case_num_range,~]=size(Res_tot);
    Allthe_Measurement=Res_tot
    Case_num=input(['Which case want to check? (1~', num2str(Case_num_range), ') \n']);
    for ii=1:7
        Res_rr{ii}=Res_tot{Case_num,ii};
    end
end
sink_o = sink('false' , 'true' , 'default' , 'true');
fprintf('\n')
disp(Res_rr{1});
fprintf(['  >>' Res_rr{7} '<<\n \n']);
Res_rr=sink_o.analyze_result(Res_rr);
end

%}