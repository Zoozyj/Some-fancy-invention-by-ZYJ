%% NR specific static properties and NR configuration parameters shall be defined here. 
% One NR object can be created for each test case configured to run.
% Note, class is a handle class and may be passed as reference without
% being copied.

classdef nr < handle
    
    % Public NR basic constants
    properties (Constant)
        N_fft = 2048                           % NR Tx basic FFT size    ddddd
        N_ofdm = 14                            % Number OFDM symbols in slot
        N_sf = 10                              % Number of subframes in frame
        N_s_slot = nr.N_fft*(nr.N_ofdm+1)     % Samples per slot
        N_s_cp_ofdm = [160, 144, 144, 144, 144, 144, 144]*(nr.N_fft/2048) % Cyclic prefix per ofdm symbol
    end
    
    % NR basic constants that not is exposed to simulator - used in this
    % object/class only
    properties(Constant, Access = private) 
        u_v = 0:4                              % Index for ofdm versions
        Tc_sc_v = 1./(15e3 * 2.^nr.u_v * nr.N_fft)% Basic NR time unit per subcarrier spacing
        Fsamp_v = 1./nr.Tc_sc_v                % Basic sampling frequency per subcarrier spacing
        Fsc_v = 2.^nr.u_v*15e3                % Subcarrier spacings in Hz
        Tu_v = 1./nr.Fsc_v                     % OFDM Symbol time in s
        T_slot_v = 1e-3*2.^-nr.u_v            % Time for one slot
        N_slots_v =  2.^nr.u_v                 % Slots in one subframe     dddd
        T_cp_v = (nr.T_slot_v - nr.Tu_v .* nr.N_ofdm)./nr.N_ofdm % Time of cyclic prefix        
    end
    
    properties 
        % Basic public NR properties dependent on object input / numerology
        Tc       % Basic time unit
        Fsamp    % Sampling frequency
        Fsc      % Subcarrier spacing
        Tu       % OFDM symbol time
        T_slot   % Time of slot 
        N_slots  % Number slots per subframe
        N_s_sf   % Number sample per subframe
        T_cp     % Time of cyclic prefix 
        MaxDLNumberRB   % Max nbr of resource blocks of the NR cells in simulation
        DLNumberRB      % Actual used nbr of resource blocks for the NR cells in simulation
        
        % Positioning specifics
        prs_seq         % {0: Gold, 1: Zaduff-Chu}
        prs_comb_n      % {1, 2, 4, 8, 12}
        prs_symbols     % list of symbol indexes with prs ex. [3:11]
        PRSNumberRB     % Nbr resource blocks allocated for PRS 
    end
    
    methods
        
        % Constructor
        function obj = nr(u, MaxDLNumber, DLNumberRB, PRSNumberRB, prs_symbols)
            if nargin == 1
                u_idx = u;                                             % Zhang: No !!??? if I want to set u_v=0, should I set u=-1 here?
                obj.Tc = obj.Tc_sc_v(u_idx);                  
                obj.Fsamp = obj.Fsamp_v(u_idx);
                obj.Fsc = obj.Fsc_v(u_idx);
                obj.Tu = obj.Tu_v(u_idx);
                obj.T_slot = obj.T_slot_v(u_idx);
                obj.N_slots = obj.N_slots_v(u_idx);
                obj.T_cp = obj.T_cp_v(u_idx);
                obj.MaxDLNumberRB = 110;
                obj.DLNumberRB = 110;
                obj.PRSNumberRB = 69;
                obj.prs_symbols = [3:11];
            elseif nargin == 5
                u_idx = u;                                             %Zhang: The same as above
                obj.N_s_sf = obj.N_s_slot*obj.N_slots
                obj.Tc = obj.Tc_sc_v(u_idx);
                obj.Fsamp = obj.Fsamp_v(u_idx);
                obj.Fsc = obj.Fsc_v(u_idx);
                obj.Tu = obj.Tu_v(u_idx);
                obj.T_slot = obj.T_slot_v(u_idx);
                obj.N_slots = obj.N_slots_v(u_idx);
                obj.T_cp = obj.T_cp_v(u_idx);
                obj.MaxDLNumberRB = MaxDLNumber;
                obj.DLNumberRB = DLNumberRB;
                obj.PRSNumberRB = PRSNumberRB;
                obj.prs_symbols = prs_symbols;
            else
                error('Wrong number of arguments for NR definition.')
            end
        end
        
    end
end
    