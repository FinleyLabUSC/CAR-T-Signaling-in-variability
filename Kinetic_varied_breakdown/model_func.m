function [err, timepoints, species_out, observables_out] = model_func( timepoints, parameters, species_init, suppress_plot )
%SUPPLEMENTAL FILE S3 Integrate reaction network and plot observables.
%   Integrates the reaction network corresponding to the BioNetGen model
%   'Supplemental File S3' and then (optionally) plots the observable trajectories,
%   or species trajectories if no observables are defined. Trajectories are
%   generated using either default or user-defined parameters and initial
%   species values. Integration is performed by the MATLAB stiff solver
%   'ode15s'. SUPPLEMENTAL FILE S3 returns an error value, a vector of timepoints,
%   species trajectories, and observable trajectories.
%   
%   [err, timepoints, species_out, observables_out]
%        = Supplemental File S3( timepoints, species_init, parameters, suppress_plot )
%
%   INPUTS:
%   -------
%   species_init    : row vector of 225 initial species populations.
%   timepoints      : column vector of time points returned by integrator.
%   parameters      : row vector of 135 model parameters.
%   suppress_plot   : 0 if a plot is desired (default), 1 if plot is suppressed.
%
%   Note: to specify default value for an input argument, pass the empty array.
%
%   OUTPUTS:
%   --------
%   err             : 0 if the integrator exits without error, non-zero otherwise.
%   timepoints      : a row vector of timepoints returned by the integrator.
%   species_out     : array of species population trajectories
%                        (columns correspond to species, rows correspond to time).
%   observables_out : array of observable trajectories
%                        (columns correspond to observables, rows correspond to time).
%
%   QUESTIONS about the BNG Mfile generator?  Email justinshogg@gmail.com



%% Process input arguments

% define any missing arguments
if ( nargin < 1 )
    timepoints = [];
end

if ( nargin < 3 )
    species_init = [];
end

if ( nargin < 2 )
    parameters = [];
end

if ( nargin < 4 )
    suppress_plot = 1;
end


% initialize outputs (to avoid error msgs if script terminates early
err = 0;
species_out     = [];
observables_out = [];


% setup default parameters, if necessary
if ( isempty(parameters) )
   parameters = [ 30, 1000, 250, 500, 500, 870, 540, 600, 300, 120, 25, 800, 60, 600, 150, 150, 250, 20, 350, 1500, 1250, 0, 400, 1000, .1, 0.023104906, 1000, 50, 40.26185694, 4.852880242, 6.373803376, .1, 416.8324872, 30.89383828, 137.0133105, 0.233721543, 0.1, 2.026709313, 0.391443224, 1.310391485, 0.101597, 96.548, 331.68, 270, 153.82, 484.9, 405.43, 0.1153, 360, 0.1, 500, 101, 28, 96, 12.8, .1, 20, 0.05, 200, 200, 20, 0.1, 10, 0.1, 75, 0.1, 142, .1, 110, .1, 15, 0.1, 55, 0.1, 250, 0.1, 65, .1, 16333, .1, 2400, 0.01, 9666, 0.03, .1, 1266, 2.28, .1, 9666, .18, .1, 66, 6, 0.1, 186, 10, .1, 2000, .6, 12, 33, 200, 33, 60, 86, 300, 33, 10000, 66, 10000, 66, 53, 60, 33, 1200, 33, 15200, 7, 2000, 12, 66, 66, 1, 0.0015, 100, 100, 100, 10, 1.248432393, 0.094342447, 0.963404796, 67.84415559, 1.459578998, 9.918819462, 44.15493881 ];
end
% check that parameters has proper dimensions
if (  size(parameters,1) ~= 1  ||  size(parameters,2) ~= 135  )
    fprintf( 1, 'Error: size of parameter argument is invalid! Correct size = [1 135].\n' );
    err = 1;
    return;
end

% setup default initial values, if necessary
if ( isempty(species_init) )
   species_init = initialize_species( parameters );
end
% check that species_init has proper dimensions
if (  size(species_init,1) ~= 1  ||  size(species_init,2) ~= 225  )
    fprintf( 1, 'Error: size of species_init argument is invalid! Correct size = [1 225].\n' );
    err = 1;
    return;
end

% setup default timepoints, if necessary
if ( isempty(timepoints) )
   timepoints = linspace(0,10,20+1)';
end
% check that timepoints has proper dimensions
if (  size(timepoints,1) < 2  ||  size(timepoints,2) ~= 1  )
    fprintf( 1, 'Error: size of timepoints argument is invalid! Correct size = [t 1], t>1.\n' );
    err = 1;
    return;
end

% setup default suppress_plot, if necessary
if ( isempty(suppress_plot) )
   suppress_plot = 0;
end
% check that suppress_plot has proper dimensions
if ( size(suppress_plot,1) ~= 1  ||  size(suppress_plot,2) ~= 1 )
    fprintf( 1, 'Error: suppress_plots argument should be a scalar!\n' );
    err = 1;
    return;
end

% define parameter labels (this is for the user's reference!)
param_labels = { 'Ant_T', 'CD45_T', 'LCK_T', 'CD3z_T', 'CD28_T', 'LAT_T', 'Gads_T',...
                 'ZAP70_T', 'SLP76_T', 'p85a_T', 'CSK_T', 'Grb2_T', 'SOS_T', 'Ras_T',...
                 'RAF_T', 'MEK_T', 'ERK_T', 'RasGAP_T', 'PLCg_T', 'RasGRP_T', 'PIP2_T',...
                 'DAG_T', 'SHP1_T', 'ant_Kd', 'ant_on', 'ktrans', 'km', 'xx', 'cat_UU_394',...
                 'cat_PU_505', 'cat_UP_394', 'CSKon', 'CSKoff_PU', 'CSKoff_UU', 'cat_CSK_PU',...
                 'cat_CSK_UU', 'on', 'off_PU_PP', 'off_UU_PP', 'off_PU_UP', 'off_UU_UP', 'KmA1',...
                 'KmA2', 'KmB1', 'KmB2', 'KmC1', 'KmC2', 'Xi', 'Kcat_LCKPU_CD3z', 'on_CD28', 'off_CD28',...
                 'cat191_LCKPU', 'cat206_LCKPU', 'cat209_LCKPU', 'cat218_LCKPU', 'ZAPu_on',...
                 'ZAPu_KD', 'ZAPp_on', 'ZAPp_KD', 'Kcat_LCKPU_ZAP315', 'Kcat_ZAP', 'CD28_p85a_on',...
                 'CD28_p85a_KD', 'CD28_GADS_on', 'CD28_GADS_KD', 'CD28_Grb2_on', 'CD28_Grb2_KD',...
                 'LAT_GADS_on', 'LAT_GADS_KD', 'GADS_SLP_on', 'GADS_SLP_KD', 'LAT_Grb2_on', 'LAT_Grb2_KD',...
                 'Grb2free_SOS_on', 'Grb2free_SOS_KD', 'Grb2bound_SOS_on', 'Grb2bound_SOS_KD',...
                 'on_SOSallo_GDP', 'KD_SOSallo_GDP', 'on_SOSallo_GTP', 'KD_SOSallo_GTP', 'on_SOScat_none',...
                 'KD_SOScat_none', 'cat_SOScat_none', 'on_SOScat_GTP', 'KD_SOScat_GTP', 'cat_SOScat_GTP',...
                 'on_SOScat_GDP', 'KD_SOScat_GDP', 'cat_SOScat_GDP', 'on_RasGAP', 'KD_RasGAP', 'cat_RasGAP',...
                 'on_PLC_LAT', 'KD_PLC_LAT', 'Kcat_PLCg', 'on_RasGRP', 'KD_RasGRP', 'Kcat_RasGRP',...
                 'Kcat1', 'Km1', 'Vmax2', 'Km2', 'Kcat3', 'Km3', 'Kcat4', 'Km4', 'Vmax5', 'Km5',...
                 'Vmax6', 'Km6', 'Ki1', 'Kcat7', 'Km7', 'Kcat8', 'Km8', 'Vmax9', 'Km9', 'Vmax10',...
                 'Km10', 'Ki2', 'Ka', 'Fa', 'SHP1_on', 'SHP1_KD', 'kcat_LCKpu_SHP1', 'kcat_SHP1',...
                 'kcat_SHP1_SHP1', 'Kcat_CD45_dephos', 'Kcat_CD45_LCK394', 'Kcat_CD45_LCK505',...
                 'Kcat_CD45_A1', 'Kcat_CD45_ZAP315', 'Kcat_CD45_ZAP493', 'Kcat_CD45_CD28191' };



%% Integrate Network Model
 
% calculate expressions
[expressions] = calc_expressions( parameters );

% set ODE integrator options
opts = odeset( 'RelTol',   1e-8,   ...11111
               'AbsTol',   0.0001,   ...
               'Stats',    'off',  ...
               'BDF',      'off',    ...
               'MaxOrder', 5   );


% define derivative function
rhs_fcn = @(t,y)( calc_species_deriv( t, y, expressions ) );

% simulate model system (stiff integrator)

[~, species_out] = ode15s( rhs_fcn, timepoints, species_init', opts );

% calculate observables
observables_out = zeros( length(timepoints), 137 );
for t = 1 : length(timepoints)
    observables_out(t,:) = calc_observables( species_out(t,:), expressions );
end


%% Plot Output, if desired

if ( ~suppress_plot )
    
    % define plot labels
    observable_labels = { 'pY_CD3z', 'pY_CD28', 'pY171_LAT', 'pY145_SLP76', 'pY394_LCK', ...
                          'pY505_LCK', 'pY315_ZAP70', 'pY493_ZAP70', 'boundPI3K', 'boundGADS_SH2', ...
                          'ActiveRAF', 'ActiveMEK', 'ActiveERK', 'LCK_U394_U505', 'LCK_P394_P505', ...
                          'LCK_P394_U505', 'LCK_U394_P505', 'UU505', 'P_A1all', 'P_A2all', 'P_B1all', ...
                          'P_B2all', 'P_C1all', 'P_C2all', 'SHP1_active', 'P_A1', 'P_A2', 'P_B1', ...
                          'P_B2', 'P_C1', 'P_C2', 'U_A1', 'U_A2', 'U_B1', 'U_B2', 'U_C1', 'U_C2', ...
                          'LCK_U394', 'LCK_U505', 'LCK_P394', 'LCK_P505', 'CD28_U191', 'CD28_U206', ...
                          'CD28_U209', 'CD28_U218', 'CD28_P191', 'CD28_P206', 'CD28_P209', ...
                          'CD28_P218', 'ZAP70_B_U315', 'ZAP70_B_U493', 'ZAP70_ub_P315', 'ZAP70_ub_P493', ...
                          'CD45', 'UU', 'PP', 'PU', 'UP', 'PU_PP', 'UU_PP', 'PU_UP', 'UU_UP', 'LCK_UU', ...
                          'LCK_PP', 'LCK_PU', 'LCK_UP', 'LCK_B_U394_U505', 'LCK_B_P394_U505', ...
                          'LCK_B_U394_P505', 'LCK_B_P394_P505', 'PU_CSK', 'PP_CSK', 'CSK', ...
                          'CSKub', 'LAT_U171', 'LAT_P171', 'LAT_P171total', 'LAT_U1', 'LAT_P1', ...
                          'LAT_P1total', 'ZAP70_ub_U315_U493', 'ZAP70_ub_P315_U493', ...
                          'ZAP70_ub_U315_P493', 'ZAP70_ub_P315_P493', 'ZAP70_B_U315_U493', ...
                          'ZAP70_B_P315_U493', 'ZAP70_B_U315_P493', 'ZAP70_B_P315_P493', ...
                          'ZAP70_B_P315', 'ZAP70_B_P493', 'SLP76_ub_U145', 'SLP76_ub_P145', ...
                          'SLP76_B_U145', 'SLP76_B_P145', 'SLP76_B_P145total', 'p85a_ub', ...
                          'p85a_B', 'boundSHP1uA1', 'boundSHP1uA2', 'boundSHP1uB1', ...
                          'boundSHP1uB2', 'boundSHP1uC1', 'boundSHP1uC2', 'SHP1_p', 'LCK_SHPu_P', ...
                          'LCK_unprotected', 'ZAP_315', 'ZAP_493', 'P_CD28_206', 'P_CD28_209', ...
                          'PIP2conc', 'PLCg_bound', 'DAG_total', 'RasGRP_bound', 'RasGDPfree', ...
                          'RasGTPfree', 'RasGDPbound', 'RasGTPbound', 'boundGrb2', 'SOS_free', ...
                          'SOS_bound_free', 'SOS_bound_inact', 'SOS_bound_act', 'Ras_GDP', ...
                          'Ras_GTP', 'RAF_u', 'RAF_p', 'MEK_u', 'MEK_p', 'MEK_pp', 'ERK_u', ...
                          'ERK_p', 'ERK_pp', 'CD3z_AntBoundA1', 'CD3z_AntUnboundA1', ...
                          'CD3z_ITAMpartial', 'CD3z_ITAMfull' };

    % construct figure
    plot(timepoints,observables_out);
    title('Supplemental File S3 observables','fontSize',14,'Interpreter','none');
    axis([0 timepoints(end) 0 inf]);
    legend(observable_labels,'fontSize',10,'Interpreter','none');
    xlabel('time','fontSize',12,'Interpreter','none');
    ylabel('number or concentration','fontSize',12,'Interpreter','none');

end


%~~~~~~~~~~~~~~~~~~~~~%
% END of main script! %
%~~~~~~~~~~~~~~~~~~~~~%

% Define if function to allow nested if statements in user-defined functions
function [val] = if__fun (cond, valT, valF)
% IF__FUN Select between two possible return values depending on the boolean
% variable COND.
    if (cond)
        val = valT;
    else
        val = valF;
    end
end

% initialize species function
function [species_init] = initialize_species( params )

    species_init = zeros(1,225);
    species_init(1) = 6*params(1);
    species_init(2) = params(4);
    species_init(3) = params(4);
    species_init(4) = params(4);
    species_init(5) = params(3)*0;
    species_init(6) = params(3);
    species_init(7) = params(3)*0;
    species_init(8) = params(11);
    species_init(9) = params(5);
    species_init(10) = params(5);
    species_init(11) = params(5);
    species_init(12) = params(6);
    species_init(13) = params(8);
    species_init(14) = params(9);
    species_init(15) = params(7);
    species_init(16) = params(10);
    species_init(17) = params(2);
    species_init(18) = params(12);
    species_init(19) = params(13);
    species_init(20) = params(18);
    species_init(21) = params(14);
    species_init(22) = params(15);
    species_init(23) = params(16);
    species_init(24) = params(17);
    species_init(25) = params(19);
    species_init(26) = params(20);
    species_init(27) = params(21);
    species_init(28) = params(22);
    species_init(29) = params(23);
    species_init(30) = 0;
    species_init(31) = 0;
    species_init(32) = 0;
    species_init(33) = 0;
    species_init(34) = 0;
    species_init(35) = 0;
    species_init(36) = 0;
    species_init(37) = 0;
    species_init(38) = 0;
    species_init(39) = 0;
    species_init(40) = 0;
    species_init(41) = 0;
    species_init(42) = 0;
    species_init(43) = 0;
    species_init(44) = 0;
    species_init(45) = 0;
    species_init(46) = 0;
    species_init(47) = 0;
    species_init(48) = 0;
    species_init(49) = 0;
    species_init(50) = 0;
    species_init(51) = 0;
    species_init(52) = 0;
    species_init(53) = 0;
    species_init(54) = 0;
    species_init(55) = 0;
    species_init(56) = 0;
    species_init(57) = 0;
    species_init(58) = 0;
    species_init(59) = 0;
    species_init(60) = 0;
    species_init(61) = 0;
    species_init(62) = 0;
    species_init(63) = 0;
    species_init(64) = 0;
    species_init(65) = 0;
    species_init(66) = 0;
    species_init(67) = 0;
    species_init(68) = 0;
    species_init(69) = 0;
    species_init(70) = 0;
    species_init(71) = 0;
    species_init(72) = 0;
    species_init(73) = 0;
    species_init(74) = 0;
    species_init(75) = 0;
    species_init(76) = 0;
    species_init(77) = 0;
    species_init(78) = 0;
    species_init(79) = 0;
    species_init(80) = 0;
    species_init(81) = 0;
    species_init(82) = 0;
    species_init(83) = 0;
    species_init(84) = 0;
    species_init(85) = 0;
    species_init(86) = 0;
    species_init(87) = 0;
    species_init(88) = 0;
    species_init(89) = 0;
    species_init(90) = 0;
    species_init(91) = 0;
    species_init(92) = 0;
    species_init(93) = 0;
    species_init(94) = 0;
    species_init(95) = 0;
    species_init(96) = 0;
    species_init(97) = 0;
    species_init(98) = 0;
    species_init(99) = 0;
    species_init(100) = 0;
    species_init(101) = 0;
    species_init(102) = 0;
    species_init(103) = 0;
    species_init(104) = 0;
    species_init(105) = 0;
    species_init(106) = 0;
    species_init(107) = 0;
    species_init(108) = 0;
    species_init(109) = 0;
    species_init(110) = 0;
    species_init(111) = 0;
    species_init(112) = 0;
    species_init(113) = 0;
    species_init(114) = 0;
    species_init(115) = 0;
    species_init(116) = 0;
    species_init(117) = 0;
    species_init(118) = 0;
    species_init(119) = 0;
    species_init(120) = 0;
    species_init(121) = 0;
    species_init(122) = 0;
    species_init(123) = 0;
    species_init(124) = 0;
    species_init(125) = 0;
    species_init(126) = 0;
    species_init(127) = 0;
    species_init(128) = 0;
    species_init(129) = 0;
    species_init(130) = 0;
    species_init(131) = 0;
    species_init(132) = 0;
    species_init(133) = 0;
    species_init(134) = 0;
    species_init(135) = 0;
    species_init(136) = 0;
    species_init(137) = 0;
    species_init(138) = 0;
    species_init(139) = 0;
    species_init(140) = 0;
    species_init(141) = 0;
    species_init(142) = 0;
    species_init(143) = 0;
    species_init(144) = 0;
    species_init(145) = 0;
    species_init(146) = 0;
    species_init(147) = 0;
    species_init(148) = 0;
    species_init(149) = 0;
    species_init(150) = 0;
    species_init(151) = 0;
    species_init(152) = 0;
    species_init(153) = 0;
    species_init(154) = 0;
    species_init(155) = 0;
    species_init(156) = 0;
    species_init(157) = 0;
    species_init(158) = 0;
    species_init(159) = 0;
    species_init(160) = 0;
    species_init(161) = 0;
    species_init(162) = 0;
    species_init(163) = 0;
    species_init(164) = 0;
    species_init(165) = 0;
    species_init(166) = 0;
    species_init(167) = 0;
    species_init(168) = 0;
    species_init(169) = 0;
    species_init(170) = 0;
    species_init(171) = 0;
    species_init(172) = 0;
    species_init(173) = 0;
    species_init(174) = 0;
    species_init(175) = 0;
    species_init(176) = 0;
    species_init(177) = 0;
    species_init(178) = 0;
    species_init(179) = 0;
    species_init(180) = 0;
    species_init(181) = 0;
    species_init(182) = 0;
    species_init(183) = 0;
    species_init(184) = 0;
    species_init(185) = 0;
    species_init(186) = 0;
    species_init(187) = 0;
    species_init(188) = 0;
    species_init(189) = 0;
    species_init(190) = 0;
    species_init(191) = 0;
    species_init(192) = 0;
    species_init(193) = 0;
    species_init(194) = 0;
    species_init(195) = 0;
    species_init(196) = 0;
    species_init(197) = 0;
    species_init(198) = 0;
    species_init(199) = 0;
    species_init(200) = 0;
    species_init(201) = 0;
    species_init(202) = 0;
    species_init(203) = 0;
    species_init(204) = 0;
    species_init(205) = 0;
    species_init(206) = 0;
    species_init(207) = 0;
    species_init(208) = 0;
    species_init(209) = 0;
    species_init(210) = 0;
    species_init(211) = 0;
    species_init(212) = 0;
    species_init(213) = 0;
    species_init(214) = 0;
    species_init(215) = 0;
    species_init(216) = 0;
    species_init(217) = 0;
    species_init(218) = 0;
    species_init(219) = 0;
    species_init(220) = 0;
    species_init(221) = 0;
    species_init(222) = 0;
    species_init(223) = 0;
    species_init(224) = 0;
    species_init(225) = 0;

end


% user-defined functions
% function rateLaw1
function [val] = rateLaw1(expressions, observables)
    val = (observables(134)*expressions(28));
end

% function rateLaw2
function [val] = rateLaw2(expressions, observables)
    val = (((expressions(31)*observables(55))/(expressions(39)+observables(55)))/2);
end

% function rateLaw3
function [val] = rateLaw3(expressions, observables)
    val = (((expressions(32)*observables(55))/(expressions(39)+observables(57)))/2);
end

% function rateLaw4
function [val] = rateLaw4(expressions, observables)
    val = (((expressions(33)*observables(55))/(expressions(39)+observables(58)))/2);
end

% function rateLaw5
function [val] = rateLaw5(expressions, observables)
    val = ((expressions(31)*observables(57))/(expressions(40)+observables(55)));
end

% function rateLaw6
function [val] = rateLaw6(expressions, observables)
    val = ((expressions(32)*observables(57))/(expressions(40)+observables(57)));
end

% function rateLaw7
function [val] = rateLaw7(expressions, observables)
    val = ((expressions(33)*observables(57))/(expressions(40)+observables(58)));
end

% function rateLaw8
function [val] = rateLaw8(expressions, observables)
    val = (((expressions(31)*observables(56))/(expressions(41)+observables(55)))/2);
end

% function rateLaw9
function [val] = rateLaw9(expressions, observables)
    val = (((expressions(32)*observables(56))/(expressions(41)+observables(57)))/2);
end

% function rateLaw10
function [val] = rateLaw10(expressions, observables)
    val = (((expressions(33)*observables(56))/(expressions(41)+observables(58)))/2);
end

% function rateLaw11
function [val] = rateLaw11(expressions, observables)
    val = ((expressions(60)/expressions(47))/((((((((((((1+(observables(32)/expressions(47)))+(observables(33)/expressions(48)))+(observables(34)/expressions(49)))+(observables(35)/expressions(50)))+(observables(36)/expressions(51)))+(observables(37)/expressions(52)))+(observables(26)/expressions(54)))+(observables(27)/expressions(55)))+(observables(28)/expressions(56)))+(observables(29)/expressions(57)))+(observables(30)/expressions(58)))+(observables(31)/expressions(59))));
end

% function rateLaw12
function [val] = rateLaw12(expressions, observables)
    val = ((expressions(60)/expressions(48))/((((((((((((1+(observables(32)/expressions(47)))+(observables(33)/expressions(48)))+(observables(34)/expressions(49)))+(observables(35)/expressions(50)))+(observables(36)/expressions(51)))+(observables(37)/expressions(52)))+(observables(26)/expressions(54)))+(observables(27)/expressions(55)))+(observables(28)/expressions(56)))+(observables(29)/expressions(57)))+(observables(30)/expressions(58)))+(observables(31)/expressions(59))));
end

% function rateLaw13
function [val] = rateLaw13(expressions, observables)
    val = ((expressions(60)/expressions(49))/((((((((((((1+(observables(32)/expressions(47)))+(observables(33)/expressions(48)))+(observables(34)/expressions(49)))+(observables(35)/expressions(50)))+(observables(36)/expressions(51)))+(observables(37)/expressions(52)))+(observables(26)/expressions(54)))+(observables(27)/expressions(55)))+(observables(28)/expressions(56)))+(observables(29)/expressions(57)))+(observables(30)/expressions(58)))+(observables(31)/expressions(59))));
end

% function rateLaw14
function [val] = rateLaw14(expressions, observables)
    val = ((expressions(60)/expressions(50))/((((((((((((1+(observables(32)/expressions(47)))+(observables(33)/expressions(48)))+(observables(34)/expressions(49)))+(observables(35)/expressions(50)))+(observables(36)/expressions(51)))+(observables(37)/expressions(52)))+(observables(26)/expressions(54)))+(observables(27)/expressions(55)))+(observables(28)/expressions(56)))+(observables(29)/expressions(57)))+(observables(30)/expressions(58)))+(observables(31)/expressions(59))));
end

% function rateLaw15
function [val] = rateLaw15(expressions, observables)
    val = ((expressions(60)/expressions(51))/((((((((((((1+(observables(32)/expressions(47)))+(observables(33)/expressions(48)))+(observables(34)/expressions(49)))+(observables(35)/expressions(50)))+(observables(36)/expressions(51)))+(observables(37)/expressions(52)))+(observables(26)/expressions(54)))+(observables(27)/expressions(55)))+(observables(28)/expressions(56)))+(observables(29)/expressions(57)))+(observables(30)/expressions(58)))+(observables(31)/expressions(59))));
end

% function rateLaw16
function [val] = rateLaw16(expressions, observables)
    val = ((expressions(60)/expressions(52))/((((((((((((1+(observables(32)/expressions(47)))+(observables(33)/expressions(48)))+(observables(34)/expressions(49)))+(observables(35)/expressions(50)))+(observables(36)/expressions(51)))+(observables(37)/expressions(52)))+(observables(26)/expressions(54)))+(observables(27)/expressions(55)))+(observables(28)/expressions(56)))+(observables(29)/expressions(57)))+(observables(30)/expressions(58)))+(observables(31)/expressions(59))));
end

% function rateLaw17
function [val] = rateLaw17(expressions, observables)
    val = ((expressions(61)/expressions(47))/((((((((((((1+(observables(32)/expressions(47)))+(observables(33)/expressions(48)))+(observables(34)/expressions(49)))+(observables(35)/expressions(50)))+(observables(36)/expressions(51)))+(observables(37)/expressions(52)))+(observables(26)/expressions(54)))+(observables(27)/expressions(55)))+(observables(28)/expressions(56)))+(observables(29)/expressions(57)))+(observables(30)/expressions(58)))+(observables(31)/expressions(59))));
end

% function rateLaw18
function [val] = rateLaw18(expressions, observables)
    val = ((expressions(61)/expressions(48))/((((((((((((1+(observables(32)/expressions(47)))+(observables(33)/expressions(48)))+(observables(34)/expressions(49)))+(observables(35)/expressions(50)))+(observables(36)/expressions(51)))+(observables(37)/expressions(52)))+(observables(26)/expressions(54)))+(observables(27)/expressions(55)))+(observables(28)/expressions(56)))+(observables(29)/expressions(57)))+(observables(30)/expressions(58)))+(observables(31)/expressions(59))));
end

% function rateLaw19
function [val] = rateLaw19(expressions, observables)
    val = ((expressions(61)/expressions(49))/((((((((((((1+(observables(32)/expressions(47)))+(observables(33)/expressions(48)))+(observables(34)/expressions(49)))+(observables(35)/expressions(50)))+(observables(36)/expressions(51)))+(observables(37)/expressions(52)))+(observables(26)/expressions(54)))+(observables(27)/expressions(55)))+(observables(28)/expressions(56)))+(observables(29)/expressions(57)))+(observables(30)/expressions(58)))+(observables(31)/expressions(59))));
end

% function rateLaw20
function [val] = rateLaw20(expressions, observables)
    val = ((expressions(61)/expressions(50))/((((((((((((1+(observables(32)/expressions(47)))+(observables(33)/expressions(48)))+(observables(34)/expressions(49)))+(observables(35)/expressions(50)))+(observables(36)/expressions(51)))+(observables(37)/expressions(52)))+(observables(26)/expressions(54)))+(observables(27)/expressions(55)))+(observables(28)/expressions(56)))+(observables(29)/expressions(57)))+(observables(30)/expressions(58)))+(observables(31)/expressions(59))));
end

% function rateLaw21
function [val] = rateLaw21(expressions, observables)
    val = ((expressions(61)/expressions(51))/((((((((((((1+(observables(32)/expressions(47)))+(observables(33)/expressions(48)))+(observables(34)/expressions(49)))+(observables(35)/expressions(50)))+(observables(36)/expressions(51)))+(observables(37)/expressions(52)))+(observables(26)/expressions(54)))+(observables(27)/expressions(55)))+(observables(28)/expressions(56)))+(observables(29)/expressions(57)))+(observables(30)/expressions(58)))+(observables(31)/expressions(59))));
end

% function rateLaw22
function [val] = rateLaw22(expressions, observables)
    val = ((expressions(61)/expressions(52))/((((((((((((1+(observables(32)/expressions(47)))+(observables(33)/expressions(48)))+(observables(34)/expressions(49)))+(observables(35)/expressions(50)))+(observables(36)/expressions(51)))+(observables(37)/expressions(52)))+(observables(26)/expressions(54)))+(observables(27)/expressions(55)))+(observables(28)/expressions(56)))+(observables(29)/expressions(57)))+(observables(30)/expressions(58)))+(observables(31)/expressions(59))));
end

% function rateLaw23
function [val] = rateLaw23(expressions, observables)
    val = (expressions(79)/(observables(50)+expressions(81)));
end

% function rateLaw24
function [val] = rateLaw24(expressions, observables)
    val = (expressions(78)/(observables(50)+expressions(80)));
end

% function rateLaw25
function [val] = rateLaw25(expressions, observables)
    val = (expressions(83)/(observables(51)+expressions(85)));
end

% function rateLaw26
function [val] = rateLaw26(expressions, observables)
    val = (expressions(82)/(observables(51)+expressions(84)));
end

% function rateLaw27
function [val] = rateLaw27(expressions, observables)
    val = ((observables(84)*expressions(87))/(observables(75)+expressions(88)));
end

% function rateLaw28
function [val] = rateLaw28(expressions, observables)
    val = ((observables(84)*expressions(87))/(observables(78)+expressions(88)));
end

% function rateLaw29
function [val] = rateLaw29(expressions, observables)
    val = ((observables(84)*expressions(89))/(observables(93)+expressions(90)));
end

% function rateLaw30
function [val] = rateLaw30(expressions, observables)
    val = (expressions(191)/(observables(56)+expressions(193)));
end

% function rateLaw31
function [val] = rateLaw31(expressions, observables)
    val = (expressions(192)/(observables(56)+expressions(194)));
end

% function rateLaw32
function [val] = rateLaw32(expressions, observables)
    val = (expressions(191)/(observables(57)+expressions(195)));
end

% function rateLaw33
function [val] = rateLaw33(expressions, observables)
    val = (expressions(192)/(observables(58)+expressions(196)));
end

% function rateLaw34
function [val] = rateLaw34(expressions, observables)
    val = (expressions(191)/(observables(59)+expressions(195)));
end

% function rateLaw35
function [val] = rateLaw35(expressions, observables)
    val = (expressions(191)/(observables(59)+expressions(193)));
end

% function rateLaw36
function [val] = rateLaw36(expressions, observables)
    val = (expressions(192)/(observables(59)+expressions(194)));
end

% function rateLaw37
function [val] = rateLaw37(expressions, observables)
    val = (expressions(191)/(observables(60)+expressions(193)));
end

% function rateLaw38
function [val] = rateLaw38(expressions, observables)
    val = (expressions(192)/(observables(60)+expressions(194)));
end

% function rateLaw39
function [val] = rateLaw39(expressions, observables)
    val = (expressions(191)/(observables(61)+expressions(195)));
end

% function rateLaw40
function [val] = rateLaw40(expressions, observables)
    val = (expressions(192)/(observables(61)+expressions(196)));
end

% function rateLaw41
function [val] = rateLaw41(expressions, observables)
    val = (expressions(192)/(observables(62)+expressions(196)));
end

% function rateLaw42
function [val] = rateLaw42(expressions, observables)
    val = (expressions(191)/(observables(71)+expressions(195)));
end

% function rateLaw43
function [val] = rateLaw43(expressions, observables)
    val = (expressions(191)/(observables(72)+expressions(193)));
end

% function rateLaw44
function [val] = rateLaw44(expressions, observables)
    val = (expressions(192)/(observables(72)+expressions(194)));
end

% function rateLaw45
function [val] = rateLaw45(expressions, observables)
    val = (expressions(197)/(observables(26)+expressions(198)));
end

% function rateLaw46
function [val] = rateLaw46(expressions, observables)
    val = (expressions(199)/(observables(27)+expressions(200)));
end

% function rateLaw47
function [val] = rateLaw47(expressions, observables)
    val = (expressions(201)/(observables(28)+expressions(202)));
end

% function rateLaw48
function [val] = rateLaw48(expressions, observables)
    val = (expressions(203)/(observables(29)+expressions(204)));
end

% function rateLaw49
function [val] = rateLaw49(expressions, observables)
    val = (expressions(205)/(observables(30)+expressions(206)));
end

% function rateLaw50
function [val] = rateLaw50(expressions, observables)
    val = (expressions(207)/(observables(31)+expressions(208)));
end

% function rateLaw51
function [val] = rateLaw51(expressions, observables)
    val = (expressions(209)/(observables(89)+expressions(210)));
end

% function rateLaw52
function [val] = rateLaw52(expressions, observables)
    val = (expressions(211)/(observables(90)+expressions(212)));
end

% function rateLaw53
function [val] = rateLaw53(expressions, observables)
    val = (expressions(209)/(observables(52)+expressions(210)));
end

% function rateLaw54
function [val] = rateLaw54(expressions, observables)
    val = (expressions(211)/(observables(53)+expressions(212)));
end

% function rateLaw55
function [val] = rateLaw55(expressions, observables)
    val = (expressions(213)/(observables(76)+expressions(214)));
end

% function rateLaw56
function [val] = rateLaw56(expressions, observables)
    val = (expressions(213)/(observables(79)+expressions(214)));
end

% function rateLaw57
function [val] = rateLaw57(expressions, observables)
    val = (expressions(215)/(observables(92)+expressions(216)));
end

% function rateLaw58
function [val] = rateLaw58(expressions, observables)
    val = (expressions(217)/(observables(46)+expressions(218)));
end

% function rateLaw59
function [val] = rateLaw59(expressions, observables)
    val = (expressions(219)/(observables(47)+expressions(220)));
end

% function rateLaw60
function [val] = rateLaw60(expressions, observables)
    val = (expressions(219)/(observables(47)+expressions(220)));
end

% function rateLaw61
function [val] = rateLaw61(expressions, observables)
    val = (expressions(221)/(observables(48)+expressions(222)));
end

% function rateLaw62
function [val] = rateLaw62(expressions, observables)
    val = (expressions(223)/(observables(49)+expressions(224)));
end

% function rateLaw63
function [val] = rateLaw63(expressions, observables)
    val = (expressions(140)/(observables(111)+expressions(141)));
end

% function rateLaw64
function [val] = rateLaw64(expressions, observables)
    val = (expressions(145)/(observables(115)+expressions(146)));
end

% function rateLaw65
function [val] = rateLaw65(expressions, observables)
    val = (((expressions(147)/expressions(148))/(1+(observables(126)/expressions(148))))*((1+(expressions(170)*((observables(133)/expressions(169))^2)))/(1+((observables(133)/expressions(169))^2))));
end

% function rateLaw66
function [val] = rateLaw66(expressions, observables)
    val = (expressions(149)/(expressions(150)+observables(127)));
end

% function rateLaw67
function [val] = rateLaw67(expressions, observables)
    val = ((expressions(151)/expressions(152))/((1+(observables(128)/expressions(152)))+(observables(129)/expressions(154))));
end

% function rateLaw68
function [val] = rateLaw68(expressions, observables)
    val = ((expressions(153)/expressions(154))/((1+(observables(128)/expressions(152)))+(observables(129)/expressions(154))));
end

% function rateLaw69
function [val] = rateLaw69(expressions, observables)
    val = ((expressions(155)/expressions(156))/(((1+(observables(130)/expressions(156)))+(observables(129)/expressions(158)))+(observables(128)/expressions(159))));
end

% function rateLaw70
function [val] = rateLaw70(expressions, observables)
    val = ((expressions(157)/expressions(158))/(((1+(observables(130)/expressions(156)))+(observables(129)/expressions(158)))+(observables(128)/expressions(159))));
end

% function rateLaw71
function [val] = rateLaw71(expressions, observables)
    val = ((expressions(160)/expressions(161))/((1+(observables(131)/expressions(161)))+(observables(132)/expressions(163))));
end

% function rateLaw72
function [val] = rateLaw72(expressions, observables)
    val = ((expressions(162)/expressions(163))/((1+(observables(131)/expressions(161)))+(observables(132)/expressions(163))));
end

% function rateLaw73
function [val] = rateLaw73(expressions, observables)
    val = ((expressions(164)/expressions(165))/(((1+(observables(133)/expressions(165)))+(observables(132)/expressions(167)))+(observables(131)/expressions(168))));
end

% function rateLaw74
function [val] = rateLaw74(expressions, observables)
    val = ((expressions(166)/expressions(167))/(((1+(observables(133)/expressions(165)))+(observables(132)/expressions(167)))+(observables(131)/expressions(168))));
end

% function rateLaw75
function [val] = rateLaw75(expressions, observables)
    val = (expressions(174)/(observables(98)+expressions(175)));
end

% function rateLaw76
function [val] = rateLaw76(expressions, observables)
    val = (expressions(174)/(observables(99)+expressions(175)));
end

% function rateLaw77
function [val] = rateLaw77(expressions, observables)
    val = (expressions(174)/(observables(100)+expressions(175)));
end

% function rateLaw78
function [val] = rateLaw78(expressions, observables)
    val = (expressions(174)/(observables(101)+expressions(175)));
end

% function rateLaw79
function [val] = rateLaw79(expressions, observables)
    val = (expressions(174)/(observables(102)+expressions(175)));
end

% function rateLaw80
function [val] = rateLaw80(expressions, observables)
    val = (expressions(174)/(observables(103)+expressions(175)));
end

% function rateLaw81
function [val] = rateLaw81(expressions, observables)
    val = (expressions(176)/(observables(98)+expressions(177)));
end

% function rateLaw82
function [val] = rateLaw82(expressions, observables)
    val = (expressions(176)/(observables(99)+expressions(177)));
end

% function rateLaw83
function [val] = rateLaw83(expressions, observables)
    val = (expressions(176)/(observables(100)+expressions(177)));
end

% function rateLaw84
function [val] = rateLaw84(expressions, observables)
    val = (expressions(176)/(observables(101)+expressions(177)));
end

% function rateLaw85
function [val] = rateLaw85(expressions, observables)
    val = (expressions(176)/(observables(102)+expressions(177)));
end

% function rateLaw86
function [val] = rateLaw86(expressions, observables)
    val = (expressions(176)/(observables(103)+expressions(177)));
end

% function rateLaw87
function [val] = rateLaw87(expressions, observables)
    val = (expressions(181)/(observables(109)+expressions(182)));
end

% function rateLaw88
function [val] = rateLaw88(expressions, observables)
    val = (expressions(181)/(observables(110)+expressions(182)));
end

% function rateLaw89
function [val] = rateLaw89(expressions, observables)
    val = (expressions(179)/(observables(105)+expressions(180)));
end

% function rateLaw90
function [val] = rateLaw90(expressions, observables)
    val = (expressions(179)/(observables(105)+expressions(180)));
end

% function rateLaw91
function [val] = rateLaw91(expressions, observables)
    val = (expressions(181)/(observables(26)+expressions(182)));
end

% function rateLaw92
function [val] = rateLaw92(expressions, observables)
    val = (expressions(181)/(observables(27)+expressions(182)));
end

% function rateLaw93
function [val] = rateLaw93(expressions, observables)
    val = (expressions(181)/(observables(28)+expressions(182)));
end

% function rateLaw94
function [val] = rateLaw94(expressions, observables)
    val = (expressions(181)/(observables(29)+expressions(182)));
end

% function rateLaw95
function [val] = rateLaw95(expressions, observables)
    val = (expressions(181)/(observables(30)+expressions(182)));
end

% function rateLaw96
function [val] = rateLaw96(expressions, observables)
    val = (expressions(181)/(observables(31)+expressions(182)));
end

% function rateLaw97
function [val] = rateLaw97(expressions, observables)
    val = (expressions(183)/(observables(107)+expressions(184)));
end

% function rateLaw98
function [val] = rateLaw98(expressions, observables)
    val = (expressions(183)/(observables(108)+expressions(184)));
end

% function rateLaw99
function [val] = rateLaw99(expressions, observables)
    val = (expressions(185)/(observables(104)+expressions(186)));
end

% function rateLaw100
function [val] = rateLaw100(expressions, observables)
    val = (expressions(187)/(observables(106)+expressions(188)));
end




% Calculate expressions
function [ expressions ] = calc_expressions ( parameters )

    expressions = zeros(1,226);
    expressions(1) = parameters(1);
    expressions(2) = (6*expressions(1));
    expressions(3) = parameters(2);
    expressions(4) = parameters(3);
    expressions(5) = parameters(4);
    expressions(6) = parameters(5);
    expressions(7) = parameters(6);
    expressions(8) = parameters(7);
    expressions(9) = parameters(8);
    expressions(10) = parameters(9);
    expressions(11) = parameters(10);
    expressions(12) = parameters(11);
    expressions(13) = parameters(12);
    expressions(14) = parameters(13);
    expressions(15) = parameters(14);
    expressions(16) = parameters(15);
    expressions(17) = parameters(16);
    expressions(18) = parameters(17);
    expressions(19) = parameters(18);
    expressions(20) = parameters(19);
    expressions(21) = parameters(20);
    expressions(22) = parameters(21);
    expressions(23) = parameters(22);
    expressions(24) = parameters(23);
    expressions(25) = parameters(24);
    expressions(26) = parameters(25);
    expressions(27) = (expressions(25)*expressions(26));
    expressions(28) = parameters(26);
    expressions(29) = parameters(27);
    expressions(30) = parameters(28);
    expressions(31) = parameters(29);
    expressions(32) = parameters(30);
    expressions(33) = parameters(31);
    expressions(34) = parameters(32);
    expressions(35) = parameters(33);
    expressions(36) = parameters(34);
    expressions(37) = parameters(35);
    expressions(38) = parameters(36);
    expressions(39) = expressions(29);
    expressions(40) = expressions(29);
    expressions(41) = expressions(29);
    expressions(42) = parameters(37);
    expressions(43) = parameters(38);
    expressions(44) = parameters(39);
    expressions(45) = parameters(40);
    expressions(46) = parameters(41);
    expressions(47) = parameters(42);
    expressions(48) = parameters(43);
    expressions(49) = parameters(44);
    expressions(50) = parameters(45);
    expressions(51) = parameters(46);
    expressions(52) = parameters(47);
    expressions(53) = parameters(48);
    expressions(54) = (expressions(47)*expressions(53));
    expressions(55) = (expressions(48)*expressions(53));
    expressions(56) = (expressions(49)*expressions(53));
    expressions(57) = (expressions(50)*expressions(53));
    expressions(58) = (expressions(51)*expressions(53));
    expressions(59) = (expressions(52)*expressions(53));
    expressions(60) = parameters(49);
    expressions(61) = (expressions(60)/expressions(30));
    expressions(62) = parameters(50);
    expressions(63) = parameters(51);
    expressions(64) = parameters(52);
    expressions(65) = parameters(53);
    expressions(66) = parameters(54);
    expressions(67) = parameters(55);
    expressions(68) = (expressions(64)/expressions(30));
    expressions(69) = (expressions(65)/expressions(30));
    expressions(70) = (expressions(66)/expressions(30));
    expressions(71) = (expressions(67)/expressions(30));
    expressions(72) = parameters(56);
    expressions(73) = parameters(57);
    expressions(74) = (expressions(72)*expressions(73));
    expressions(75) = parameters(58);
    expressions(76) = parameters(59);
    expressions(77) = (expressions(75)*expressions(76));
    expressions(78) = parameters(60);
    expressions(79) = (expressions(78)/expressions(30));
    expressions(80) = expressions(29);
    expressions(81) = expressions(80);
    expressions(82) = expressions(78);
    expressions(83) = (expressions(78)/expressions(30));
    expressions(84) = expressions(29);
    expressions(85) = expressions(80);
    expressions(86) = parameters(61);
    expressions(87) = expressions(86);
    expressions(88) = expressions(29);
    expressions(89) = expressions(86);
    expressions(90) = expressions(29);
    expressions(91) = parameters(62);
    expressions(92) = parameters(63);
    expressions(93) = (expressions(92)*expressions(91));
    expressions(94) = parameters(64);
    expressions(95) = parameters(65);
    expressions(96) = (expressions(95)*expressions(94));
    expressions(97) = parameters(66);
    expressions(98) = parameters(67);
    expressions(99) = (expressions(98)*expressions(97));
    expressions(100) = parameters(68);
    expressions(101) = parameters(69);
    expressions(102) = (expressions(101)*expressions(100));
    expressions(103) = parameters(70);
    expressions(104) = parameters(71);
    expressions(105) = (expressions(104)*expressions(103));
    expressions(106) = parameters(72);
    expressions(107) = parameters(73);
    expressions(108) = (expressions(107)*expressions(106));
    expressions(109) = parameters(74);
    expressions(110) = parameters(75);
    expressions(111) = (expressions(110)*expressions(109));
    expressions(112) = parameters(76);
    expressions(113) = parameters(77);
    expressions(114) = (expressions(113)*expressions(112));
    expressions(115) = parameters(78);
    expressions(116) = parameters(79);
    expressions(117) = (expressions(116)*expressions(115));
    expressions(118) = parameters(80);
    expressions(119) = parameters(81);
    expressions(120) = (expressions(119)*expressions(118));
    expressions(121) = parameters(82);
    expressions(122) = parameters(83);
    expressions(123) = (expressions(122)*expressions(121));
    expressions(124) = parameters(84);
    expressions(125) = parameters(85);
    expressions(126) = parameters(86);
    expressions(127) = (expressions(126)*expressions(125));
    expressions(128) = parameters(87);
    expressions(129) = parameters(88);
    expressions(130) = parameters(89);
    expressions(131) = (expressions(130)*expressions(129));
    expressions(132) = parameters(90);
    expressions(133) = parameters(91);
    expressions(134) = parameters(92);
    expressions(135) = (expressions(134)*expressions(133));
    expressions(136) = parameters(93);
    expressions(137) = parameters(94);
    expressions(138) = parameters(95);
    expressions(139) = (expressions(138)*expressions(137));
    expressions(140) = parameters(96);
    expressions(141) = expressions(29);
    expressions(142) = parameters(97);
    expressions(143) = parameters(98);
    expressions(144) = (expressions(143)*expressions(142));
    expressions(145) = parameters(99);
    expressions(146) = expressions(29);
    expressions(147) = parameters(100);
    expressions(148) = parameters(101);
    expressions(149) = parameters(102);
    expressions(150) = parameters(103);
    expressions(151) = parameters(104);
    expressions(152) = parameters(105);
    expressions(153) = parameters(106);
    expressions(154) = parameters(107);
    expressions(155) = parameters(108);
    expressions(156) = parameters(109);
    expressions(157) = parameters(110);
    expressions(158) = parameters(111);
    expressions(159) = parameters(112);
    expressions(160) = parameters(113);
    expressions(161) = parameters(114);
    expressions(162) = parameters(115);
    expressions(163) = parameters(116);
    expressions(164) = parameters(117);
    expressions(165) = parameters(118);
    expressions(166) = parameters(119);
    expressions(167) = parameters(120);
    expressions(168) = parameters(121);
    expressions(169) = parameters(122);
    expressions(170) = parameters(123);
    expressions(171) = parameters(124);
    expressions(172) = parameters(125);
    expressions(173) = (expressions(172)*expressions(171));
    expressions(174) = parameters(126);
    expressions(175) = expressions(29);
    expressions(176) = (expressions(174)/expressions(30));
    expressions(177) = expressions(175);
    expressions(178) = parameters(127);
    expressions(179) = expressions(178);
    expressions(180) = expressions(29);
    expressions(181) = expressions(178);
    expressions(182) = expressions(29);
    expressions(183) = expressions(178);
    expressions(184) = expressions(29);
    expressions(185) = parameters(128);
    expressions(186) = expressions(29);
    expressions(187) = expressions(178);
    expressions(188) = expressions(29);
    expressions(189) = parameters(129);
    expressions(190) = expressions(29);
    expressions(191) = parameters(130);
    expressions(192) = parameters(131);
    expressions(193) = expressions(190);
    expressions(194) = expressions(190);
    expressions(195) = expressions(190);
    expressions(196) = expressions(190);
    expressions(197) = parameters(132);
    expressions(198) = expressions(29);
    expressions(199) = expressions(197);
    expressions(200) = expressions(198);
    expressions(201) = expressions(197);
    expressions(202) = expressions(198);
    expressions(203) = expressions(197);
    expressions(204) = expressions(198);
    expressions(205) = expressions(197);
    expressions(206) = expressions(198);
    expressions(207) = expressions(197);
    expressions(208) = expressions(198);
    expressions(209) = parameters(133);
    expressions(210) = expressions(29);
    expressions(211) = parameters(134);
    expressions(212) = expressions(29);
    expressions(213) = expressions(189);
    expressions(214) = expressions(190);
    expressions(215) = expressions(189);
    expressions(216) = expressions(190);
    expressions(217) = parameters(135);
    expressions(218) = expressions(29);
    expressions(219) = expressions(217);
    expressions(220) = expressions(218);
    expressions(221) = expressions(217);
    expressions(222) = expressions(220);
    expressions(223) = expressions(217);
    expressions(224) = expressions(220);
    expressions(225) = (expressions(4)*0);
    expressions(226) = (expressions(4)*0);
   
end



% Calculate observables
function [ observables ] = calc_observables ( species, expressions )

    observables = zeros(1,137);
    observables(1) = species(49) +species(50) +species(51) +species(52) +species(53) +species(54) +species(66) +species(67) +species(68) +species(69) +species(70) +species(71) +2*species(78) +2*species(79) +2*species(80) +species(93) +species(94) +species(95) +species(96) +species(97) +species(98) +2*species(99) +species(100) +species(101) +2*species(102) +species(103) +species(104) +2*species(105) +species(106) +species(107) +2*species(113) +2*species(114) +2*species(115) +2*species(129) +2*species(130) +2*species(131) +2*species(136) +2*species(137) +2*species(138) +2*species(147) +2*species(148) +2*species(149) +2*species(152) +2*species(153) +2*species(154) +2*species(166) +2*species(167) +2*species(168) +2*species(181) +2*species(182) +2*species(183) +2*species(185) +2*species(186) +2*species(187);
    observables(2) = species(81) +species(82) +species(83) +species(108) +species(109) +species(110) +species(111) +species(112) +species(121) +species(132) +species(133) +species(134) +2*species(135) +2*species(150) +2*species(156) +2*species(157) +2*species(158) +2*species(159) +2*species(160) +2*species(161) +2*species(162) +2*species(163) +species(165) +2*species(169) +2*species(170) +2*species(171) +2*species(172) +2*species(173) +2*species(174) +2*species(175) +2*species(176) +species(177) +2*species(179) +2*species(180) +2*species(188) +2*species(189) +species(212) +species(216) +species(221) +species(223);
    observables(3) = species(42) +species(61) +species(62) +species(64) +species(65) +species(84) +species(85) +species(86) +species(87) +species(88) +species(89) +species(90) +species(116) +species(117) +species(118) +species(119) +species(120) +species(125) +species(126) +species(139) +species(140) +species(141) +species(143) +species(155);
    observables(4) = species(120) +species(139) +species(142) +species(143) +species(155) +species(180) +species(189);
    observables(5) = species(5) +species(7) +2*species(37) +species(38) +species(39) +species(41) +species(47) +species(55) +species(56) +species(57) +species(58) +species(59) +species(60) +species(72) +species(73) +species(74) +species(75) +species(76) +species(77) +species(111) +species(112) +species(133) +species(134) +species(191) +species(193) +2*species(195) +2*species(196) +2*species(197) +species(198) +species(199) +species(200) +species(201) +species(202) +species(203) +species(207) +species(209) +species(210) +species(211) +species(212) +species(213) +species(214) +species(215) +species(216) +species(217) +species(218) +species(219) +species(220) +species(221) +species(222) +species(223) +species(224) +species(225);
    observables(6) = species(7) +species(37) +species(38) +species(41) +species(44) +species(47) +species(48) +species(58) +species(59) +species(60) +species(73) +species(75) +species(77) +species(112) +species(134) +species(193) +species(194) +species(195) +species(196) +species(197) +species(198) +species(199) +species(200) +species(201) +species(202) +species(203) +species(204) +species(205) +species(206) +species(209) +species(214) +species(215) +species(216) +species(217) +species(219) +species(222) +species(223) +species(225);
    observables(7) = species(136) +species(137) +species(138) +species(147) +species(148) +species(149) +species(151) +species(152) +species(153) +species(154) +species(166) +species(167) +species(168) +species(178);
    observables(8) = species(152) +species(153) +species(154) +species(166) +species(167) +species(168) +species(178) +species(181) +species(182) +species(183) +species(185) +species(186) +species(187) +species(190);
    observables(9) = species(121) +species(132);
    observables(10) = species(62) +species(84) +species(87) +species(116) +species(120) +species(139) +species(143) +species(155) +species(156) +species(169) +species(179) +species(180) +species(188) +species(189);
    observables(11) = species(92);
    observables(12) = species(146);
    observables(13) = species(184);
    observables(14) = species(6) +species(38) +species(40) +species(48) +species(192) +species(198) +species(199) +species(200) +species(204) +species(205) +species(206) +species(208);
    observables(15) = species(7) +species(37) +species(38) +species(41) +species(58) +species(59) +species(60) +species(73) +species(75) +species(77) +species(112) +species(134) +species(193) +species(195) +species(196) +species(197) +species(198) +species(199) +species(200) +species(209) +species(214) +species(215) +species(216) +species(217) +species(219) +species(222) +species(223) +species(225);
    observables(16) = species(5) +species(37) +species(39) +species(47) +species(55) +species(56) +species(57) +species(72) +species(74) +species(76) +species(111) +species(133) +species(191) +species(195) +species(196) +species(197) +species(201) +species(202) +species(203) +species(207) +species(210) +species(211) +species(212) +species(213) +species(218) +species(220) +species(221) +species(224);
    observables(17) = species(44) +species(47) +species(48) +species(194) +species(201) +species(202) +species(203) +species(204) +species(205) +species(206);
    observables(18) = species(6) +species(192);
    observables(19) = species(49) +species(66) +species(78) +species(93) +species(99) +species(100) +species(113) +species(129) +species(136) +species(147) +species(152) +species(166) +species(181) +species(185);
    observables(20) = species(50) +species(67) +species(78) +species(94) +species(99) +species(101) +species(113) +species(129) +species(136) +species(147) +species(152) +species(166) +species(181) +species(185);
    observables(21) = species(51) +species(68) +species(79) +species(95) +species(102) +species(103) +species(114) +species(130) +species(137) +species(148) +species(153) +species(167) +species(182) +species(186);
    observables(22) = species(52) +species(69) +species(79) +species(96) +species(102) +species(104) +species(114) +species(130) +species(137) +species(148) +species(153) +species(167) +species(182) +species(186);
    observables(23) = species(53) +species(70) +species(80) +species(97) +species(105) +species(106) +species(115) +species(131) +species(138) +species(149) +species(154) +species(168) +species(183) +species(187);
    observables(24) = species(54) +species(71) +species(80) +species(98) +species(105) +species(107) +species(115) +species(131) +species(138) +species(149) +species(154) +species(168) +species(183) +species(187);
    observables(25) = species(128);
    observables(26) = species(49) +species(66) +species(78) +species(99);
    observables(27) = species(50) +species(67) +species(78) +species(99);
    observables(28) = species(51) +species(68) +species(79) +species(102);
    observables(29) = species(52) +species(69) +species(79) +species(102);
    observables(30) = species(53) +species(70) +species(80) +species(105);
    observables(31) = species(54) +species(71) +species(80) +species(105);
    observables(32) = species(2) +species(30) +species(50) +species(67);
    observables(33) = species(2) +species(30) +species(49) +species(66);
    observables(34) = species(3) +species(31) +species(52) +species(69);
    observables(35) = species(3) +species(31) +species(51) +species(68);
    observables(36) = species(4) +species(32) +species(54) +species(71);
    observables(37) = species(4) +species(32) +species(53) +species(70);
    observables(38) = species(6) +species(38) +species(40) +species(44) +species(47) +2*species(48) +species(192) +species(194) +species(198) +species(199) +species(200) +species(201) +species(202) +species(203) +2*species(204) +2*species(205) +2*species(206) +species(208);
    observables(39) = species(5) +species(6) +species(37) +species(38) +species(39) +species(40) +species(47) +species(48) +species(55) +species(56) +species(57) +species(72) +species(74) +species(76) +species(111) +species(133) +species(191) +species(192) +species(195) +species(196) +species(197) +species(198) +species(199) +species(200) +species(201) +species(202) +species(203) +species(204) +species(205) +species(206) +species(207) +species(208) +species(210) +species(211) +species(212) +species(213) +species(218) +species(220) +species(221) +species(224);
    observables(40) = species(5) +species(7) +2*species(37) +species(38) +species(39) +species(41) +species(47) +species(55) +species(56) +species(57) +species(58) +species(59) +species(60) +species(72) +species(73) +species(74) +species(75) +species(76) +species(77) +species(111) +species(112) +species(133) +species(134) +species(191) +species(193) +2*species(195) +2*species(196) +2*species(197) +species(198) +species(199) +species(200) +species(201) +species(202) +species(203) +species(207) +species(209) +species(210) +species(211) +species(212) +species(213) +species(214) +species(215) +species(216) +species(217) +species(218) +species(219) +species(220) +species(221) +species(222) +species(223) +species(224) +species(225);
    observables(41) = species(7) +species(37) +species(38) +species(41) +species(44) +species(47) +species(48) +species(58) +species(59) +species(60) +species(73) +species(75) +species(77) +species(112) +species(134) +species(193) +species(194) +species(195) +species(196) +species(197) +species(198) +species(199) +species(200) +species(201) +species(202) +species(203) +species(204) +species(205) +species(206) +species(209) +species(214) +species(215) +species(216) +species(217) +species(219) +species(222) +species(223) +species(225);
    observables(42) = species(9) +species(33);
    observables(43) = species(10) +species(34) +species(56) +species(59) +species(74) +species(75) +species(82) +species(109) +species(211) +species(215) +species(220) +species(222);
    observables(44) = species(10) +species(34) +species(165) +species(177);
    observables(45) = species(11) +species(35);
    observables(46) = species(81) +species(108);
    observables(47) = species(135) +species(150) +species(165) +species(177);
    observables(48) = species(82) +species(109) +species(111) +species(112) +species(133) +species(134) +species(135) +species(150) +species(212) +species(216) +species(221) +species(223);
    observables(49) = species(83) +species(110);
    observables(50) = species(113) +species(114) +species(115) +species(129) +species(130) +species(131) +species(181) +species(182) +species(183) +species(185) +species(186) +species(187);
    observables(51) = species(113) +species(114) +species(115) +species(129) +species(130) +species(131) +species(136) +species(137) +species(138) +species(147) +species(148) +species(149);
    observables(52) = species(151) +species(178);
    observables(53) = species(178) +species(190);
    observables(54) = species(17);
    observables(55) = species(6) +species(192);
    observables(56) = species(7) +species(193);
    observables(57) = species(5) +species(191);
    observables(58) = species(44) +species(194);
    observables(59) = species(37) +species(195) +species(196) +species(197);
    observables(60) = species(38) +species(198) +species(199) +species(200);
    observables(61) = species(47) +species(201) +species(202) +species(203);
    observables(62) = species(48) +species(204) +species(205) +species(206);
    observables(63) = species(6) +species(192);
    observables(64) = species(7) +species(193);
    observables(65) = species(5) +species(191);
    observables(66) = species(44) +species(194);
    observables(67) = species(38) +species(40) +species(48) +species(198) +species(199) +species(200) +species(204) +species(205) +species(206) +species(208);
    observables(68) = species(37) +species(39) +species(47) +species(55) +species(56) +species(57) +species(72) +species(74) +species(76) +species(111) +species(133) +species(195) +species(196) +species(197) +species(201) +species(202) +species(203) +species(207) +species(210) +species(211) +species(212) +species(213) +species(218) +species(220) +species(221) +species(224);
    observables(69) = species(47) +species(48) +species(201) +species(202) +species(203) +species(204) +species(205) +species(206);
    observables(70) = species(37) +species(38) +species(41) +species(58) +species(59) +species(60) +species(73) +species(75) +species(77) +species(112) +species(134) +species(195) +species(196) +species(197) +species(198) +species(199) +species(200) +species(209) +species(214) +species(215) +species(216) +species(217) +species(219) +species(222) +species(223) +species(225);
    observables(71) = species(39) +species(207);
    observables(72) = species(41) +species(209);
    observables(73) = species(8) +species(39) +species(40) +species(41) +species(207) +species(208) +species(209);
    observables(74) = species(8);
    observables(75) = species(12) +species(43);
    observables(76) = species(42) +species(61);
    observables(77) = species(42) +species(61) +species(62) +species(64) +species(65) +species(84) +species(85) +species(86) +species(87) +species(88) +species(89) +species(90) +species(116) +species(117) +species(118) +species(119) +species(120) +species(125) +species(126) +species(139) +species(140) +species(141) +species(143) +species(155);
    observables(78) = species(12) +species(42) +species(62) +species(64) +species(65) +species(87) +species(88) +species(89) +species(90) +species(120) +species(125) +species(126) +species(143);
    observables(79) = species(43) +species(61) +species(84) +species(85) +species(86) +species(116) +species(117) +species(118) +species(119) +species(139) +species(140) +species(141) +species(155);
    observables(80) = species(43) +species(61) +species(84) +species(85) +species(86) +species(116) +species(117) +species(118) +species(119) +species(139) +species(140) +species(141) +species(155);
    observables(81) = species(13);
    observables(82) = species(151);
    observables(83) = species(190);
    observables(84) = species(178);
    observables(85) = species(113) +species(114) +species(115) +species(129) +species(130) +species(131);
    observables(86) = species(136) +species(137) +species(138) +species(147) +species(148) +species(149);
    observables(87) = species(181) +species(182) +species(183) +species(185) +species(186) +species(187);
    observables(88) = species(152) +species(153) +species(154) +species(166) +species(167) +species(168);
    observables(89) = species(136) +species(137) +species(138) +species(147) +species(148) +species(149) +species(152) +species(153) +species(154) +species(166) +species(167) +species(168);
    observables(90) = species(152) +species(153) +species(154) +species(166) +species(167) +species(168) +species(181) +species(182) +species(183) +species(185) +species(186) +species(187);
    observables(91) = species(14);
    observables(92) = species(142);
    observables(93) = species(87) +species(116) +species(179) +species(188);
    observables(94) = species(120) +species(139) +species(180) +species(189);
    observables(95) = species(120) +species(139) +species(143) +species(155) +species(180) +species(189);
    observables(96) = species(16);
    observables(97) = species(121) +species(132);
    observables(98) = species(93) +species(100);
    observables(99) = species(94) +species(101);
    observables(100) = species(95) +species(103);
    observables(101) = species(96) +species(104);
    observables(102) = species(97) +species(106);
    observables(103) = species(98) +species(107);
    observables(104) = species(128);
    observables(105) = species(5) +species(7);
    observables(106) = species(5) +species(6) +species(7) +species(44);
    observables(107) = species(136) +species(137) +species(138) +species(147) +species(148) +species(149) +species(151) +species(152) +species(153) +species(154) +species(166) +species(167) +species(168) +species(178);
    observables(108) = species(152) +species(153) +species(154) +species(166) +species(167) +species(168) +species(178) +species(181) +species(182) +species(183) +species(185) +species(186) +species(187) +species(190);
    observables(109) = species(135) +species(150) +species(165) +species(177);
    observables(110) = species(82) +species(109) +species(135) +species(150);
    observables(111) = species(27);
    observables(112) = species(143) +species(155);
    observables(113) = species(28) +species(45);
    observables(114) = species(45);
    observables(115) = species(21);
    observables(116) = species(63);
    observables(117) = species(88) +species(90) +species(117) +species(119) +species(122) +species(124) +2*species(125) +species(126) +2*species(140) +species(141) +2*species(144) +species(145) +species(159) +species(161) +2*species(162) +species(163) +species(172) +species(174) +2*species(175) +species(176);
    observables(118) = species(89) +species(91) +species(118) +species(123) +species(126) +species(141) +species(145) +species(160) +species(163) +species(173) +species(176);
    observables(119) = species(64) +species(65) +species(85) +species(86) +species(88) +species(89) +species(90) +species(117) +species(118) +species(119) +species(125) +species(126) +species(140) +species(141) +species(157) +species(158) +species(159) +species(160) +species(161) +species(162) +species(163) +species(170) +species(171) +species(172) +species(173) +species(174) +species(175) +species(176);
    observables(120) = species(19) +species(46) +species(122) +species(123) +species(124) +species(144) +species(145);
    observables(121) = species(65) +species(86) +species(90) +species(119) +species(158) +species(161) +species(171) +species(174);
    observables(122) = species(88) +species(117) +species(125) +species(140) +species(159) +species(162) +species(172) +species(175);
    observables(123) = species(89) +species(118) +species(126) +species(141) +species(160) +species(163) +species(173) +species(176);
    observables(124) = species(21);
    observables(125) = species(63);
    observables(126) = species(22);
    observables(127) = species(92);
    observables(128) = species(23);
    observables(129) = species(127);
    observables(130) = species(146);
    observables(131) = species(24);
    observables(132) = species(164);
    observables(133) = species(184);
    observables(134) = species(30) +species(49) +species(50) +species(78) +species(93) +species(94) +species(113) +species(136) +species(152) +species(181);
    observables(135) = species(2) +species(66) +species(67) +species(99) +species(100) +species(101) +species(129) +species(147) +species(166) +species(185);
    observables(136) = species(49) +species(50) +species(51) +species(52) +species(53) +species(54) +species(66) +species(67) +species(68) +species(69) +species(70) +species(71) +species(93) +species(94) +species(95) +species(96) +species(97) +species(98) +species(100) +species(101) +species(103) +species(104) +species(106) +species(107);
    observables(137) = species(78) +species(79) +species(80) +species(99) +species(102) +species(105) +species(113) +species(114) +species(115) +species(129) +species(130) +species(131) +species(136) +species(137) +species(138) +species(147) +species(148) +species(149) +species(152) +species(153) +species(154) +species(166) +species(167) +species(168) +species(181) +species(182) +species(183) +species(185) +species(186) +species(187);

end


% Calculate ratelaws
function [ ratelaws ] = calcratelaws ( species, expressions, observables )

    ratelaws = zeros(1,137);
    ratelaws(1) = expressions(26)*species(2)*species(1);
    ratelaws(2) = expressions(26)*species(3)*species(1);
    ratelaws(3) = expressions(26)*species(4)*species(1);
    ratelaws(4) = expressions(26)*species(9)*species(1);
    ratelaws(5) = expressions(26)*species(10)*species(1);
    ratelaws(6) = expressions(26)*species(11)*species(1);
    ratelaws(7) = rateLaw1(expressions,observables)*species(17);
    ratelaws(8) = rateLaw2(expressions,observables)*species(6);
    ratelaws(9) = rateLaw3(expressions,observables)*species(5);
    ratelaws(10) = rateLaw5(expressions,observables)*species(6);
    ratelaws(11) = rateLaw6(expressions,observables)*species(5);
    ratelaws(12) = rateLaw8(expressions,observables)*species(6);
    ratelaws(13) = rateLaw9(expressions,observables)*species(5);
    ratelaws(14) = expressions(42)*species(5)*species(7);
    ratelaws(15) = expressions(42)*species(6)*species(7);
    ratelaws(16) = expressions(34)*species(5)*species(8);
    ratelaws(17) = expressions(34)*species(6)*species(8);
    ratelaws(18) = expressions(34)*species(7)*species(8);
    ratelaws(19) = rateLaw27(expressions,observables)*species(12);
    ratelaws(20) = rateLaw28(expressions,observables)*species(12);
    ratelaws(21) = rateLaw30(expressions,observables)*species(7)*species(17);
    ratelaws(22) = rateLaw31(expressions,observables)*species(7)*species(17);
    ratelaws(23) = rateLaw32(expressions,observables)*species(5)*species(17);
    ratelaws(24) = expressions(142)*species(28)*species(26);
    ratelaws(25) = expressions(109)*species(18)*species(19);
    ratelaws(26) = (expressions(25)*expressions(26))*species(30);
    ratelaws(27) = (expressions(25)*expressions(26))*species(31);
    ratelaws(28) = (expressions(25)*expressions(26))*species(32);
    ratelaws(29) = (expressions(25)*expressions(26))*species(33);
    ratelaws(30) = (expressions(25)*expressions(26))*species(34);
    ratelaws(31) = (expressions(25)*expressions(26))*species(35);
    ratelaws(32) = expressions(28)*species(36);
    ratelaws(33) = rateLaw4(expressions,observables)*species(44);
    ratelaws(34) = rateLaw7(expressions,observables)*species(44);
    ratelaws(35) = rateLaw10(expressions,observables)*species(44);
    ratelaws(36) = expressions(43)*species(37);
    ratelaws(37) = expressions(44)*species(38);
    ratelaws(38) = expressions(42)*species(5)*species(44);
    ratelaws(39) = expressions(42)*species(6)*species(44);
    ratelaws(40) = expressions(35)*species(39);
    ratelaws(41) = expressions(37)*species(39);
    ratelaws(42) = expressions(36)*species(40);
    ratelaws(43) = expressions(38)*species(40);
    ratelaws(44) = expressions(35)*species(41);
    ratelaws(45) = rateLaw11(expressions,observables)*species(30)*species(5);
    ratelaws(46) = rateLaw12(expressions,observables)*species(30)*species(5);
    ratelaws(47) = rateLaw13(expressions,observables)*species(31)*species(5);
    ratelaws(48) = rateLaw14(expressions,observables)*species(31)*species(5);
    ratelaws(49) = rateLaw15(expressions,observables)*species(32)*species(5);
    ratelaws(50) = rateLaw16(expressions,observables)*species(32)*species(5);
    ratelaws(51) = rateLaw17(expressions,observables)*species(30)*species(7);
    ratelaws(52) = rateLaw18(expressions,observables)*species(30)*species(7);
    ratelaws(53) = rateLaw19(expressions,observables)*species(31)*species(7);
    ratelaws(54) = rateLaw20(expressions,observables)*species(31)*species(7);
    ratelaws(55) = rateLaw21(expressions,observables)*species(32)*species(7);
    ratelaws(56) = rateLaw22(expressions,observables)*species(32)*species(7);
    ratelaws(57) = expressions(62)*species(33)*species(5);
    ratelaws(58) = 2*expressions(62)*species(34)*species(5);
    ratelaws(59) = expressions(62)*species(35)*species(5);
    ratelaws(60) = expressions(62)*species(33)*species(7);
    ratelaws(61) = 2*expressions(62)*species(34)*species(7);
    ratelaws(62) = expressions(62)*species(35)*species(7);
    ratelaws(63) = rateLaw27(expressions,observables)*species(43);
    ratelaws(64) = rateLaw28(expressions,observables)*species(42);
    ratelaws(65) = expressions(100)*species(42)*species(15);
    ratelaws(66) = rateLaw33(expressions,observables)*species(44)*species(17);
    ratelaws(67) = rateLaw34(expressions,observables)*species(37)*species(17);
    ratelaws(68) = rateLaw35(expressions,observables)*species(37)*species(17);
    ratelaws(69) = rateLaw36(expressions,observables)*species(37)*species(17);
    ratelaws(70) = rateLaw37(expressions,observables)*species(38)*species(17);
    ratelaws(71) = rateLaw38(expressions,observables)*species(38)*species(17);
    ratelaws(72) = rateLaw42(expressions,observables)*species(39)*species(17);
    ratelaws(73) = rateLaw43(expressions,observables)*species(41)*species(17);
    ratelaws(74) = rateLaw44(expressions,observables)*species(41)*species(17);
    ratelaws(75) = rateLaw55(expressions,observables)*species(42)*species(17);
    ratelaws(76) = rateLaw56(expressions,observables)*species(43)*species(17);
    ratelaws(77) = (expressions(143)*expressions(142))*species(45);
    ratelaws(78) = rateLaw64(expressions,observables)*species(45)*species(21);
    ratelaws(79) = expressions(106)*species(42)*species(18);
    ratelaws(80) = expressions(106)*species(42)*species(46);
    ratelaws(81) = (expressions(110)*expressions(109))*species(46);
    ratelaws(82) = (expressions(25)*expressions(26))*species(49);
    ratelaws(83) = (expressions(25)*expressions(26))*species(50);
    ratelaws(84) = (expressions(25)*expressions(26))*species(51);
    ratelaws(85) = (expressions(25)*expressions(26))*species(52);
    ratelaws(86) = (expressions(25)*expressions(26))*species(53);
    ratelaws(87) = (expressions(25)*expressions(26))*species(54);
    ratelaws(88) = (expressions(25)*expressions(26))*species(55);
    ratelaws(89) = (expressions(25)*expressions(26))*species(58);
    ratelaws(90) = (expressions(25)*expressions(26))*species(56);
    ratelaws(91) = (expressions(25)*expressions(26))*species(59);
    ratelaws(92) = (expressions(25)*expressions(26))*species(57);
    ratelaws(93) = (expressions(25)*expressions(26))*species(60);
    ratelaws(94) = expressions(45)*species(47);
    ratelaws(95) = expressions(46)*species(48);
    ratelaws(96) = rateLaw11(expressions,observables)*species(50)*species(5);
    ratelaws(97) = rateLaw12(expressions,observables)*species(49)*species(5);
    ratelaws(98) = rateLaw13(expressions,observables)*species(52)*species(5);
    ratelaws(99) = rateLaw14(expressions,observables)*species(51)*species(5);
    ratelaws(100) = rateLaw15(expressions,observables)*species(54)*species(5);
    ratelaws(101) = rateLaw16(expressions,observables)*species(53)*species(5);
    ratelaws(102) = rateLaw17(expressions,observables)*species(50)*species(7);
    ratelaws(103) = rateLaw18(expressions,observables)*species(49)*species(7);
    ratelaws(104) = rateLaw19(expressions,observables)*species(52)*species(7);
    ratelaws(105) = rateLaw20(expressions,observables)*species(51)*species(7);
    ratelaws(106) = rateLaw21(expressions,observables)*species(54)*species(7);
    ratelaws(107) = rateLaw22(expressions,observables)*species(53)*species(7);
    ratelaws(108) = expressions(63)*species(55);
    ratelaws(109) = expressions(64)*species(55);
    ratelaws(110) = 2*expressions(63)*species(56);
    ratelaws(111) = 2*expressions(66)*species(56);
    ratelaws(112) = expressions(63)*species(57);
    ratelaws(113) = expressions(67)*species(57);
    ratelaws(114) = expressions(63)*species(58);
    ratelaws(115) = (expressions(64)/expressions(30))*species(58);
    ratelaws(116) = 2*expressions(63)*species(59);
    ratelaws(117) = 2*(expressions(66)/expressions(30))*species(59);
    ratelaws(118) = expressions(63)*species(60);
    ratelaws(119) = (expressions(67)/expressions(30))*species(60);
    ratelaws(120) = rateLaw28(expressions,observables)*species(62);
    ratelaws(121) = rateLaw28(expressions,observables)*species(64);
    ratelaws(122) = rateLaw28(expressions,observables)*species(65);
    ratelaws(123) = expressions(100)*species(61)*species(15);
    ratelaws(124) = (expressions(101)*expressions(100))*species(62);
    ratelaws(125) = expressions(103)*species(62)*species(14);
    ratelaws(126) = rateLaw39(expressions,observables)*species(47)*species(17);
    ratelaws(127) = rateLaw40(expressions,observables)*species(47)*species(17);
    ratelaws(128) = rateLaw41(expressions,observables)*species(48)*species(17);
    ratelaws(129) = rateLaw55(expressions,observables)*species(61)*species(17);
    ratelaws(130) = rateLaw56(expressions,observables)*species(61)*species(17);
    ratelaws(131) = expressions(106)*species(61)*species(18);
    ratelaws(132) = expressions(106)*species(61)*species(46);
    ratelaws(133) = (expressions(107)*expressions(106))*species(64);
    ratelaws(134) = (expressions(107)*expressions(106))*species(65);
    ratelaws(135) = expressions(112)*species(64)*species(19);
    ratelaws(136) = (expressions(113)*expressions(112))*species(65);
    ratelaws(137) = expressions(115)*species(65)*species(21);
    ratelaws(138) = expressions(118)*species(65)*species(63);
    ratelaws(139) = expressions(121)*species(65)*species(21);
    ratelaws(140) = expressions(133)*species(63)*species(20);
    ratelaws(141) = rateLaw65(expressions,observables)*species(63)*species(22);
    ratelaws(142) = expressions(171)*species(49)*species(29);
    ratelaws(143) = expressions(171)*species(50)*species(29);
    ratelaws(144) = expressions(171)*species(51)*species(29);
    ratelaws(145) = expressions(171)*species(52)*species(29);
    ratelaws(146) = expressions(171)*species(53)*species(29);
    ratelaws(147) = expressions(171)*species(54)*species(29);
    ratelaws(148) = expressions(26)*species(66)*species(1);
    ratelaws(149) = expressions(26)*species(67)*species(1);
    ratelaws(150) = (expressions(25)*expressions(26))*species(78);
    ratelaws(151) = (expressions(25)*expressions(26))*species(93);
    ratelaws(152) = (expressions(25)*expressions(26))*species(94);
    ratelaws(153) = expressions(26)*species(68)*species(1);
    ratelaws(154) = expressions(26)*species(69)*species(1);
    ratelaws(155) = (expressions(25)*expressions(26))*species(79);
    ratelaws(156) = (expressions(25)*expressions(26))*species(95);
    ratelaws(157) = (expressions(25)*expressions(26))*species(96);
    ratelaws(158) = expressions(26)*species(70)*species(1);
    ratelaws(159) = expressions(26)*species(71)*species(1);
    ratelaws(160) = (expressions(25)*expressions(26))*species(80);
    ratelaws(161) = (expressions(25)*expressions(26))*species(97);
    ratelaws(162) = (expressions(25)*expressions(26))*species(98);
    ratelaws(163) = expressions(26)*species(72)*species(1);
    ratelaws(164) = expressions(26)*species(73)*species(1);
    ratelaws(165) = (expressions(25)*expressions(26))*species(81);
    ratelaws(166) = expressions(26)*species(74)*species(1);
    ratelaws(167) = expressions(26)*species(75)*species(1);
    ratelaws(168) = (expressions(25)*expressions(26))*species(82);
    ratelaws(169) = expressions(26)*species(76)*species(1);
    ratelaws(170) = expressions(26)*species(77)*species(1);
    ratelaws(171) = (expressions(25)*expressions(26))*species(83);
    ratelaws(172) = expressions(63)*species(72);
    ratelaws(173) = expressions(64)*species(72);
    ratelaws(174) = 2*expressions(63)*species(74);
    ratelaws(175) = 2*expressions(66)*species(74);
    ratelaws(176) = expressions(62)*species(82)*species(5);
    ratelaws(177) = expressions(63)*species(76);
    ratelaws(178) = expressions(67)*species(76);
    ratelaws(179) = expressions(63)*species(73);
    ratelaws(180) = (expressions(64)/expressions(30))*species(73);
    ratelaws(181) = 2*expressions(63)*species(75);
    ratelaws(182) = 2*(expressions(66)/expressions(30))*species(75);
    ratelaws(183) = expressions(62)*species(82)*species(7);
    ratelaws(184) = expressions(63)*species(77);
    ratelaws(185) = (expressions(67)/expressions(30))*species(77);
    ratelaws(186) = expressions(72)*species(78)*species(13);
    ratelaws(187) = expressions(72)*species(79)*species(13);
    ratelaws(188) = expressions(72)*species(80)*species(13);
    ratelaws(189) = rateLaw28(expressions,observables)*species(87);
    ratelaws(190) = rateLaw28(expressions,observables)*species(88);
    ratelaws(191) = rateLaw28(expressions,observables)*species(89);
    ratelaws(192) = rateLaw28(expressions,observables)*species(90);
    ratelaws(193) = (expressions(101)*expressions(100))*species(84);
    ratelaws(194) = expressions(103)*species(84)*species(14);
    ratelaws(195) = (expressions(104)*expressions(103))*species(87);
    ratelaws(196) = rateLaw29(expressions,observables)*species(87);
    ratelaws(197) = expressions(91)*species(81)*species(16);
    ratelaws(198) = rateLaw45(expressions,observables)*species(66)*species(17);
    ratelaws(199) = rateLaw46(expressions,observables)*species(67)*species(17);
    ratelaws(200) = rateLaw47(expressions,observables)*species(68)*species(17);
    ratelaws(201) = rateLaw48(expressions,observables)*species(69)*species(17);
    ratelaws(202) = rateLaw49(expressions,observables)*species(70)*species(17);
    ratelaws(203) = rateLaw50(expressions,observables)*species(71)*species(17);
    ratelaws(204) = rateLaw56(expressions,observables)*species(84)*species(17);
    ratelaws(205) = rateLaw56(expressions,observables)*species(85)*species(17);
    ratelaws(206) = rateLaw56(expressions,observables)*species(86)*species(17);
    ratelaws(207) = (expressions(107)*expressions(106))*species(85);
    ratelaws(208) = (expressions(107)*expressions(106))*species(86);
    ratelaws(209) = (expressions(107)*expressions(106))*species(88);
    ratelaws(210) = (expressions(107)*expressions(106))*species(89);
    ratelaws(211) = (expressions(107)*expressions(106))*species(90);
    ratelaws(212) = expressions(112)*species(85)*species(19);
    ratelaws(213) = (expressions(113)*expressions(112))*species(86);
    ratelaws(214) = expressions(115)*species(86)*species(21);
    ratelaws(215) = (expressions(116)*expressions(115))*species(88);
    ratelaws(216) = expressions(118)*species(86)*species(63);
    ratelaws(217) = (expressions(119)*expressions(118))*species(89);
    ratelaws(218) = expressions(121)*species(86)*species(21);
    ratelaws(219) = (expressions(122)*expressions(121))*species(90);
    ratelaws(220) = expressions(124)*species(90);
    ratelaws(221) = expressions(129)*species(88)*species(21);
    ratelaws(222) = expressions(125)*species(89)*species(21);
    ratelaws(223) = (expressions(134)*expressions(133))*species(91);
    ratelaws(224) = expressions(136)*species(91);
    ratelaws(225) = rateLaw66(expressions,observables)*species(92);
    ratelaws(226) = rateLaw67(expressions,observables)*species(92)*species(23);
    ratelaws(227) = (expressions(172)*expressions(171))*species(93);
    ratelaws(228) = (expressions(172)*expressions(171))*species(94);
    ratelaws(229) = (expressions(172)*expressions(171))*species(95);
    ratelaws(230) = (expressions(172)*expressions(171))*species(96);
    ratelaws(231) = (expressions(172)*expressions(171))*species(97);
    ratelaws(232) = (expressions(172)*expressions(171))*species(98);
    ratelaws(233) = rateLaw75(expressions,observables)*species(93)*species(5);
    ratelaws(234) = rateLaw76(expressions,observables)*species(94)*species(5);
    ratelaws(235) = rateLaw77(expressions,observables)*species(95)*species(5);
    ratelaws(236) = rateLaw78(expressions,observables)*species(96)*species(5);
    ratelaws(237) = rateLaw79(expressions,observables)*species(97)*species(5);
    ratelaws(238) = rateLaw80(expressions,observables)*species(98)*species(5);
    ratelaws(239) = rateLaw81(expressions,observables)*species(93)*species(7);
    ratelaws(240) = rateLaw82(expressions,observables)*species(94)*species(7);
    ratelaws(241) = rateLaw83(expressions,observables)*species(95)*species(7);
    ratelaws(242) = rateLaw84(expressions,observables)*species(96)*species(7);
    ratelaws(243) = rateLaw85(expressions,observables)*species(97)*species(7);
    ratelaws(244) = rateLaw86(expressions,observables)*species(98)*species(7);
    ratelaws(245) = expressions(26)*species(99)*species(1);
    ratelaws(246) = expressions(26)*species(100)*species(1);
    ratelaws(247) = expressions(26)*species(101)*species(1);
    ratelaws(248) = (expressions(25)*expressions(26))*species(113);
    ratelaws(249) = expressions(26)*species(102)*species(1);
    ratelaws(250) = expressions(26)*species(103)*species(1);
    ratelaws(251) = expressions(26)*species(104)*species(1);
    ratelaws(252) = (expressions(25)*expressions(26))*species(114);
    ratelaws(253) = expressions(26)*species(105)*species(1);
    ratelaws(254) = expressions(26)*species(106)*species(1);
    ratelaws(255) = expressions(26)*species(107)*species(1);
    ratelaws(256) = (expressions(25)*expressions(26))*species(115);
    ratelaws(257) = expressions(26)*species(108)*species(1);
    ratelaws(258) = (expressions(25)*expressions(26))*species(121);
    ratelaws(259) = expressions(26)*species(109)*species(1);
    ratelaws(260) = (expressions(25)*expressions(26))*species(111);
    ratelaws(261) = (expressions(25)*expressions(26))*species(112);
    ratelaws(262) = expressions(26)*species(110)*species(1);
    ratelaws(263) = expressions(63)*species(111);
    ratelaws(264) = expressions(65)*species(111);
    ratelaws(265) = expressions(63)*species(112);
    ratelaws(266) = (expressions(65)/expressions(30))*species(112);
    ratelaws(267) = expressions(72)*species(99)*species(13);
    ratelaws(268) = (expressions(72)*expressions(73))*species(113);
    ratelaws(269) = expressions(72)*species(102)*species(13);
    ratelaws(270) = (expressions(72)*expressions(73))*species(114);
    ratelaws(271) = expressions(72)*species(105)*species(13);
    ratelaws(272) = (expressions(72)*expressions(73))*species(115);
    ratelaws(273) = rateLaw23(expressions,observables)*species(113)*species(7);
    ratelaws(274) = rateLaw23(expressions,observables)*species(114)*species(7);
    ratelaws(275) = rateLaw23(expressions,observables)*species(115)*species(7);
    ratelaws(276) = rateLaw24(expressions,observables)*species(113)*species(5);
    ratelaws(277) = rateLaw24(expressions,observables)*species(114)*species(5);
    ratelaws(278) = rateLaw24(expressions,observables)*species(115)*species(5);
    ratelaws(279) = rateLaw28(expressions,observables)*species(120);
    ratelaws(280) = rateLaw28(expressions,observables)*species(125);
    ratelaws(281) = rateLaw28(expressions,observables)*species(126);
    ratelaws(282) = (expressions(104)*expressions(103))*species(116);
    ratelaws(283) = (expressions(104)*expressions(103))*species(120);
    ratelaws(284) = rateLaw29(expressions,observables)*species(116);
    ratelaws(285) = expressions(91)*species(108)*species(16);
    ratelaws(286) = (expressions(92)*expressions(91))*species(121);
    ratelaws(287) = rateLaw45(expressions,observables)*species(99)*species(17);
    ratelaws(288) = rateLaw46(expressions,observables)*species(99)*species(17);
    ratelaws(289) = rateLaw47(expressions,observables)*species(102)*species(17);
    ratelaws(290) = rateLaw48(expressions,observables)*species(102)*species(17);
    ratelaws(291) = rateLaw49(expressions,observables)*species(105)*species(17);
    ratelaws(292) = rateLaw50(expressions,observables)*species(105)*species(17);
    ratelaws(293) = rateLaw56(expressions,observables)*species(116)*species(17);
    ratelaws(294) = rateLaw56(expressions,observables)*species(117)*species(17);
    ratelaws(295) = rateLaw56(expressions,observables)*species(118)*species(17);
    ratelaws(296) = rateLaw56(expressions,observables)*species(119)*species(17);
    ratelaws(297) = rateLaw58(expressions,observables)*species(108)*species(17);
    ratelaws(298) = rateLaw61(expressions,observables)*species(109)*species(17);
    ratelaws(299) = rateLaw62(expressions,observables)*species(110)*species(17);
    ratelaws(300) = expressions(137)*species(120)*species(25);
    ratelaws(301) = expressions(106)*species(42)*species(122);
    ratelaws(302) = expressions(106)*species(42)*species(123);
    ratelaws(303) = expressions(106)*species(42)*species(124);
    ratelaws(304) = expressions(106)*species(61)*species(122);
    ratelaws(305) = expressions(106)*species(61)*species(123);
    ratelaws(306) = expressions(106)*species(61)*species(124);
    ratelaws(307) = (expressions(107)*expressions(106))*species(117);
    ratelaws(308) = (expressions(107)*expressions(106))*species(118);
    ratelaws(309) = (expressions(107)*expressions(106))*species(119);
    ratelaws(310) = (expressions(107)*expressions(106))*species(125);
    ratelaws(311) = (expressions(107)*expressions(106))*species(126);
    ratelaws(312) = (expressions(116)*expressions(115))*species(117);
    ratelaws(313) = (expressions(119)*expressions(118))*species(118);
    ratelaws(314) = (expressions(122)*expressions(121))*species(119);
    ratelaws(315) = expressions(124)*species(119);
    ratelaws(316) = expressions(129)*species(117)*species(21);
    ratelaws(317) = (expressions(130)*expressions(129))*species(125);
    ratelaws(318) = expressions(132)*species(125);
    ratelaws(319) = expressions(125)*species(118)*species(21);
    ratelaws(320) = (expressions(126)*expressions(125))*species(126);
    ratelaws(321) = expressions(128)*species(126);
    ratelaws(322) = rateLaw68(expressions,observables)*species(92)*species(127);
    ratelaws(323) = rateLaw70(expressions,observables)*species(127);
    ratelaws(324) = (expressions(172)*expressions(171))*species(100);
    ratelaws(325) = (expressions(172)*expressions(171))*species(101);
    ratelaws(326) = (expressions(172)*expressions(171))*species(103);
    ratelaws(327) = (expressions(172)*expressions(171))*species(104);
    ratelaws(328) = (expressions(172)*expressions(171))*species(106);
    ratelaws(329) = (expressions(172)*expressions(171))*species(107);
    ratelaws(330) = rateLaw75(expressions,observables)*species(100)*species(5);
    ratelaws(331) = rateLaw76(expressions,observables)*species(101)*species(5);
    ratelaws(332) = rateLaw77(expressions,observables)*species(103)*species(5);
    ratelaws(333) = rateLaw78(expressions,observables)*species(104)*species(5);
    ratelaws(334) = rateLaw79(expressions,observables)*species(106)*species(5);
    ratelaws(335) = rateLaw80(expressions,observables)*species(107)*species(5);
    ratelaws(336) = rateLaw81(expressions,observables)*species(100)*species(7);
    ratelaws(337) = rateLaw82(expressions,observables)*species(101)*species(7);
    ratelaws(338) = rateLaw83(expressions,observables)*species(103)*species(7);
    ratelaws(339) = rateLaw84(expressions,observables)*species(104)*species(7);
    ratelaws(340) = rateLaw85(expressions,observables)*species(106)*species(7);
    ratelaws(341) = rateLaw86(expressions,observables)*species(107)*species(7);
    ratelaws(342) = rateLaw88(expressions,observables)*species(82)*species(128);
    ratelaws(343) = rateLaw88(expressions,observables)*species(109)*species(128);
    ratelaws(344) = rateLaw89(expressions,observables)*species(5)*species(128);
    ratelaws(345) = rateLaw89(expressions,observables)*species(7)*species(128);
    ratelaws(346) = rateLaw90(expressions,observables)*species(7)*species(128);
    ratelaws(347) = rateLaw90(expressions,observables)*species(44)*species(128);
    ratelaws(348) = rateLaw91(expressions,observables)*species(49)*species(128);
    ratelaws(349) = rateLaw91(expressions,observables)*species(66)*species(128);
    ratelaws(350) = rateLaw91(expressions,observables)*species(78)*species(128);
    ratelaws(351) = rateLaw91(expressions,observables)*species(99)*species(128);
    ratelaws(352) = rateLaw92(expressions,observables)*species(50)*species(128);
    ratelaws(353) = rateLaw92(expressions,observables)*species(67)*species(128);
    ratelaws(354) = rateLaw92(expressions,observables)*species(78)*species(128);
    ratelaws(355) = rateLaw92(expressions,observables)*species(99)*species(128);
    ratelaws(356) = rateLaw93(expressions,observables)*species(51)*species(128);
    ratelaws(357) = rateLaw93(expressions,observables)*species(68)*species(128);
    ratelaws(358) = rateLaw93(expressions,observables)*species(79)*species(128);
    ratelaws(359) = rateLaw93(expressions,observables)*species(102)*species(128);
    ratelaws(360) = rateLaw94(expressions,observables)*species(52)*species(128);
    ratelaws(361) = rateLaw94(expressions,observables)*species(69)*species(128);
    ratelaws(362) = rateLaw94(expressions,observables)*species(79)*species(128);
    ratelaws(363) = rateLaw94(expressions,observables)*species(102)*species(128);
    ratelaws(364) = rateLaw95(expressions,observables)*species(53)*species(128);
    ratelaws(365) = rateLaw95(expressions,observables)*species(70)*species(128);
    ratelaws(366) = rateLaw95(expressions,observables)*species(80)*species(128);
    ratelaws(367) = rateLaw95(expressions,observables)*species(105)*species(128);
    ratelaws(368) = rateLaw96(expressions,observables)*species(54)*species(128);
    ratelaws(369) = rateLaw96(expressions,observables)*species(71)*species(128);
    ratelaws(370) = rateLaw96(expressions,observables)*species(80)*species(128);
    ratelaws(371) = rateLaw96(expressions,observables)*species(105)*species(128);
    ratelaws(372) = rateLaw99(expressions,observables)*species(128)*species(128);
    ratelaws(373) = expressions(26)*species(129)*species(1);
    ratelaws(374) = (expressions(25)*expressions(26))*species(136);
    ratelaws(375) = expressions(26)*species(130)*species(1);
    ratelaws(376) = (expressions(25)*expressions(26))*species(137);
    ratelaws(377) = expressions(26)*species(131)*species(1);
    ratelaws(378) = (expressions(25)*expressions(26))*species(138);
    ratelaws(379) = expressions(26)*species(132)*species(1);
    ratelaws(380) = expressions(26)*species(133)*species(1);
    ratelaws(381) = expressions(26)*species(134)*species(1);
    ratelaws(382) = (expressions(25)*expressions(26))*species(135);
    ratelaws(383) = expressions(63)*species(133);
    ratelaws(384) = expressions(65)*species(133);
    ratelaws(385) = expressions(63)*species(134);
    ratelaws(386) = (expressions(65)/expressions(30))*species(134);
    ratelaws(387) = (expressions(72)*expressions(73))*species(129);
    ratelaws(388) = (expressions(72)*expressions(73))*species(130);
    ratelaws(389) = (expressions(72)*expressions(73))*species(131);
    ratelaws(390) = (expressions(72)*expressions(73))*species(136);
    ratelaws(391) = (expressions(72)*expressions(73))*species(137);
    ratelaws(392) = (expressions(72)*expressions(73))*species(138);
    ratelaws(393) = rateLaw23(expressions,observables)*species(129)*species(7);
    ratelaws(394) = rateLaw23(expressions,observables)*species(130)*species(7);
    ratelaws(395) = rateLaw23(expressions,observables)*species(131)*species(7);
    ratelaws(396) = rateLaw24(expressions,observables)*species(129)*species(5);
    ratelaws(397) = rateLaw24(expressions,observables)*species(130)*species(5);
    ratelaws(398) = rateLaw24(expressions,observables)*species(131)*species(5);
    ratelaws(399) = rateLaw25(expressions,observables)*species(136)*species(7);
    ratelaws(400) = rateLaw25(expressions,observables)*species(137)*species(7);
    ratelaws(401) = rateLaw25(expressions,observables)*species(138)*species(7);
    ratelaws(402) = rateLaw26(expressions,observables)*species(136)*species(5);
    ratelaws(403) = rateLaw26(expressions,observables)*species(137)*species(5);
    ratelaws(404) = rateLaw26(expressions,observables)*species(138)*species(5);
    ratelaws(405) = rateLaw28(expressions,observables)*species(143);
    ratelaws(406) = expressions(103)*species(62)*species(142);
    ratelaws(407) = expressions(103)*species(84)*species(142);
    ratelaws(408) = (expressions(104)*expressions(103))*species(139);
    ratelaws(409) = (expressions(92)*expressions(91))*species(132);
    ratelaws(410) = expressions(94)*species(135)*species(15);
    ratelaws(411) = rateLaw51(expressions,observables)*species(136)*species(17);
    ratelaws(412) = rateLaw51(expressions,observables)*species(137)*species(17);
    ratelaws(413) = rateLaw51(expressions,observables)*species(138)*species(17);
    ratelaws(414) = rateLaw56(expressions,observables)*species(139)*species(17);
    ratelaws(415) = rateLaw56(expressions,observables)*species(140)*species(17);
    ratelaws(416) = rateLaw56(expressions,observables)*species(141)*species(17);
    ratelaws(417) = rateLaw57(expressions,observables)*species(142)*species(17);
    ratelaws(418) = expressions(137)*species(139)*species(25);
    ratelaws(419) = (expressions(138)*expressions(137))*species(143);
    ratelaws(420) = rateLaw63(expressions,observables)*species(143)*species(27);
    ratelaws(421) = expressions(106)*species(42)*species(144);
    ratelaws(422) = expressions(106)*species(42)*species(145);
    ratelaws(423) = expressions(106)*species(61)*species(144);
    ratelaws(424) = expressions(106)*species(61)*species(145);
    ratelaws(425) = (expressions(107)*expressions(106))*species(140);
    ratelaws(426) = (expressions(107)*expressions(106))*species(141);
    ratelaws(427) = expressions(97)*species(135)*species(18);
    ratelaws(428) = expressions(97)*species(135)*species(46);
    ratelaws(429) = expressions(97)*species(135)*species(122);
    ratelaws(430) = expressions(97)*species(135)*species(123);
    ratelaws(431) = expressions(97)*species(135)*species(124);
    ratelaws(432) = expressions(97)*species(135)*species(144);
    ratelaws(433) = expressions(97)*species(135)*species(145);
    ratelaws(434) = (expressions(130)*expressions(129))*species(140);
    ratelaws(435) = expressions(132)*species(140);
    ratelaws(436) = (expressions(126)*expressions(125))*species(141);
    ratelaws(437) = expressions(128)*species(141);
    ratelaws(438) = rateLaw69(expressions,observables)*species(146);
    ratelaws(439) = rateLaw71(expressions,observables)*species(146)*species(24);
    ratelaws(440) = rateLaw87(expressions,observables)*species(135)*species(128);
    ratelaws(441) = rateLaw88(expressions,observables)*species(135)*species(128);
    ratelaws(442) = rateLaw97(expressions,observables)*species(136)*species(128);
    ratelaws(443) = rateLaw97(expressions,observables)*species(137)*species(128);
    ratelaws(444) = rateLaw97(expressions,observables)*species(138)*species(128);
    ratelaws(445) = expressions(26)*species(147)*species(1);
    ratelaws(446) = (expressions(25)*expressions(26))*species(152);
    ratelaws(447) = expressions(26)*species(148)*species(1);
    ratelaws(448) = (expressions(25)*expressions(26))*species(153);
    ratelaws(449) = expressions(26)*species(149)*species(1);
    ratelaws(450) = (expressions(25)*expressions(26))*species(154);
    ratelaws(451) = expressions(26)*species(150)*species(1);
    ratelaws(452) = (expressions(25)*expressions(26))*species(156);
    ratelaws(453) = (expressions(25)*expressions(26))*species(157);
    ratelaws(454) = (expressions(25)*expressions(26))*species(158);
    ratelaws(455) = (expressions(25)*expressions(26))*species(159);
    ratelaws(456) = (expressions(25)*expressions(26))*species(160);
    ratelaws(457) = (expressions(25)*expressions(26))*species(161);
    ratelaws(458) = (expressions(25)*expressions(26))*species(162);
    ratelaws(459) = (expressions(25)*expressions(26))*species(163);
    ratelaws(460) = (expressions(25)*expressions(26))*species(165);
    ratelaws(461) = expressions(72)*species(78)*species(151);
    ratelaws(462) = expressions(72)*species(99)*species(151);
    ratelaws(463) = (expressions(72)*expressions(73))*species(147);
    ratelaws(464) = expressions(72)*species(79)*species(151);
    ratelaws(465) = expressions(72)*species(102)*species(151);
    ratelaws(466) = (expressions(72)*expressions(73))*species(148);
    ratelaws(467) = expressions(72)*species(80)*species(151);
    ratelaws(468) = expressions(72)*species(105)*species(151);
    ratelaws(469) = (expressions(72)*expressions(73))*species(149);
    ratelaws(470) = (expressions(75)*expressions(76))*species(152);
    ratelaws(471) = (expressions(75)*expressions(76))*species(153);
    ratelaws(472) = (expressions(75)*expressions(76))*species(154);
    ratelaws(473) = rateLaw25(expressions,observables)*species(147)*species(7);
    ratelaws(474) = rateLaw25(expressions,observables)*species(148)*species(7);
    ratelaws(475) = rateLaw25(expressions,observables)*species(149)*species(7);
    ratelaws(476) = rateLaw26(expressions,observables)*species(147)*species(5);
    ratelaws(477) = rateLaw26(expressions,observables)*species(148)*species(5);
    ratelaws(478) = rateLaw26(expressions,observables)*species(149)*species(5);
    ratelaws(479) = expressions(103)*species(156)*species(14);
    ratelaws(480) = expressions(103)*species(156)*species(142);
    ratelaws(481) = expressions(94)*species(150)*species(15);
    ratelaws(482) = (expressions(95)*expressions(94))*species(156);
    ratelaws(483) = rateLaw51(expressions,observables)*species(147)*species(17);
    ratelaws(484) = rateLaw51(expressions,observables)*species(148)*species(17);
    ratelaws(485) = rateLaw51(expressions,observables)*species(149)*species(17);
    ratelaws(486) = rateLaw51(expressions,observables)*species(152)*species(17);
    ratelaws(487) = rateLaw51(expressions,observables)*species(153)*species(17);
    ratelaws(488) = rateLaw51(expressions,observables)*species(154)*species(17);
    ratelaws(489) = rateLaw52(expressions,observables)*species(152)*species(17);
    ratelaws(490) = rateLaw52(expressions,observables)*species(153)*species(17);
    ratelaws(491) = rateLaw52(expressions,observables)*species(154)*species(17);
    ratelaws(492) = rateLaw53(expressions,observables)*species(151)*species(17);
    ratelaws(493) = rateLaw56(expressions,observables)*species(155)*species(17);
    ratelaws(494) = rateLaw59(expressions,observables)*species(150)*species(17);
    ratelaws(495) = rateLaw61(expressions,observables)*species(150)*species(17);
    ratelaws(496) = (expressions(138)*expressions(137))*species(155);
    ratelaws(497) = rateLaw63(expressions,observables)*species(155)*species(27);
    ratelaws(498) = expressions(97)*species(150)*species(18);
    ratelaws(499) = expressions(97)*species(150)*species(46);
    ratelaws(500) = expressions(97)*species(150)*species(122);
    ratelaws(501) = expressions(97)*species(150)*species(123);
    ratelaws(502) = expressions(97)*species(150)*species(124);
    ratelaws(503) = expressions(97)*species(150)*species(144);
    ratelaws(504) = expressions(97)*species(150)*species(145);
    ratelaws(505) = (expressions(98)*expressions(97))*species(157);
    ratelaws(506) = (expressions(98)*expressions(97))*species(158);
    ratelaws(507) = (expressions(98)*expressions(97))*species(159);
    ratelaws(508) = (expressions(98)*expressions(97))*species(160);
    ratelaws(509) = (expressions(98)*expressions(97))*species(161);
    ratelaws(510) = (expressions(98)*expressions(97))*species(162);
    ratelaws(511) = (expressions(98)*expressions(97))*species(163);
    ratelaws(512) = expressions(112)*species(157)*species(19);
    ratelaws(513) = (expressions(113)*expressions(112))*species(158);
    ratelaws(514) = expressions(115)*species(158)*species(21);
    ratelaws(515) = (expressions(116)*expressions(115))*species(159);
    ratelaws(516) = expressions(118)*species(158)*species(63);
    ratelaws(517) = (expressions(119)*expressions(118))*species(160);
    ratelaws(518) = expressions(121)*species(158)*species(21);
    ratelaws(519) = (expressions(122)*expressions(121))*species(161);
    ratelaws(520) = expressions(124)*species(161);
    ratelaws(521) = expressions(129)*species(159)*species(21);
    ratelaws(522) = (expressions(130)*expressions(129))*species(162);
    ratelaws(523) = expressions(132)*species(162);
    ratelaws(524) = expressions(125)*species(160)*species(21);
    ratelaws(525) = (expressions(126)*expressions(125))*species(163);
    ratelaws(526) = expressions(128)*species(163);
    ratelaws(527) = rateLaw72(expressions,observables)*species(146)*species(164);
    ratelaws(528) = rateLaw74(expressions,observables)*species(164);
    ratelaws(529) = rateLaw87(expressions,observables)*species(150)*species(128);
    ratelaws(530) = rateLaw87(expressions,observables)*species(165)*species(128);
    ratelaws(531) = rateLaw88(expressions,observables)*species(150)*species(128);
    ratelaws(532) = rateLaw97(expressions,observables)*species(147)*species(128);
    ratelaws(533) = rateLaw97(expressions,observables)*species(148)*species(128);
    ratelaws(534) = rateLaw97(expressions,observables)*species(149)*species(128);
    ratelaws(535) = rateLaw97(expressions,observables)*species(151)*species(128);
    ratelaws(536) = rateLaw97(expressions,observables)*species(152)*species(128);
    ratelaws(537) = rateLaw97(expressions,observables)*species(153)*species(128);
    ratelaws(538) = rateLaw97(expressions,observables)*species(154)*species(128);
    ratelaws(539) = rateLaw98(expressions,observables)*species(152)*species(128);
    ratelaws(540) = rateLaw98(expressions,observables)*species(153)*species(128);
    ratelaws(541) = rateLaw98(expressions,observables)*species(154)*species(128);
    ratelaws(542) = expressions(26)*species(166)*species(1);
    ratelaws(543) = (expressions(25)*expressions(26))*species(181);
    ratelaws(544) = expressions(26)*species(167)*species(1);
    ratelaws(545) = (expressions(25)*expressions(26))*species(182);
    ratelaws(546) = expressions(26)*species(168)*species(1);
    ratelaws(547) = (expressions(25)*expressions(26))*species(183);
    ratelaws(548) = expressions(26)*species(169)*species(1);
    ratelaws(549) = expressions(26)*species(170)*species(1);
    ratelaws(550) = expressions(26)*species(171)*species(1);
    ratelaws(551) = expressions(26)*species(172)*species(1);
    ratelaws(552) = expressions(26)*species(173)*species(1);
    ratelaws(553) = expressions(26)*species(174)*species(1);
    ratelaws(554) = expressions(26)*species(175)*species(1);
    ratelaws(555) = expressions(26)*species(176)*species(1);
    ratelaws(556) = expressions(26)*species(177)*species(1);
    ratelaws(557) = (expressions(25)*expressions(26))*species(179);
    ratelaws(558) = (expressions(25)*expressions(26))*species(180);
    ratelaws(559) = (expressions(72)*expressions(73))*species(181);
    ratelaws(560) = (expressions(72)*expressions(73))*species(182);
    ratelaws(561) = (expressions(72)*expressions(73))*species(183);
    ratelaws(562) = expressions(75)*species(78)*species(178);
    ratelaws(563) = expressions(75)*species(99)*species(178);
    ratelaws(564) = (expressions(75)*expressions(76))*species(166);
    ratelaws(565) = expressions(75)*species(79)*species(178);
    ratelaws(566) = expressions(75)*species(102)*species(178);
    ratelaws(567) = (expressions(75)*expressions(76))*species(167);
    ratelaws(568) = expressions(75)*species(80)*species(178);
    ratelaws(569) = expressions(75)*species(105)*species(178);
    ratelaws(570) = (expressions(75)*expressions(76))*species(168);
    ratelaws(571) = rateLaw23(expressions,observables)*species(181)*species(7);
    ratelaws(572) = rateLaw23(expressions,observables)*species(182)*species(7);
    ratelaws(573) = rateLaw23(expressions,observables)*species(183)*species(7);
    ratelaws(574) = rateLaw24(expressions,observables)*species(181)*species(5);
    ratelaws(575) = rateLaw24(expressions,observables)*species(182)*species(5);
    ratelaws(576) = rateLaw24(expressions,observables)*species(183)*species(5);
    ratelaws(577) = expressions(103)*species(169)*species(14);
    ratelaws(578) = expressions(103)*species(169)*species(142);
    ratelaws(579) = (expressions(104)*expressions(103))*species(179);
    ratelaws(580) = (expressions(104)*expressions(103))*species(180);
    ratelaws(581) = (expressions(95)*expressions(94))*species(169);
    ratelaws(582) = rateLaw51(expressions,observables)*species(166)*species(17);
    ratelaws(583) = rateLaw51(expressions,observables)*species(167)*species(17);
    ratelaws(584) = rateLaw51(expressions,observables)*species(168)*species(17);
    ratelaws(585) = rateLaw52(expressions,observables)*species(166)*species(17);
    ratelaws(586) = rateLaw52(expressions,observables)*species(167)*species(17);
    ratelaws(587) = rateLaw52(expressions,observables)*species(168)*species(17);
    ratelaws(588) = rateLaw52(expressions,observables)*species(181)*species(17);
    ratelaws(589) = rateLaw52(expressions,observables)*species(182)*species(17);
    ratelaws(590) = rateLaw52(expressions,observables)*species(183)*species(17);
    ratelaws(591) = rateLaw53(expressions,observables)*species(178)*species(17);
    ratelaws(592) = rateLaw54(expressions,observables)*species(178)*species(17);
    ratelaws(593) = rateLaw60(expressions,observables)*species(177)*species(17);
    ratelaws(594) = (expressions(98)*expressions(97))*species(170);
    ratelaws(595) = (expressions(98)*expressions(97))*species(171);
    ratelaws(596) = (expressions(98)*expressions(97))*species(172);
    ratelaws(597) = (expressions(98)*expressions(97))*species(173);
    ratelaws(598) = (expressions(98)*expressions(97))*species(174);
    ratelaws(599) = (expressions(98)*expressions(97))*species(175);
    ratelaws(600) = (expressions(98)*expressions(97))*species(176);
    ratelaws(601) = expressions(112)*species(170)*species(19);
    ratelaws(602) = (expressions(113)*expressions(112))*species(171);
    ratelaws(603) = expressions(115)*species(171)*species(21);
    ratelaws(604) = (expressions(116)*expressions(115))*species(172);
    ratelaws(605) = expressions(118)*species(171)*species(63);
    ratelaws(606) = (expressions(119)*expressions(118))*species(173);
    ratelaws(607) = expressions(121)*species(171)*species(21);
    ratelaws(608) = (expressions(122)*expressions(121))*species(174);
    ratelaws(609) = expressions(124)*species(174);
    ratelaws(610) = expressions(129)*species(172)*species(21);
    ratelaws(611) = (expressions(130)*expressions(129))*species(175);
    ratelaws(612) = expressions(132)*species(175);
    ratelaws(613) = expressions(125)*species(173)*species(21);
    ratelaws(614) = (expressions(126)*expressions(125))*species(176);
    ratelaws(615) = expressions(128)*species(176);
    ratelaws(616) = rateLaw73(expressions,observables)*species(184);
    ratelaws(617) = rateLaw87(expressions,observables)*species(177)*species(128);
    ratelaws(618) = rateLaw97(expressions,observables)*species(166)*species(128);
    ratelaws(619) = rateLaw97(expressions,observables)*species(167)*species(128);
    ratelaws(620) = rateLaw97(expressions,observables)*species(168)*species(128);
    ratelaws(621) = rateLaw97(expressions,observables)*species(178)*species(128);
    ratelaws(622) = rateLaw98(expressions,observables)*species(166)*species(128);
    ratelaws(623) = rateLaw98(expressions,observables)*species(167)*species(128);
    ratelaws(624) = rateLaw98(expressions,observables)*species(168)*species(128);
    ratelaws(625) = rateLaw98(expressions,observables)*species(178)*species(128);
    ratelaws(626) = rateLaw98(expressions,observables)*species(181)*species(128);
    ratelaws(627) = rateLaw98(expressions,observables)*species(182)*species(128);
    ratelaws(628) = rateLaw98(expressions,observables)*species(183)*species(128);
    ratelaws(629) = rateLaw100(expressions,observables)*species(184)*species(5);
    ratelaws(630) = rateLaw100(expressions,observables)*species(184)*species(6);
    ratelaws(631) = rateLaw100(expressions,observables)*species(184)*species(7);
    ratelaws(632) = rateLaw100(expressions,observables)*species(184)*species(44);
    ratelaws(633) = expressions(26)*species(185)*species(1);
    ratelaws(634) = expressions(26)*species(186)*species(1);
    ratelaws(635) = expressions(26)*species(187)*species(1);
    ratelaws(636) = expressions(26)*species(188)*species(1);
    ratelaws(637) = expressions(26)*species(189)*species(1);
    ratelaws(638) = rateLaw2(expressions,observables)*species(192);
    ratelaws(639) = rateLaw3(expressions,observables)*species(191);
    ratelaws(640) = rateLaw4(expressions,observables)*species(194);
    ratelaws(641) = rateLaw5(expressions,observables)*species(192);
    ratelaws(642) = rateLaw6(expressions,observables)*species(191);
    ratelaws(643) = rateLaw7(expressions,observables)*species(194);
    ratelaws(644) = rateLaw8(expressions,observables)*species(192);
    ratelaws(645) = rateLaw9(expressions,observables)*species(191);
    ratelaws(646) = rateLaw10(expressions,observables)*species(194);
    ratelaws(647) = expressions(42)*species(5)*species(193);
    ratelaws(648) = expressions(42)*species(191)*species(7);
    ratelaws(649) = expressions(42)*species(191)*species(193);
    ratelaws(650) = expressions(42)*species(6)*species(193);
    ratelaws(651) = expressions(42)*species(192)*species(7);
    ratelaws(652) = expressions(42)*species(192)*species(193);
    ratelaws(653) = expressions(42)*species(5)*species(194);
    ratelaws(654) = expressions(42)*species(191)*species(44);
    ratelaws(655) = expressions(42)*species(191)*species(194);
    ratelaws(656) = expressions(42)*species(6)*species(194);
    ratelaws(657) = expressions(42)*species(192)*species(44);
    ratelaws(658) = expressions(42)*species(192)*species(194);
    ratelaws(659) = expressions(34)*species(191)*species(8);
    ratelaws(660) = expressions(34)*species(192)*species(8);
    ratelaws(661) = expressions(34)*species(193)*species(8);
    ratelaws(662) = rateLaw11(expressions,observables)*species(30)*species(191);
    ratelaws(663) = rateLaw11(expressions,observables)*species(50)*species(191);
    ratelaws(664) = rateLaw12(expressions,observables)*species(30)*species(191);
    ratelaws(665) = rateLaw12(expressions,observables)*species(49)*species(191);
    ratelaws(666) = rateLaw13(expressions,observables)*species(31)*species(191);
    ratelaws(667) = rateLaw13(expressions,observables)*species(52)*species(191);
    ratelaws(668) = rateLaw14(expressions,observables)*species(31)*species(191);
    ratelaws(669) = rateLaw14(expressions,observables)*species(51)*species(191);
    ratelaws(670) = rateLaw15(expressions,observables)*species(32)*species(191);
    ratelaws(671) = rateLaw15(expressions,observables)*species(54)*species(191);
    ratelaws(672) = rateLaw16(expressions,observables)*species(32)*species(191);
    ratelaws(673) = rateLaw16(expressions,observables)*species(53)*species(191);
    ratelaws(674) = rateLaw17(expressions,observables)*species(30)*species(193);
    ratelaws(675) = rateLaw17(expressions,observables)*species(50)*species(193);
    ratelaws(676) = rateLaw18(expressions,observables)*species(30)*species(193);
    ratelaws(677) = rateLaw18(expressions,observables)*species(49)*species(193);
    ratelaws(678) = rateLaw19(expressions,observables)*species(31)*species(193);
    ratelaws(679) = rateLaw19(expressions,observables)*species(52)*species(193);
    ratelaws(680) = rateLaw20(expressions,observables)*species(31)*species(193);
    ratelaws(681) = rateLaw20(expressions,observables)*species(51)*species(193);
    ratelaws(682) = rateLaw21(expressions,observables)*species(32)*species(193);
    ratelaws(683) = rateLaw21(expressions,observables)*species(54)*species(193);
    ratelaws(684) = rateLaw22(expressions,observables)*species(32)*species(193);
    ratelaws(685) = rateLaw22(expressions,observables)*species(53)*species(193);
    ratelaws(686) = expressions(62)*species(33)*species(191);
    ratelaws(687) = 2*expressions(62)*species(34)*species(191);
    ratelaws(688) = expressions(62)*species(82)*species(191);
    ratelaws(689) = expressions(62)*species(35)*species(191);
    ratelaws(690) = expressions(62)*species(33)*species(193);
    ratelaws(691) = 2*expressions(62)*species(34)*species(193);
    ratelaws(692) = expressions(62)*species(82)*species(193);
    ratelaws(693) = expressions(62)*species(35)*species(193);
    ratelaws(694) = expressions(72)*species(78)*species(190);
    ratelaws(695) = expressions(72)*species(99)*species(190);
    ratelaws(696) = (expressions(72)*expressions(73))*species(185);
    ratelaws(697) = expressions(72)*species(79)*species(190);
    ratelaws(698) = expressions(72)*species(102)*species(190);
    ratelaws(699) = (expressions(72)*expressions(73))*species(186);
    ratelaws(700) = expressions(72)*species(80)*species(190);
    ratelaws(701) = expressions(72)*species(105)*species(190);
    ratelaws(702) = (expressions(72)*expressions(73))*species(187);
    ratelaws(703) = rateLaw23(expressions,observables)*species(113)*species(193);
    ratelaws(704) = rateLaw23(expressions,observables)*species(114)*species(193);
    ratelaws(705) = rateLaw23(expressions,observables)*species(115)*species(193);
    ratelaws(706) = rateLaw23(expressions,observables)*species(129)*species(193);
    ratelaws(707) = rateLaw23(expressions,observables)*species(130)*species(193);
    ratelaws(708) = rateLaw23(expressions,observables)*species(131)*species(193);
    ratelaws(709) = rateLaw23(expressions,observables)*species(181)*species(193);
    ratelaws(710) = rateLaw23(expressions,observables)*species(182)*species(193);
    ratelaws(711) = rateLaw23(expressions,observables)*species(183)*species(193);
    ratelaws(712) = rateLaw23(expressions,observables)*species(185)*species(7);
    ratelaws(713) = rateLaw23(expressions,observables)*species(185)*species(193);
    ratelaws(714) = rateLaw23(expressions,observables)*species(186)*species(7);
    ratelaws(715) = rateLaw23(expressions,observables)*species(186)*species(193);
    ratelaws(716) = rateLaw23(expressions,observables)*species(187)*species(7);
    ratelaws(717) = rateLaw23(expressions,observables)*species(187)*species(193);
    ratelaws(718) = rateLaw24(expressions,observables)*species(113)*species(191);
    ratelaws(719) = rateLaw24(expressions,observables)*species(114)*species(191);
    ratelaws(720) = rateLaw24(expressions,observables)*species(115)*species(191);
    ratelaws(721) = rateLaw24(expressions,observables)*species(129)*species(191);
    ratelaws(722) = rateLaw24(expressions,observables)*species(130)*species(191);
    ratelaws(723) = rateLaw24(expressions,observables)*species(131)*species(191);
    ratelaws(724) = rateLaw24(expressions,observables)*species(181)*species(191);
    ratelaws(725) = rateLaw24(expressions,observables)*species(182)*species(191);
    ratelaws(726) = rateLaw24(expressions,observables)*species(183)*species(191);
    ratelaws(727) = rateLaw24(expressions,observables)*species(185)*species(5);
    ratelaws(728) = rateLaw24(expressions,observables)*species(185)*species(191);
    ratelaws(729) = rateLaw24(expressions,observables)*species(186)*species(5);
    ratelaws(730) = rateLaw24(expressions,observables)*species(186)*species(191);
    ratelaws(731) = rateLaw24(expressions,observables)*species(187)*species(5);
    ratelaws(732) = rateLaw24(expressions,observables)*species(187)*species(191);
    ratelaws(733) = rateLaw25(expressions,observables)*species(136)*species(193);
    ratelaws(734) = rateLaw25(expressions,observables)*species(137)*species(193);
    ratelaws(735) = rateLaw25(expressions,observables)*species(138)*species(193);
    ratelaws(736) = rateLaw25(expressions,observables)*species(147)*species(193);
    ratelaws(737) = rateLaw25(expressions,observables)*species(148)*species(193);
    ratelaws(738) = rateLaw25(expressions,observables)*species(149)*species(193);
    ratelaws(739) = rateLaw26(expressions,observables)*species(136)*species(191);
    ratelaws(740) = rateLaw26(expressions,observables)*species(137)*species(191);
    ratelaws(741) = rateLaw26(expressions,observables)*species(138)*species(191);
    ratelaws(742) = rateLaw26(expressions,observables)*species(147)*species(191);
    ratelaws(743) = rateLaw26(expressions,observables)*species(148)*species(191);
    ratelaws(744) = rateLaw26(expressions,observables)*species(149)*species(191);
    ratelaws(745) = (expressions(104)*expressions(103))*species(188);
    ratelaws(746) = (expressions(104)*expressions(103))*species(189);
    ratelaws(747) = rateLaw30(expressions,observables)*species(193)*species(17);
    ratelaws(748) = rateLaw31(expressions,observables)*species(193)*species(17);
    ratelaws(749) = rateLaw32(expressions,observables)*species(191)*species(17);
    ratelaws(750) = rateLaw33(expressions,observables)*species(194)*species(17);
    ratelaws(751) = rateLaw52(expressions,observables)*species(185)*species(17);
    ratelaws(752) = rateLaw52(expressions,observables)*species(186)*species(17);
    ratelaws(753) = rateLaw52(expressions,observables)*species(187)*species(17);
    ratelaws(754) = rateLaw54(expressions,observables)*species(190)*species(17);
    ratelaws(755) = rateLaw98(expressions,observables)*species(185)*species(128);
    ratelaws(756) = rateLaw98(expressions,observables)*species(186)*species(128);
    ratelaws(757) = rateLaw98(expressions,observables)*species(187)*species(128);
    ratelaws(758) = rateLaw98(expressions,observables)*species(190)*species(128);
    ratelaws(759) = (expressions(25)*expressions(26))*species(210);
    ratelaws(760) = (expressions(25)*expressions(26))*species(214);
    ratelaws(761) = (expressions(25)*expressions(26))*species(211);
    ratelaws(762) = (expressions(25)*expressions(26))*species(212);
    ratelaws(763) = (expressions(25)*expressions(26))*species(215);
    ratelaws(764) = (expressions(25)*expressions(26))*species(216);
    ratelaws(765) = (expressions(25)*expressions(26))*species(213);
    ratelaws(766) = (expressions(25)*expressions(26))*species(217);
    ratelaws(767) = expressions(43)*species(195);
    ratelaws(768) = expressions(43)*species(196);
    ratelaws(769) = expressions(43)*species(197);
    ratelaws(770) = expressions(44)*species(198);
    ratelaws(771) = expressions(44)*species(199);
    ratelaws(772) = expressions(44)*species(200);
    ratelaws(773) = expressions(45)*species(201);
    ratelaws(774) = expressions(45)*species(202);
    ratelaws(775) = expressions(45)*species(203);
    ratelaws(776) = expressions(46)*species(204);
    ratelaws(777) = expressions(46)*species(205);
    ratelaws(778) = expressions(46)*species(206);
    ratelaws(779) = expressions(35)*species(207);
    ratelaws(780) = expressions(37)*species(207);
    ratelaws(781) = expressions(36)*species(208);
    ratelaws(782) = expressions(38)*species(208);
    ratelaws(783) = expressions(35)*species(209);
    ratelaws(784) = expressions(63)*species(210);
    ratelaws(785) = expressions(64)*species(210);
    ratelaws(786) = 2*expressions(63)*species(211);
    ratelaws(787) = 2*expressions(66)*species(211);
    ratelaws(788) = expressions(63)*species(212);
    ratelaws(789) = expressions(65)*species(212);
    ratelaws(790) = expressions(63)*species(213);
    ratelaws(791) = expressions(67)*species(213);
    ratelaws(792) = expressions(63)*species(214);
    ratelaws(793) = (expressions(64)/expressions(30))*species(214);
    ratelaws(794) = 2*expressions(63)*species(215);
    ratelaws(795) = 2*(expressions(66)/expressions(30))*species(215);
    ratelaws(796) = expressions(63)*species(216);
    ratelaws(797) = (expressions(65)/expressions(30))*species(216);
    ratelaws(798) = expressions(63)*species(217);
    ratelaws(799) = (expressions(67)/expressions(30))*species(217);
    ratelaws(800) = rateLaw34(expressions,observables)*species(195)*species(17);
    ratelaws(801) = rateLaw34(expressions,observables)*species(196)*species(17);
    ratelaws(802) = rateLaw34(expressions,observables)*species(197)*species(17);
    ratelaws(803) = rateLaw35(expressions,observables)*species(195)*species(17);
    ratelaws(804) = rateLaw35(expressions,observables)*species(196)*species(17);
    ratelaws(805) = rateLaw35(expressions,observables)*species(197)*species(17);
    ratelaws(806) = rateLaw36(expressions,observables)*species(195)*species(17);
    ratelaws(807) = rateLaw36(expressions,observables)*species(196)*species(17);
    ratelaws(808) = rateLaw36(expressions,observables)*species(197)*species(17);
    ratelaws(809) = rateLaw37(expressions,observables)*species(198)*species(17);
    ratelaws(810) = rateLaw37(expressions,observables)*species(199)*species(17);
    ratelaws(811) = rateLaw37(expressions,observables)*species(200)*species(17);
    ratelaws(812) = rateLaw38(expressions,observables)*species(198)*species(17);
    ratelaws(813) = rateLaw38(expressions,observables)*species(199)*species(17);
    ratelaws(814) = rateLaw38(expressions,observables)*species(200)*species(17);
    ratelaws(815) = rateLaw39(expressions,observables)*species(201)*species(17);
    ratelaws(816) = rateLaw39(expressions,observables)*species(202)*species(17);
    ratelaws(817) = rateLaw39(expressions,observables)*species(203)*species(17);
    ratelaws(818) = rateLaw40(expressions,observables)*species(201)*species(17);
    ratelaws(819) = rateLaw40(expressions,observables)*species(202)*species(17);
    ratelaws(820) = rateLaw40(expressions,observables)*species(203)*species(17);
    ratelaws(821) = rateLaw41(expressions,observables)*species(204)*species(17);
    ratelaws(822) = rateLaw41(expressions,observables)*species(205)*species(17);
    ratelaws(823) = rateLaw41(expressions,observables)*species(206)*species(17);
    ratelaws(824) = rateLaw42(expressions,observables)*species(207)*species(17);
    ratelaws(825) = rateLaw43(expressions,observables)*species(209)*species(17);
    ratelaws(826) = rateLaw44(expressions,observables)*species(209)*species(17);
    ratelaws(827) = expressions(26)*species(218)*species(1);
    ratelaws(828) = expressions(26)*species(219)*species(1);
    ratelaws(829) = expressions(26)*species(220)*species(1);
    ratelaws(830) = expressions(26)*species(221)*species(1);
    ratelaws(831) = expressions(26)*species(222)*species(1);
    ratelaws(832) = expressions(26)*species(223)*species(1);
    ratelaws(833) = expressions(26)*species(224)*species(1);
    ratelaws(834) = expressions(26)*species(225)*species(1);
    ratelaws(835) = expressions(63)*species(218);
    ratelaws(836) = expressions(64)*species(218);
    ratelaws(837) = 2*expressions(63)*species(220);
    ratelaws(838) = 2*expressions(66)*species(220);
    ratelaws(839) = expressions(63)*species(221);
    ratelaws(840) = expressions(65)*species(221);
    ratelaws(841) = expressions(63)*species(224);
    ratelaws(842) = expressions(67)*species(224);
    ratelaws(843) = expressions(63)*species(219);
    ratelaws(844) = (expressions(64)/expressions(30))*species(219);
    ratelaws(845) = 2*expressions(63)*species(222);
    ratelaws(846) = 2*(expressions(66)/expressions(30))*species(222);
    ratelaws(847) = expressions(63)*species(223);
    ratelaws(848) = (expressions(65)/expressions(30))*species(223);
    ratelaws(849) = expressions(63)*species(225);
    ratelaws(850) = (expressions(67)/expressions(30))*species(225);

end

% Calculate species derivatives
function [ Dspecies ] = calc_species_deriv ( time, species, expressions )
    
    % initialize derivative vector
    Dspecies = zeros(225,1);
    
    % update observables
    [ observables ] = calc_observables( species, expressions );
    
    % update ratelaws
    [ ratelaws ] = calcratelaws( species, expressions, observables );
                        
    % calculate derivatives
    Dspecies(1) = -ratelaws(1) -ratelaws(2) -ratelaws(3) -ratelaws(4) -ratelaws(5) -ratelaws(6) +ratelaws(26) +ratelaws(27) +ratelaws(28) +ratelaws(29) +ratelaws(30) +ratelaws(31) +ratelaws(82) +ratelaws(83) +ratelaws(84) +ratelaws(85) +ratelaws(86) +ratelaws(87) +ratelaws(88) +ratelaws(89) +ratelaws(90) +ratelaws(91) +ratelaws(92) +ratelaws(93) -ratelaws(148) -ratelaws(149) +ratelaws(150) +ratelaws(151) +ratelaws(152) -ratelaws(153) -ratelaws(154) +ratelaws(155) +ratelaws(156) +ratelaws(157) -ratelaws(158) -ratelaws(159) +ratelaws(160) +ratelaws(161) +ratelaws(162) -ratelaws(163) -ratelaws(164) +ratelaws(165) -ratelaws(166) -ratelaws(167) +ratelaws(168) -ratelaws(169) -ratelaws(170) +ratelaws(171) -ratelaws(245) -ratelaws(246) -ratelaws(247) +ratelaws(248) -ratelaws(249) -ratelaws(250) -ratelaws(251) +ratelaws(252) -ratelaws(253) -ratelaws(254) -ratelaws(255) +ratelaws(256) -ratelaws(257) +ratelaws(258) -ratelaws(259) +ratelaws(260) +ratelaws(261) -ratelaws(262) -ratelaws(373) +ratelaws(374) -ratelaws(375) +ratelaws(376) -ratelaws(377) +ratelaws(378) -ratelaws(379) -ratelaws(380) -ratelaws(381) +ratelaws(382) -ratelaws(445) +ratelaws(446) -ratelaws(447) +ratelaws(448) -ratelaws(449) +ratelaws(450) -ratelaws(451) +ratelaws(452) +ratelaws(453) +ratelaws(454) +ratelaws(455) +ratelaws(456) +ratelaws(457) +ratelaws(458) +ratelaws(459) +ratelaws(460) -ratelaws(542) +ratelaws(543) -ratelaws(544) +ratelaws(545) -ratelaws(546) +ratelaws(547) -ratelaws(548) -ratelaws(549) -ratelaws(550) -ratelaws(551) -ratelaws(552) -ratelaws(553) -ratelaws(554) -ratelaws(555) -ratelaws(556) +ratelaws(557) +ratelaws(558) -ratelaws(633) -ratelaws(634) -ratelaws(635) -ratelaws(636) -ratelaws(637) +ratelaws(759) +ratelaws(760) +ratelaws(761) +ratelaws(762) +ratelaws(763) +ratelaws(764) +ratelaws(765) +ratelaws(766) -ratelaws(827) -ratelaws(828) -ratelaws(829) -ratelaws(830) -ratelaws(831) -ratelaws(832) -ratelaws(833) -ratelaws(834);
    Dspecies(2) = -ratelaws(1) +ratelaws(26) +ratelaws(198) +ratelaws(199) +ratelaws(349) +ratelaws(353);
    Dspecies(3) = -ratelaws(2) +ratelaws(27) +ratelaws(200) +ratelaws(201) +ratelaws(357) +ratelaws(361);
    Dspecies(4) = -ratelaws(3) +ratelaws(28) +ratelaws(202) +ratelaws(203) +ratelaws(365) +ratelaws(369);
    Dspecies(5) = ratelaws(8) -ratelaws(9) +ratelaws(10) -ratelaws(11) +ratelaws(12) -ratelaws(13) -ratelaws(14) -ratelaws(16) +ratelaws(22) -ratelaws(23) +ratelaws(36) -ratelaws(38) +ratelaws(40) -ratelaws(57) -ratelaws(58) -ratelaws(59) +ratelaws(68) +2.0*ratelaws(69) +ratelaws(71) +ratelaws(74) +ratelaws(94) +ratelaws(108) +ratelaws(109) +ratelaws(110) +ratelaws(111) +ratelaws(112) +ratelaws(113) +ratelaws(127) +ratelaws(172) +ratelaws(173) +ratelaws(174) +ratelaws(175) -ratelaws(176) +ratelaws(177) +ratelaws(178) +ratelaws(263) +ratelaws(264) -ratelaws(344) +ratelaws(346) +ratelaws(383) +ratelaws(384) -ratelaws(629) -ratelaws(647) -ratelaws(653) +ratelaws(767) +ratelaws(773) +ratelaws(803) +ratelaws(806) +ratelaws(807) +ratelaws(813) +ratelaws(818);
    Dspecies(6) = -ratelaws(8) -ratelaws(10) -ratelaws(12) -ratelaws(15) -ratelaws(17) +ratelaws(23) +ratelaws(37) -ratelaws(39) +ratelaws(42) +ratelaws(66) +ratelaws(67) +ratelaws(70) +ratelaws(71) +ratelaws(72) +ratelaws(95) +ratelaws(126) +ratelaws(127) +2.0*ratelaws(128) +ratelaws(344) +ratelaws(347) -ratelaws(630) -ratelaws(650) -ratelaws(656) +ratelaws(770) +ratelaws(776) +ratelaws(800) +ratelaws(809) +ratelaws(812) +ratelaws(815) +ratelaws(819) +ratelaws(821) +ratelaws(822);
    Dspecies(7) = ratelaws(9) +ratelaws(11) +ratelaws(13) -ratelaws(14) -ratelaws(15) -ratelaws(18) -ratelaws(21) -ratelaws(22) +ratelaws(33) +ratelaws(34) +ratelaws(35) +ratelaws(36) +ratelaws(37) +ratelaws(41) +ratelaws(44) -ratelaws(60) -ratelaws(61) -ratelaws(62) +ratelaws(67) +ratelaws(114) +ratelaws(115) +ratelaws(116) +ratelaws(117) +ratelaws(118) +ratelaws(119) +ratelaws(179) +ratelaws(180) +ratelaws(181) +ratelaws(182) -ratelaws(183) +ratelaws(184) +ratelaws(185) +ratelaws(265) +ratelaws(266) -ratelaws(345) -ratelaws(346) +ratelaws(385) +ratelaws(386) -ratelaws(631) -ratelaws(648) -ratelaws(651) +ratelaws(768) +ratelaws(771) +ratelaws(801);
    Dspecies(8) = -ratelaws(16) -ratelaws(17) -ratelaws(18) +ratelaws(40) +ratelaws(41) +ratelaws(42) +ratelaws(43) +ratelaws(44) +ratelaws(72) +ratelaws(73) +ratelaws(74) -ratelaws(659) -ratelaws(660) -ratelaws(661) +ratelaws(779) +ratelaws(780) +ratelaws(781) +ratelaws(782) +ratelaws(783) +ratelaws(824) +ratelaws(825) +ratelaws(826);
    Dspecies(9) = -ratelaws(4) +ratelaws(29) +ratelaws(172) +ratelaws(179) +ratelaws(297) +ratelaws(835) +ratelaws(843);
    Dspecies(10) = -ratelaws(5) +ratelaws(30) +ratelaws(174) +ratelaws(181) +ratelaws(298) +ratelaws(343) +ratelaws(593) +ratelaws(617) +ratelaws(837) +ratelaws(845);
    Dspecies(11) = -ratelaws(6) +ratelaws(31) +ratelaws(177) +ratelaws(184) +ratelaws(299) +ratelaws(841) +ratelaws(849);
    Dspecies(12) = -ratelaws(19) -ratelaws(20) +ratelaws(75) +ratelaws(76);
    Dspecies(13) = -ratelaws(186) -ratelaws(187) -ratelaws(188) -ratelaws(267) +ratelaws(268) -ratelaws(269) +ratelaws(270) -ratelaws(271) +ratelaws(272) +ratelaws(387) +ratelaws(388) +ratelaws(389) +ratelaws(492) +ratelaws(535) +ratelaws(754) +ratelaws(758);
    Dspecies(14) = -ratelaws(125) -ratelaws(194) +ratelaws(195) +ratelaws(282) +ratelaws(417) -ratelaws(479) -ratelaws(577) +ratelaws(579) +ratelaws(745);
    Dspecies(15) = -ratelaws(65) -ratelaws(123) +ratelaws(124) +ratelaws(193) -ratelaws(410) -ratelaws(481) +ratelaws(482) +ratelaws(581);
    Dspecies(16) = -ratelaws(197) -ratelaws(285) +ratelaws(286) +ratelaws(409);
    Dspecies(17) = -ratelaws(7) +ratelaws(32);
    Dspecies(18) = -ratelaws(25) -ratelaws(79) +ratelaws(81) -ratelaws(131) +ratelaws(133) +ratelaws(207) -ratelaws(427) -ratelaws(498) +ratelaws(505) +ratelaws(594);
    Dspecies(19) = -ratelaws(25) +ratelaws(81) -ratelaws(135) +ratelaws(136) -ratelaws(212) +ratelaws(213) -ratelaws(512) +ratelaws(513) -ratelaws(601) +ratelaws(602);
    Dspecies(20) = -ratelaws(140) +ratelaws(223) +ratelaws(224);
    Dspecies(21) = -ratelaws(78) -ratelaws(137) -ratelaws(139) -ratelaws(214) +ratelaws(215) -ratelaws(218) +ratelaws(219) -ratelaws(221) -ratelaws(222) +ratelaws(224) +ratelaws(312) +ratelaws(314) -ratelaws(316) +ratelaws(317) -ratelaws(319) +ratelaws(320) +ratelaws(434) +ratelaws(436) -ratelaws(514) +ratelaws(515) -ratelaws(518) +ratelaws(519) -ratelaws(521) +ratelaws(522) -ratelaws(524) +ratelaws(525) -ratelaws(603) +ratelaws(604) -ratelaws(607) +ratelaws(608) -ratelaws(610) +ratelaws(611) -ratelaws(613) +ratelaws(614);
    Dspecies(22) = -ratelaws(141) +ratelaws(225);
    Dspecies(23) = -ratelaws(226) +ratelaws(323);
    Dspecies(24) = -ratelaws(439) +ratelaws(528);
    Dspecies(25) = -ratelaws(300) -ratelaws(418) +ratelaws(419) +ratelaws(496);
    Dspecies(26) = -ratelaws(24) +ratelaws(77);
    Dspecies(27) = -ratelaws(420) -ratelaws(497);
    Dspecies(28) = -ratelaws(24) +ratelaws(77) +ratelaws(420) +ratelaws(497);
    Dspecies(29) = -ratelaws(142) -ratelaws(143) -ratelaws(144) -ratelaws(145) -ratelaws(146) -ratelaws(147) +ratelaws(227) +ratelaws(228) +ratelaws(229) +ratelaws(230) +ratelaws(231) +ratelaws(232) +ratelaws(324) +ratelaws(325) +ratelaws(326) +ratelaws(327) +ratelaws(328) +ratelaws(329) +ratelaws(372);
    Dspecies(30) = ratelaws(1) -ratelaws(26) -ratelaws(45) -ratelaws(46) -ratelaws(51) -ratelaws(52) +ratelaws(348) +ratelaws(352) -ratelaws(662) -ratelaws(664) -ratelaws(674) -ratelaws(676);
    Dspecies(31) = ratelaws(2) -ratelaws(27) -ratelaws(47) -ratelaws(48) -ratelaws(53) -ratelaws(54) +ratelaws(356) +ratelaws(360) -ratelaws(666) -ratelaws(668) -ratelaws(678) -ratelaws(680);
    Dspecies(32) = ratelaws(3) -ratelaws(28) -ratelaws(49) -ratelaws(50) -ratelaws(55) -ratelaws(56) +ratelaws(364) +ratelaws(368) -ratelaws(670) -ratelaws(672) -ratelaws(682) -ratelaws(684);
    Dspecies(33) = ratelaws(4) -ratelaws(29) -ratelaws(57) -ratelaws(60) +ratelaws(108) +ratelaws(114) -ratelaws(686) -ratelaws(690) +ratelaws(784) +ratelaws(792);
    Dspecies(34) = ratelaws(5) -ratelaws(30) -ratelaws(58) -ratelaws(61) +ratelaws(110) +ratelaws(116) +ratelaws(342) +ratelaws(530) -ratelaws(687) -ratelaws(691) +ratelaws(786) +ratelaws(794);
    Dspecies(35) = ratelaws(6) -ratelaws(31) -ratelaws(59) -ratelaws(62) +ratelaws(112) +ratelaws(118) -ratelaws(689) -ratelaws(693) +ratelaws(790) +ratelaws(798);
    Dspecies(36) = ratelaws(7) -ratelaws(32);
    Dspecies(37) = ratelaws(14) -ratelaws(36) -ratelaws(67) -ratelaws(68) -ratelaws(69);
    Dspecies(38) = ratelaws(15) -ratelaws(37) -ratelaws(70) -ratelaws(71);
    Dspecies(39) = ratelaws(16) -ratelaws(40) -ratelaws(41) -ratelaws(72);
    Dspecies(40) = ratelaws(17) -ratelaws(42) -ratelaws(43);
    Dspecies(41) = ratelaws(18) -ratelaws(44) -ratelaws(73) -ratelaws(74);
    Dspecies(42) = ratelaws(19) -ratelaws(64) -ratelaws(65) -ratelaws(75) -ratelaws(79) -ratelaws(80) +ratelaws(124) +ratelaws(130) +ratelaws(133) +ratelaws(134) +ratelaws(209) +ratelaws(210) +ratelaws(211) -ratelaws(301) -ratelaws(302) -ratelaws(303) +ratelaws(310) +ratelaws(311) -ratelaws(421) -ratelaws(422);
    Dspecies(43) = ratelaws(20) -ratelaws(63) -ratelaws(76) +ratelaws(129);
    Dspecies(44) = ratelaws(21) -ratelaws(33) -ratelaws(34) -ratelaws(35) -ratelaws(38) -ratelaws(39) +ratelaws(43) -ratelaws(66) +ratelaws(68) +ratelaws(70) +ratelaws(73) +ratelaws(94) +ratelaws(95) +ratelaws(126) +ratelaws(345) -ratelaws(347) -ratelaws(632) -ratelaws(654) -ratelaws(657) +ratelaws(774) +ratelaws(777) +ratelaws(804) +ratelaws(810) +ratelaws(816);
    Dspecies(45) = ratelaws(24) -ratelaws(77);
    Dspecies(46) = ratelaws(25) -ratelaws(80) -ratelaws(81) -ratelaws(132) +ratelaws(134) +ratelaws(208) -ratelaws(428) -ratelaws(499) +ratelaws(506) +ratelaws(595);
    Dspecies(47) = ratelaws(38) -ratelaws(94) -ratelaws(126) -ratelaws(127);
    Dspecies(48) = ratelaws(39) -ratelaws(95) -ratelaws(128);
    Dspecies(49) = ratelaws(45) +ratelaws(51) -ratelaws(82) -ratelaws(97) -ratelaws(103) -ratelaws(142) +ratelaws(148) +ratelaws(227) +ratelaws(233) +ratelaws(239) -ratelaws(348) +ratelaws(354) +ratelaws(662) -ratelaws(665) +ratelaws(674) -ratelaws(677);
    Dspecies(50) = ratelaws(46) +ratelaws(52) -ratelaws(83) -ratelaws(96) -ratelaws(102) -ratelaws(143) +ratelaws(149) +ratelaws(228) +ratelaws(234) +ratelaws(240) +ratelaws(350) -ratelaws(352) -ratelaws(663) +ratelaws(664) -ratelaws(675) +ratelaws(676);
    Dspecies(51) = ratelaws(47) +ratelaws(53) -ratelaws(84) -ratelaws(99) -ratelaws(105) -ratelaws(144) +ratelaws(153) +ratelaws(229) +ratelaws(235) +ratelaws(241) -ratelaws(356) +ratelaws(362) +ratelaws(666) -ratelaws(669) +ratelaws(678) -ratelaws(681);
    Dspecies(52) = ratelaws(48) +ratelaws(54) -ratelaws(85) -ratelaws(98) -ratelaws(104) -ratelaws(145) +ratelaws(154) +ratelaws(230) +ratelaws(236) +ratelaws(242) +ratelaws(358) -ratelaws(360) -ratelaws(667) +ratelaws(668) -ratelaws(679) +ratelaws(680);
    Dspecies(53) = ratelaws(49) +ratelaws(55) -ratelaws(86) -ratelaws(101) -ratelaws(107) -ratelaws(146) +ratelaws(158) +ratelaws(231) +ratelaws(237) +ratelaws(243) -ratelaws(364) +ratelaws(370) +ratelaws(670) -ratelaws(673) +ratelaws(682) -ratelaws(685);
    Dspecies(54) = ratelaws(50) +ratelaws(56) -ratelaws(87) -ratelaws(100) -ratelaws(106) -ratelaws(147) +ratelaws(159) +ratelaws(232) +ratelaws(238) +ratelaws(244) +ratelaws(366) -ratelaws(368) -ratelaws(671) +ratelaws(672) -ratelaws(683) +ratelaws(684);
    Dspecies(55) = ratelaws(57) -ratelaws(88) -ratelaws(108) -ratelaws(109) +ratelaws(163);
    Dspecies(56) = ratelaws(58) -ratelaws(90) -ratelaws(110) -ratelaws(111) +ratelaws(166);
    Dspecies(57) = ratelaws(59) -ratelaws(92) -ratelaws(112) -ratelaws(113) +ratelaws(169);
    Dspecies(58) = ratelaws(60) -ratelaws(89) -ratelaws(114) -ratelaws(115) +ratelaws(164);
    Dspecies(59) = ratelaws(61) -ratelaws(91) -ratelaws(116) -ratelaws(117) +ratelaws(167);
    Dspecies(60) = ratelaws(62) -ratelaws(93) -ratelaws(118) -ratelaws(119) +ratelaws(170);
    Dspecies(61) = ratelaws(63) +ratelaws(64) -ratelaws(123) -ratelaws(129) -ratelaws(130) -ratelaws(131) -ratelaws(132) +ratelaws(193) +ratelaws(207) +ratelaws(208) -ratelaws(304) -ratelaws(305) -ratelaws(306) +ratelaws(307) +ratelaws(308) +ratelaws(309) -ratelaws(423) -ratelaws(424) +ratelaws(425) +ratelaws(426);
    Dspecies(62) = ratelaws(65) -ratelaws(120) -ratelaws(124) -ratelaws(125) +ratelaws(195) +ratelaws(204) +ratelaws(283) -ratelaws(406);
    Dspecies(63) = ratelaws(78) -ratelaws(138) -ratelaws(140) -ratelaws(216) +ratelaws(217) +ratelaws(220) +ratelaws(223) +ratelaws(313) +ratelaws(315) +ratelaws(318) +ratelaws(321) +ratelaws(435) +ratelaws(437) -ratelaws(516) +ratelaws(517) +ratelaws(520) +ratelaws(523) +ratelaws(526) -ratelaws(605) +ratelaws(606) +ratelaws(609) +ratelaws(612) +ratelaws(615);
    Dspecies(64) = ratelaws(79) -ratelaws(121) -ratelaws(133) -ratelaws(135) +ratelaws(136) +ratelaws(205);
    Dspecies(65) = ratelaws(80) -ratelaws(122) -ratelaws(134) +ratelaws(135) -ratelaws(136) -ratelaws(137) -ratelaws(138) -ratelaws(139) +ratelaws(206) +ratelaws(215) +ratelaws(217) +ratelaws(219) +ratelaws(220);
    Dspecies(66) = ratelaws(82) -ratelaws(148) -ratelaws(198) +ratelaws(288) +ratelaws(324) +ratelaws(330) +ratelaws(336) -ratelaws(349) +ratelaws(355);
    Dspecies(67) = ratelaws(83) -ratelaws(149) -ratelaws(199) +ratelaws(287) +ratelaws(325) +ratelaws(331) +ratelaws(337) +ratelaws(351) -ratelaws(353);
    Dspecies(68) = ratelaws(84) -ratelaws(153) -ratelaws(200) +ratelaws(290) +ratelaws(326) +ratelaws(332) +ratelaws(338) -ratelaws(357) +ratelaws(363);
    Dspecies(69) = ratelaws(85) -ratelaws(154) -ratelaws(201) +ratelaws(289) +ratelaws(327) +ratelaws(333) +ratelaws(339) +ratelaws(359) -ratelaws(361);
    Dspecies(70) = ratelaws(86) -ratelaws(158) -ratelaws(202) +ratelaws(292) +ratelaws(328) +ratelaws(334) +ratelaws(340) -ratelaws(365) +ratelaws(371);
    Dspecies(71) = ratelaws(87) -ratelaws(159) -ratelaws(203) +ratelaws(291) +ratelaws(329) +ratelaws(335) +ratelaws(341) +ratelaws(367) -ratelaws(369);
    Dspecies(72) = ratelaws(88) -ratelaws(163) -ratelaws(172) -ratelaws(173);
    Dspecies(73) = ratelaws(89) -ratelaws(164) -ratelaws(179) -ratelaws(180);
    Dspecies(74) = ratelaws(90) -ratelaws(166) -ratelaws(174) -ratelaws(175);
    Dspecies(75) = ratelaws(91) -ratelaws(167) -ratelaws(181) -ratelaws(182);
    Dspecies(76) = ratelaws(92) -ratelaws(169) -ratelaws(177) -ratelaws(178);
    Dspecies(77) = ratelaws(93) -ratelaws(170) -ratelaws(184) -ratelaws(185);
    Dspecies(78) = ratelaws(96) +ratelaws(97) +ratelaws(102) +ratelaws(103) -ratelaws(150) -ratelaws(186) +ratelaws(245) +ratelaws(268) -ratelaws(350) -ratelaws(354) +ratelaws(390) -ratelaws(461) +ratelaws(470) +ratelaws(559) -ratelaws(562) +ratelaws(663) +ratelaws(665) +ratelaws(675) +ratelaws(677) -ratelaws(694);
    Dspecies(79) = ratelaws(98) +ratelaws(99) +ratelaws(104) +ratelaws(105) -ratelaws(155) -ratelaws(187) +ratelaws(249) +ratelaws(270) -ratelaws(358) -ratelaws(362) +ratelaws(391) -ratelaws(464) +ratelaws(471) +ratelaws(560) -ratelaws(565) +ratelaws(667) +ratelaws(669) +ratelaws(679) +ratelaws(681) -ratelaws(697);
    Dspecies(80) = ratelaws(100) +ratelaws(101) +ratelaws(106) +ratelaws(107) -ratelaws(160) -ratelaws(188) +ratelaws(253) +ratelaws(272) -ratelaws(366) -ratelaws(370) +ratelaws(392) -ratelaws(467) +ratelaws(472) +ratelaws(561) -ratelaws(568) +ratelaws(671) +ratelaws(673) +ratelaws(683) +ratelaws(685) -ratelaws(700);
    Dspecies(81) = ratelaws(109) +ratelaws(115) -ratelaws(165) -ratelaws(197) +ratelaws(257) +ratelaws(286) +ratelaws(785) +ratelaws(793);
    Dspecies(82) = ratelaws(111) +ratelaws(117) -ratelaws(168) -ratelaws(176) -ratelaws(183) +ratelaws(259) +ratelaws(263) +ratelaws(265) -ratelaws(342) +ratelaws(440) -ratelaws(688) -ratelaws(692) +ratelaws(787) +ratelaws(788) +ratelaws(795) +ratelaws(796);
    Dspecies(83) = ratelaws(113) +ratelaws(119) -ratelaws(171) +ratelaws(262) +ratelaws(791) +ratelaws(799);
    Dspecies(84) = ratelaws(120) +ratelaws(123) -ratelaws(193) -ratelaws(194) -ratelaws(204) +ratelaws(282) -ratelaws(407) +ratelaws(408);
    Dspecies(85) = ratelaws(121) +ratelaws(131) -ratelaws(205) -ratelaws(207) -ratelaws(212) +ratelaws(213);
    Dspecies(86) = ratelaws(122) +ratelaws(132) -ratelaws(206) -ratelaws(208) +ratelaws(212) -ratelaws(213) -ratelaws(214) -ratelaws(216) -ratelaws(218) +ratelaws(312) +ratelaws(313) +ratelaws(314) +ratelaws(315);
    Dspecies(87) = ratelaws(125) -ratelaws(189) -ratelaws(195) -ratelaws(196) +ratelaws(293);
    Dspecies(88) = ratelaws(137) -ratelaws(190) -ratelaws(209) -ratelaws(215) -ratelaws(221) +ratelaws(294) +ratelaws(301) +ratelaws(317) +ratelaws(318);
    Dspecies(89) = ratelaws(138) -ratelaws(191) -ratelaws(210) -ratelaws(217) -ratelaws(222) +ratelaws(295) +ratelaws(302) +ratelaws(320) +ratelaws(321);
    Dspecies(90) = ratelaws(139) -ratelaws(192) -ratelaws(211) -ratelaws(219) -ratelaws(220) +ratelaws(296) +ratelaws(303);
    Dspecies(91) = ratelaws(140) -ratelaws(223) -ratelaws(224);
    Dspecies(92) = ratelaws(141) -ratelaws(225);
    Dspecies(93) = ratelaws(142) -ratelaws(151) -ratelaws(227) -ratelaws(233) -ratelaws(239) +ratelaws(246);
    Dspecies(94) = ratelaws(143) -ratelaws(152) -ratelaws(228) -ratelaws(234) -ratelaws(240) +ratelaws(247);
    Dspecies(95) = ratelaws(144) -ratelaws(156) -ratelaws(229) -ratelaws(235) -ratelaws(241) +ratelaws(250);
    Dspecies(96) = ratelaws(145) -ratelaws(157) -ratelaws(230) -ratelaws(236) -ratelaws(242) +ratelaws(251);
    Dspecies(97) = ratelaws(146) -ratelaws(161) -ratelaws(231) -ratelaws(237) -ratelaws(243) +ratelaws(254);
    Dspecies(98) = ratelaws(147) -ratelaws(162) -ratelaws(232) -ratelaws(238) -ratelaws(244) +ratelaws(255);
    Dspecies(99) = ratelaws(150) -ratelaws(245) -ratelaws(267) -ratelaws(287) -ratelaws(288) -ratelaws(351) -ratelaws(355) +ratelaws(387) -ratelaws(462) +ratelaws(463) -ratelaws(563) +ratelaws(564) -ratelaws(695) +ratelaws(696);
    Dspecies(100) = ratelaws(151) -ratelaws(246) -ratelaws(324) -ratelaws(330) -ratelaws(336);
    Dspecies(101) = ratelaws(152) -ratelaws(247) -ratelaws(325) -ratelaws(331) -ratelaws(337);
    Dspecies(102) = ratelaws(155) -ratelaws(249) -ratelaws(269) -ratelaws(289) -ratelaws(290) -ratelaws(359) -ratelaws(363) +ratelaws(388) -ratelaws(465) +ratelaws(466) -ratelaws(566) +ratelaws(567) -ratelaws(698) +ratelaws(699);
    Dspecies(103) = ratelaws(156) -ratelaws(250) -ratelaws(326) -ratelaws(332) -ratelaws(338);
    Dspecies(104) = ratelaws(157) -ratelaws(251) -ratelaws(327) -ratelaws(333) -ratelaws(339);
    Dspecies(105) = ratelaws(160) -ratelaws(253) -ratelaws(271) -ratelaws(291) -ratelaws(292) -ratelaws(367) -ratelaws(371) +ratelaws(389) -ratelaws(468) +ratelaws(469) -ratelaws(569) +ratelaws(570) -ratelaws(701) +ratelaws(702);
    Dspecies(106) = ratelaws(161) -ratelaws(254) -ratelaws(328) -ratelaws(334) -ratelaws(340);
    Dspecies(107) = ratelaws(162) -ratelaws(255) -ratelaws(329) -ratelaws(335) -ratelaws(341);
    Dspecies(108) = ratelaws(165) +ratelaws(173) +ratelaws(180) -ratelaws(257) -ratelaws(285) -ratelaws(297) +ratelaws(409) +ratelaws(836) +ratelaws(844);
    Dspecies(109) = ratelaws(168) +ratelaws(175) +ratelaws(182) -ratelaws(259) -ratelaws(298) -ratelaws(343) +ratelaws(383) +ratelaws(385) +ratelaws(494) +ratelaws(529) +ratelaws(838) +ratelaws(839) +ratelaws(846) +ratelaws(847);
    Dspecies(110) = ratelaws(171) +ratelaws(178) +ratelaws(185) -ratelaws(262) -ratelaws(299) +ratelaws(842) +ratelaws(850);
    Dspecies(111) = ratelaws(176) -ratelaws(260) -ratelaws(263) -ratelaws(264) +ratelaws(380);
    Dspecies(112) = ratelaws(183) -ratelaws(261) -ratelaws(265) -ratelaws(266) +ratelaws(381);
    Dspecies(113) = ratelaws(186) -ratelaws(248) -ratelaws(268) -ratelaws(273) -ratelaws(276) +ratelaws(373) +ratelaws(411) +ratelaws(442) +ratelaws(588) +ratelaws(626) -ratelaws(703) -ratelaws(718);
    Dspecies(114) = ratelaws(187) -ratelaws(252) -ratelaws(270) -ratelaws(274) -ratelaws(277) +ratelaws(375) +ratelaws(412) +ratelaws(443) +ratelaws(589) +ratelaws(627) -ratelaws(704) -ratelaws(719);
    Dspecies(115) = ratelaws(188) -ratelaws(256) -ratelaws(272) -ratelaws(275) -ratelaws(278) +ratelaws(377) +ratelaws(413) +ratelaws(444) +ratelaws(590) +ratelaws(628) -ratelaws(705) -ratelaws(720);
    Dspecies(116) = ratelaws(189) +ratelaws(194) -ratelaws(282) -ratelaws(284) -ratelaws(293);
    Dspecies(117) = ratelaws(190) +ratelaws(214) -ratelaws(294) +ratelaws(304) -ratelaws(307) -ratelaws(312) -ratelaws(316) +ratelaws(434) +ratelaws(435);
    Dspecies(118) = ratelaws(191) +ratelaws(216) -ratelaws(295) +ratelaws(305) -ratelaws(308) -ratelaws(313) -ratelaws(319) +ratelaws(436) +ratelaws(437);
    Dspecies(119) = ratelaws(192) +ratelaws(218) -ratelaws(296) +ratelaws(306) -ratelaws(309) -ratelaws(314) -ratelaws(315);
    Dspecies(120) = ratelaws(196) -ratelaws(279) -ratelaws(283) -ratelaws(300) +ratelaws(406) +ratelaws(414) +ratelaws(419);
    Dspecies(121) = ratelaws(197) -ratelaws(258) -ratelaws(286) +ratelaws(379);
    Dspecies(122) = ratelaws(209) -ratelaws(301) -ratelaws(304) +ratelaws(307) -ratelaws(429) -ratelaws(500) +ratelaws(507) +ratelaws(596);
    Dspecies(123) = ratelaws(210) -ratelaws(302) -ratelaws(305) +ratelaws(308) -ratelaws(430) -ratelaws(501) +ratelaws(508) +ratelaws(597);
    Dspecies(124) = ratelaws(211) -ratelaws(303) -ratelaws(306) +ratelaws(309) -ratelaws(431) -ratelaws(502) +ratelaws(509) +ratelaws(598);
    Dspecies(125) = ratelaws(221) -ratelaws(280) -ratelaws(310) -ratelaws(317) -ratelaws(318) +ratelaws(415) +ratelaws(421);
    Dspecies(126) = ratelaws(222) -ratelaws(281) -ratelaws(311) -ratelaws(320) -ratelaws(321) +ratelaws(416) +ratelaws(422);
    Dspecies(127) = ratelaws(226) -ratelaws(322) -ratelaws(323) +ratelaws(438);
    Dspecies(128) = ratelaws(233) +ratelaws(234) +ratelaws(235) +ratelaws(236) +ratelaws(237) +ratelaws(238) +ratelaws(239) +ratelaws(240) +ratelaws(241) +ratelaws(242) +ratelaws(243) +ratelaws(244) +ratelaws(330) +ratelaws(331) +ratelaws(332) +ratelaws(333) +ratelaws(334) +ratelaws(335) +ratelaws(336) +ratelaws(337) +ratelaws(338) +ratelaws(339) +ratelaws(340) +ratelaws(341) -ratelaws(372);
    Dspecies(129) = ratelaws(248) +ratelaws(267) -ratelaws(373) -ratelaws(387) -ratelaws(393) -ratelaws(396) +ratelaws(483) +ratelaws(532) -ratelaws(706) -ratelaws(721) +ratelaws(751) +ratelaws(755);
    Dspecies(130) = ratelaws(252) +ratelaws(269) -ratelaws(375) -ratelaws(388) -ratelaws(394) -ratelaws(397) +ratelaws(484) +ratelaws(533) -ratelaws(707) -ratelaws(722) +ratelaws(752) +ratelaws(756);
    Dspecies(131) = ratelaws(256) +ratelaws(271) -ratelaws(377) -ratelaws(389) -ratelaws(395) -ratelaws(398) +ratelaws(485) +ratelaws(534) -ratelaws(708) -ratelaws(723) +ratelaws(753) +ratelaws(757);
    Dspecies(132) = ratelaws(258) +ratelaws(285) -ratelaws(379) -ratelaws(409);
    Dspecies(133) = ratelaws(260) -ratelaws(380) -ratelaws(383) -ratelaws(384);
    Dspecies(134) = ratelaws(261) -ratelaws(381) -ratelaws(385) -ratelaws(386);
    Dspecies(135) = ratelaws(264) +ratelaws(266) -ratelaws(382) -ratelaws(410) -ratelaws(427) -ratelaws(428) -ratelaws(429) -ratelaws(430) -ratelaws(431) -ratelaws(432) -ratelaws(433) -ratelaws(440) -ratelaws(441) +ratelaws(451) +ratelaws(482) +ratelaws(505) +ratelaws(506) +ratelaws(507) +ratelaws(508) +ratelaws(509) +ratelaws(510) +ratelaws(511) +ratelaws(789) +ratelaws(797);
    Dspecies(136) = ratelaws(273) +ratelaws(276) -ratelaws(374) -ratelaws(390) -ratelaws(399) -ratelaws(402) -ratelaws(411) -ratelaws(442) +ratelaws(445) +ratelaws(461) +ratelaws(489) +ratelaws(539) +ratelaws(703) +ratelaws(718) -ratelaws(733) -ratelaws(739);
    Dspecies(137) = ratelaws(274) +ratelaws(277) -ratelaws(376) -ratelaws(391) -ratelaws(400) -ratelaws(403) -ratelaws(412) -ratelaws(443) +ratelaws(447) +ratelaws(464) +ratelaws(490) +ratelaws(540) +ratelaws(704) +ratelaws(719) -ratelaws(734) -ratelaws(740);
    Dspecies(138) = ratelaws(275) +ratelaws(278) -ratelaws(378) -ratelaws(392) -ratelaws(401) -ratelaws(404) -ratelaws(413) -ratelaws(444) +ratelaws(449) +ratelaws(467) +ratelaws(491) +ratelaws(541) +ratelaws(705) +ratelaws(720) -ratelaws(735) -ratelaws(741);
    Dspecies(139) = ratelaws(279) +ratelaws(284) +ratelaws(407) -ratelaws(408) -ratelaws(414) -ratelaws(418) +ratelaws(496);
    Dspecies(140) = ratelaws(280) +ratelaws(316) -ratelaws(415) +ratelaws(423) -ratelaws(425) -ratelaws(434) -ratelaws(435);
    Dspecies(141) = ratelaws(281) +ratelaws(319) -ratelaws(416) +ratelaws(424) -ratelaws(426) -ratelaws(436) -ratelaws(437);
    Dspecies(142) = ratelaws(283) -ratelaws(406) -ratelaws(407) +ratelaws(408) -ratelaws(417) -ratelaws(480) -ratelaws(578) +ratelaws(580) +ratelaws(746);
    Dspecies(143) = ratelaws(300) -ratelaws(405) -ratelaws(419) +ratelaws(493);
    Dspecies(144) = ratelaws(310) -ratelaws(421) -ratelaws(423) +ratelaws(425) -ratelaws(432) -ratelaws(503) +ratelaws(510) +ratelaws(599);
    Dspecies(145) = ratelaws(311) -ratelaws(422) -ratelaws(424) +ratelaws(426) -ratelaws(433) -ratelaws(504) +ratelaws(511) +ratelaws(600);
    Dspecies(146) = ratelaws(322) -ratelaws(438);
    Dspecies(147) = ratelaws(374) +ratelaws(393) +ratelaws(396) -ratelaws(445) +ratelaws(462) -ratelaws(463) -ratelaws(473) -ratelaws(476) -ratelaws(483) -ratelaws(532) +ratelaws(585) +ratelaws(622) +ratelaws(706) +ratelaws(721) -ratelaws(736) -ratelaws(742);
    Dspecies(148) = ratelaws(376) +ratelaws(394) +ratelaws(397) -ratelaws(447) +ratelaws(465) -ratelaws(466) -ratelaws(474) -ratelaws(477) -ratelaws(484) -ratelaws(533) +ratelaws(586) +ratelaws(623) +ratelaws(707) +ratelaws(722) -ratelaws(737) -ratelaws(743);
    Dspecies(149) = ratelaws(378) +ratelaws(395) +ratelaws(398) -ratelaws(449) +ratelaws(468) -ratelaws(469) -ratelaws(475) -ratelaws(478) -ratelaws(485) -ratelaws(534) +ratelaws(587) +ratelaws(624) +ratelaws(708) +ratelaws(723) -ratelaws(738) -ratelaws(744);
    Dspecies(150) = ratelaws(382) +ratelaws(384) +ratelaws(386) -ratelaws(451) -ratelaws(481) -ratelaws(494) -ratelaws(495) -ratelaws(498) -ratelaws(499) -ratelaws(500) -ratelaws(501) -ratelaws(502) -ratelaws(503) -ratelaws(504) -ratelaws(529) -ratelaws(531) +ratelaws(581) +ratelaws(594) +ratelaws(595) +ratelaws(596) +ratelaws(597) +ratelaws(598) +ratelaws(599) +ratelaws(600) +ratelaws(840) +ratelaws(848);
    Dspecies(151) = ratelaws(390) +ratelaws(391) +ratelaws(392) -ratelaws(461) -ratelaws(462) +ratelaws(463) -ratelaws(464) -ratelaws(465) +ratelaws(466) -ratelaws(467) -ratelaws(468) +ratelaws(469) -ratelaws(492) -ratelaws(535) +ratelaws(592) +ratelaws(625);
    Dspecies(152) = ratelaws(399) +ratelaws(402) -ratelaws(446) -ratelaws(470) -ratelaws(486) -ratelaws(489) -ratelaws(536) -ratelaws(539) +ratelaws(542) +ratelaws(562) +ratelaws(571) +ratelaws(574) +ratelaws(709) +ratelaws(724) +ratelaws(733) +ratelaws(739);
    Dspecies(153) = ratelaws(400) +ratelaws(403) -ratelaws(448) -ratelaws(471) -ratelaws(487) -ratelaws(490) -ratelaws(537) -ratelaws(540) +ratelaws(544) +ratelaws(565) +ratelaws(572) +ratelaws(575) +ratelaws(710) +ratelaws(725) +ratelaws(734) +ratelaws(740);
    Dspecies(154) = ratelaws(401) +ratelaws(404) -ratelaws(450) -ratelaws(472) -ratelaws(488) -ratelaws(491) -ratelaws(538) -ratelaws(541) +ratelaws(546) +ratelaws(568) +ratelaws(573) +ratelaws(576) +ratelaws(711) +ratelaws(726) +ratelaws(735) +ratelaws(741);
    Dspecies(155) = ratelaws(405) +ratelaws(418) -ratelaws(493) -ratelaws(496);
    Dspecies(156) = ratelaws(410) -ratelaws(452) -ratelaws(479) -ratelaws(480) -ratelaws(482) +ratelaws(548) +ratelaws(579) +ratelaws(580);
    Dspecies(157) = ratelaws(427) -ratelaws(453) -ratelaws(505) -ratelaws(512) +ratelaws(513) +ratelaws(549);
    Dspecies(158) = ratelaws(428) -ratelaws(454) -ratelaws(506) +ratelaws(512) -ratelaws(513) -ratelaws(514) +ratelaws(515) -ratelaws(516) +ratelaws(517) -ratelaws(518) +ratelaws(519) +ratelaws(520) +ratelaws(550);
    Dspecies(159) = ratelaws(429) -ratelaws(455) -ratelaws(507) +ratelaws(514) -ratelaws(515) -ratelaws(521) +ratelaws(522) +ratelaws(523) +ratelaws(551);
    Dspecies(160) = ratelaws(430) -ratelaws(456) -ratelaws(508) +ratelaws(516) -ratelaws(517) -ratelaws(524) +ratelaws(525) +ratelaws(526) +ratelaws(552);
    Dspecies(161) = ratelaws(431) -ratelaws(457) -ratelaws(509) +ratelaws(518) -ratelaws(519) -ratelaws(520) +ratelaws(553);
    Dspecies(162) = ratelaws(432) -ratelaws(458) -ratelaws(510) +ratelaws(521) -ratelaws(522) -ratelaws(523) +ratelaws(554);
    Dspecies(163) = ratelaws(433) -ratelaws(459) -ratelaws(511) +ratelaws(524) -ratelaws(525) -ratelaws(526) +ratelaws(555);
    Dspecies(164) = ratelaws(439) -ratelaws(527) -ratelaws(528) +ratelaws(616);
    Dspecies(165) = ratelaws(441) -ratelaws(460) -ratelaws(530) +ratelaws(556);
    Dspecies(166) = ratelaws(446) +ratelaws(473) +ratelaws(476) -ratelaws(542) +ratelaws(563) -ratelaws(564) -ratelaws(582) -ratelaws(585) -ratelaws(618) -ratelaws(622) +ratelaws(712) +ratelaws(713) +ratelaws(727) +ratelaws(728) +ratelaws(736) +ratelaws(742);
    Dspecies(167) = ratelaws(448) +ratelaws(474) +ratelaws(477) -ratelaws(544) +ratelaws(566) -ratelaws(567) -ratelaws(583) -ratelaws(586) -ratelaws(619) -ratelaws(623) +ratelaws(714) +ratelaws(715) +ratelaws(729) +ratelaws(730) +ratelaws(737) +ratelaws(743);
    Dspecies(168) = ratelaws(450) +ratelaws(475) +ratelaws(478) -ratelaws(546) +ratelaws(569) -ratelaws(570) -ratelaws(584) -ratelaws(587) -ratelaws(620) -ratelaws(624) +ratelaws(716) +ratelaws(717) +ratelaws(731) +ratelaws(732) +ratelaws(738) +ratelaws(744);
    Dspecies(169) = ratelaws(452) +ratelaws(481) -ratelaws(548) -ratelaws(577) -ratelaws(578) -ratelaws(581) +ratelaws(745) +ratelaws(746);
    Dspecies(170) = ratelaws(453) +ratelaws(498) -ratelaws(549) -ratelaws(594) -ratelaws(601) +ratelaws(602);
    Dspecies(171) = ratelaws(454) +ratelaws(499) -ratelaws(550) -ratelaws(595) +ratelaws(601) -ratelaws(602) -ratelaws(603) +ratelaws(604) -ratelaws(605) +ratelaws(606) -ratelaws(607) +ratelaws(608) +ratelaws(609);
    Dspecies(172) = ratelaws(455) +ratelaws(500) -ratelaws(551) -ratelaws(596) +ratelaws(603) -ratelaws(604) -ratelaws(610) +ratelaws(611) +ratelaws(612);
    Dspecies(173) = ratelaws(456) +ratelaws(501) -ratelaws(552) -ratelaws(597) +ratelaws(605) -ratelaws(606) -ratelaws(613) +ratelaws(614) +ratelaws(615);
    Dspecies(174) = ratelaws(457) +ratelaws(502) -ratelaws(553) -ratelaws(598) +ratelaws(607) -ratelaws(608) -ratelaws(609);
    Dspecies(175) = ratelaws(458) +ratelaws(503) -ratelaws(554) -ratelaws(599) +ratelaws(610) -ratelaws(611) -ratelaws(612);
    Dspecies(176) = ratelaws(459) +ratelaws(504) -ratelaws(555) -ratelaws(600) +ratelaws(613) -ratelaws(614) -ratelaws(615);
    Dspecies(177) = ratelaws(460) +ratelaws(495) +ratelaws(531) -ratelaws(556) -ratelaws(593) -ratelaws(617);
    Dspecies(178) = ratelaws(470) +ratelaws(471) +ratelaws(472) -ratelaws(562) -ratelaws(563) +ratelaws(564) -ratelaws(565) -ratelaws(566) +ratelaws(567) -ratelaws(568) -ratelaws(569) +ratelaws(570) -ratelaws(591) -ratelaws(592) -ratelaws(621) -ratelaws(625);
    Dspecies(179) = ratelaws(479) -ratelaws(557) -ratelaws(579) +ratelaws(636);
    Dspecies(180) = ratelaws(480) -ratelaws(558) -ratelaws(580) +ratelaws(637);
    Dspecies(181) = ratelaws(486) +ratelaws(536) -ratelaws(543) -ratelaws(559) -ratelaws(571) -ratelaws(574) -ratelaws(588) -ratelaws(626) +ratelaws(633) +ratelaws(694) -ratelaws(709) -ratelaws(724);
    Dspecies(182) = ratelaws(487) +ratelaws(537) -ratelaws(545) -ratelaws(560) -ratelaws(572) -ratelaws(575) -ratelaws(589) -ratelaws(627) +ratelaws(634) +ratelaws(697) -ratelaws(710) -ratelaws(725);
    Dspecies(183) = ratelaws(488) +ratelaws(538) -ratelaws(547) -ratelaws(561) -ratelaws(573) -ratelaws(576) -ratelaws(590) -ratelaws(628) +ratelaws(635) +ratelaws(700) -ratelaws(711) -ratelaws(726);
    Dspecies(184) = ratelaws(527) -ratelaws(616);
    Dspecies(185) = ratelaws(543) +ratelaws(582) +ratelaws(618) -ratelaws(633) +ratelaws(695) -ratelaws(696) -ratelaws(712) -ratelaws(713) -ratelaws(727) -ratelaws(728) -ratelaws(751) -ratelaws(755);
    Dspecies(186) = ratelaws(545) +ratelaws(583) +ratelaws(619) -ratelaws(634) +ratelaws(698) -ratelaws(699) -ratelaws(714) -ratelaws(715) -ratelaws(729) -ratelaws(730) -ratelaws(752) -ratelaws(756);
    Dspecies(187) = ratelaws(547) +ratelaws(584) +ratelaws(620) -ratelaws(635) +ratelaws(701) -ratelaws(702) -ratelaws(716) -ratelaws(717) -ratelaws(731) -ratelaws(732) -ratelaws(753) -ratelaws(757);
    Dspecies(188) = ratelaws(557) +ratelaws(577) -ratelaws(636) -ratelaws(745);
    Dspecies(189) = ratelaws(558) +ratelaws(578) -ratelaws(637) -ratelaws(746);
    Dspecies(190) = ratelaws(559) +ratelaws(560) +ratelaws(561) +ratelaws(591) +ratelaws(621) -ratelaws(694) -ratelaws(695) +ratelaws(696) -ratelaws(697) -ratelaws(698) +ratelaws(699) -ratelaws(700) -ratelaws(701) +ratelaws(702) -ratelaws(754) -ratelaws(758);
    Dspecies(191) = ratelaws(629) +ratelaws(638) -ratelaws(639) +ratelaws(641) -ratelaws(642) +ratelaws(644) -ratelaws(645) -ratelaws(648) -ratelaws(649) -ratelaws(654) -ratelaws(655) -ratelaws(659) -ratelaws(686) -ratelaws(687) -ratelaws(688) -ratelaws(689) +ratelaws(748) -ratelaws(749) +ratelaws(768) +ratelaws(769) +ratelaws(774) +ratelaws(775) +ratelaws(779) +ratelaws(784) +ratelaws(785) +ratelaws(786) +ratelaws(787) +ratelaws(788) +ratelaws(789) +ratelaws(790) +ratelaws(791) +ratelaws(804) +ratelaws(805) +ratelaws(806) +ratelaws(807) +2.0*ratelaws(808) +ratelaws(812) +ratelaws(814) +ratelaws(819) +ratelaws(820) +ratelaws(826) +ratelaws(835) +ratelaws(836) +ratelaws(837) +ratelaws(838) +ratelaws(839) +ratelaws(840) +ratelaws(841) +ratelaws(842);
    Dspecies(192) = ratelaws(630) -ratelaws(638) -ratelaws(641) -ratelaws(644) -ratelaws(651) -ratelaws(652) -ratelaws(657) -ratelaws(658) -ratelaws(660) +ratelaws(749) +ratelaws(750) +ratelaws(771) +ratelaws(772) +ratelaws(777) +ratelaws(778) +ratelaws(781) +ratelaws(801) +ratelaws(802) +ratelaws(810) +ratelaws(811) +ratelaws(813) +ratelaws(814) +ratelaws(816) +ratelaws(817) +ratelaws(818) +ratelaws(820) +ratelaws(821) +ratelaws(822) +2.0*ratelaws(823) +ratelaws(824);
    Dspecies(193) = ratelaws(631) +ratelaws(639) +ratelaws(640) +ratelaws(642) +ratelaws(643) +ratelaws(645) +ratelaws(646) -ratelaws(647) -ratelaws(649) -ratelaws(650) -ratelaws(652) -ratelaws(661) -ratelaws(690) -ratelaws(691) -ratelaws(692) -ratelaws(693) -ratelaws(747) -ratelaws(748) +ratelaws(767) +ratelaws(769) +ratelaws(770) +ratelaws(772) +ratelaws(780) +ratelaws(783) +ratelaws(792) +ratelaws(793) +ratelaws(794) +ratelaws(795) +ratelaws(796) +ratelaws(797) +ratelaws(798) +ratelaws(799) +ratelaws(800) +ratelaws(802) +ratelaws(843) +ratelaws(844) +ratelaws(845) +ratelaws(846) +ratelaws(847) +ratelaws(848) +ratelaws(849) +ratelaws(850);
    Dspecies(194) = ratelaws(632) -ratelaws(640) -ratelaws(643) -ratelaws(646) -ratelaws(653) -ratelaws(655) -ratelaws(656) -ratelaws(658) +ratelaws(747) -ratelaws(750) +ratelaws(773) +ratelaws(775) +ratelaws(776) +ratelaws(778) +ratelaws(782) +ratelaws(803) +ratelaws(805) +ratelaws(809) +ratelaws(811) +ratelaws(815) +ratelaws(817) +ratelaws(825);
    Dspecies(195) = ratelaws(647) -ratelaws(767) -ratelaws(800) -ratelaws(803) -ratelaws(806);
    Dspecies(196) = ratelaws(648) -ratelaws(768) -ratelaws(801) -ratelaws(804) -ratelaws(807);
    Dspecies(197) = ratelaws(649) -ratelaws(769) -ratelaws(802) -ratelaws(805) -ratelaws(808);
    Dspecies(198) = ratelaws(650) -ratelaws(770) -ratelaws(809) -ratelaws(812);
    Dspecies(199) = ratelaws(651) -ratelaws(771) -ratelaws(810) -ratelaws(813);
    Dspecies(200) = ratelaws(652) -ratelaws(772) -ratelaws(811) -ratelaws(814);
    Dspecies(201) = ratelaws(653) -ratelaws(773) -ratelaws(815) -ratelaws(818);
    Dspecies(202) = ratelaws(654) -ratelaws(774) -ratelaws(816) -ratelaws(819);
    Dspecies(203) = ratelaws(655) -ratelaws(775) -ratelaws(817) -ratelaws(820);
    Dspecies(204) = ratelaws(656) -ratelaws(776) -ratelaws(821);
    Dspecies(205) = ratelaws(657) -ratelaws(777) -ratelaws(822);
    Dspecies(206) = ratelaws(658) -ratelaws(778) -ratelaws(823);
    Dspecies(207) = ratelaws(659) -ratelaws(779) -ratelaws(780) -ratelaws(824);
    Dspecies(208) = ratelaws(660) -ratelaws(781) -ratelaws(782);
    Dspecies(209) = ratelaws(661) -ratelaws(783) -ratelaws(825) -ratelaws(826);
    Dspecies(210) = ratelaws(686) -ratelaws(759) -ratelaws(784) -ratelaws(785) +ratelaws(827);
    Dspecies(211) = ratelaws(687) -ratelaws(761) -ratelaws(786) -ratelaws(787) +ratelaws(829);
    Dspecies(212) = ratelaws(688) -ratelaws(762) -ratelaws(788) -ratelaws(789) +ratelaws(830);
    Dspecies(213) = ratelaws(689) -ratelaws(765) -ratelaws(790) -ratelaws(791) +ratelaws(833);
    Dspecies(214) = ratelaws(690) -ratelaws(760) -ratelaws(792) -ratelaws(793) +ratelaws(828);
    Dspecies(215) = ratelaws(691) -ratelaws(763) -ratelaws(794) -ratelaws(795) +ratelaws(831);
    Dspecies(216) = ratelaws(692) -ratelaws(764) -ratelaws(796) -ratelaws(797) +ratelaws(832);
    Dspecies(217) = ratelaws(693) -ratelaws(766) -ratelaws(798) -ratelaws(799) +ratelaws(834);
    Dspecies(218) = ratelaws(759) -ratelaws(827) -ratelaws(835) -ratelaws(836);
    Dspecies(219) = ratelaws(760) -ratelaws(828) -ratelaws(843) -ratelaws(844);
    Dspecies(220) = ratelaws(761) -ratelaws(829) -ratelaws(837) -ratelaws(838);
    Dspecies(221) = ratelaws(762) -ratelaws(830) -ratelaws(839) -ratelaws(840);
    Dspecies(222) = ratelaws(763) -ratelaws(831) -ratelaws(845) -ratelaws(846);
    Dspecies(223) = ratelaws(764) -ratelaws(832) -ratelaws(847) -ratelaws(848);
    Dspecies(224) = ratelaws(765) -ratelaws(833) -ratelaws(841) -ratelaws(842);
    Dspecies(225) = ratelaws(766) -ratelaws(834) -ratelaws(849) -ratelaws(850);

end


end
