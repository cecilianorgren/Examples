%% Go to your working folder
cd /Users/Cecilia/Data/Cluster/20080422
% or make it:
% mkdir /Users/Cecilia/Data/Cluster/2008Mar22

%% Define time intervals
% Give start and stop times
tint = irf.tint('2008-04-22T17:50:00.00Z/2008-04-22T18:15:00.00Z')
% or give start time and duration in seconds
tint = irf.tint('2008-04-22T17:50:00.00Z',25*60)

%% Download data, if you don't have already
% This will download data into a subfolder of your current directory, so
% check where you are before you download it so you don't have to move it
% around.
% Do only once, otherwise the loading functions will find several files and
% ask you each time which one to download.
% I sort my events like this:
% ../Data/Cluster/Date1/CAA/
% ../Data/Cluster/Date2/CAA/
% ../Data/Cluster/Date3/CAA/
% etc...
% ls /Users/Cecilia/Data/Cluster/

caa_download(tint,'C?_CP_CIS_HIA_ONBOARD_MOMENTS');
caa_download(tint,'C?_CP_FGM_FULL');

%% Load data into a matrix with time in first column and date in subsequent columns
% In this case, the data (magnetic field,temperature, density, velocity,
% etc) and the time are the same data type. This has caused some precision
% problem with the time in prior missions.
gseB1mat  = c_caa_var_get('B_vec_xyz_gse__C1_CP_FGM_FULL','mat');
Ti1mat    = c_caa_var_get('temperature__C1_CP_CIS_HIA_ONBOARD_MOMENTS','mat');
ni1mat    = c_caa_var_get('density__C1_CP_CIS_HIA_ONBOARD_MOMENTS','mat');
gseVi1mat = c_caa_var_get('velocity_gse__C1_CP_CIS_HIA_ONBOARD_MOMENTS','mat');
if 0 % Use c_eval to load data for all four spacecraft 
  %%
  sclist = 2:4;
  c_eval('gseB?mat  = c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');',sclist)
  c_eval('Ti?mat    = c_caa_var_get(''temperature__C?_CP_CIS_HIA_ONBOARD_MOMENTS'',''mat'');',sclist)
  c_eval('ni?mat    = c_caa_var_get(''density__C?_CP_CIS_HIA_ONBOARD_MOMENTS'',''mat'');',sclist)
  c_eval('gseVi?mat = c_caa_var_get(''velocity_gse__C?_CP_CIS_HIA_ONBOARD_MOMENTS'',''mat'');',sclist)
end

%% Load data into a TSeries class object
gseB1  = c_caa_var_get('B_vec_xyz_gse__C1_CP_FGM_FULL','ts');
Ti1    = c_caa_var_get('temperature__C1_CP_CIS_HIA_ONBOARD_MOMENTS','ts');
ni1    = c_caa_var_get('density__C1_CP_CIS_HIA_ONBOARD_MOMENTS','ts');
gseVi1 = c_caa_var_get('velocity_gse__C1_CP_CIS_HIA_ONBOARD_MOMENTS','ts');

%% Check the difference and compare them
whos gseB1 gseB1mat Ti1 Ti1mat ni1 ni1mat
irf_plot({gseB1,gseB1mat},'comp') % plots different components in different panels
irf_legend({'TSeries','matrix'},[0.95 0.95])

%% Make time series from matrix double
c_eval('gseB?mat2ts=irf.ts_vec_xyz(gseB?mat(:,1),gseB?mat(:,2:4));',1)

%% Make a simple overview plot
irf_plot({gseB1,gseVi1,ni1,Ti1}) % plots different components in different panels
% ylabels are not very pretty, but they display the information contained
% within the TSeries object.

%% Check what information we have there
% gseB1.tensorBasis
% gseB1.coordinateSystem
% gseB1.userData.GlobalAttributes.OBSERVATORY_DESCRIPTION
% gseB1.x

%% Make nicer figure
% Initialize plot
h = irf_plot(4); % four panels, h contains axes handles

hca = irf_panel('B');
irf_plot(hca,gseB1);
hca.YLabel.String = 'B (nT)';

hca = irf_panel('Vi');
irf_plot(hca,gseVi1);
%hca.YLabel.String = 'v_{i,1} (km/s)'; 
ylabel(hca,'v_{i} (km/s)','interpreter','tex')

hca = irf_panel('ni');
irf_plot(hca,ni1);
%hca.YLabel.String = 'n_i (cm^{-3})';
ylabel(hca,'n_i (cm^{-3})','interpreter','tex')

hca = irf_panel('Ti');
% I want to change unist of temperature from MK to eV. This can be made
% directly on the TSeries object. With matrix, we had to do.
% [Ti1(:,1) Ti1(:,2)*units.MK*units.kB/units.e]
irf_plot(Ti1*units.MK*units.kB/units.e)
%hca.YLabel.String = 'T_i (eV)';
ylabel(hca,'T_i (eV)','interpreter','tex')

%% Adjust the figure a little bit
h(1).Title.String = 'Cluster 1 (GSE)';
tintZoom = irf.tint('2008-04-22T18:00:00.00Z/2008-04-22T18:10:00.00Z');
irf_zoom(h,'x',tintZoom); % zoom in figure
irf_zoom(h,'y'); % adjusts y limits so it's good for the new time interval
irf_plot_axis_align % align labels on y axis

%% Marklight some time interval
tintZoom = irf.tint('2008-04-22T18:04:55.00Z/2008-04-22T18:10:05.05Z');
irf_pl_mark(h,tintMark.epochUnix')

%% I want to check the parallel ion velocity
parVi1 = gseVi1.dot(gseB1/gseB1.abs);

%% Resampling gseB1 to the time of gseVi1
parVi1 = gseVi1.dot(gseB1.resample(gseVi1)/gseB1.abs);
% Note: resample didn't need to be applied to gseB1.abs, since the
% .resample has been built into the TSeries.mrdivide (same as '/'). The
% package is under development...

%% Check B unit vector
hatB = gseB1.resample(gseVi1)/gseB1.abs;
irf_plot(hatB)
% irf_plot({hatB.x,hatB.y,hatB.z,hatB.abs},'comp')

%% Almost same figure as before, but including the parallel ion velocity
% Initialize plot
h = irf_plot(5); % four panels, h contains axes handles

hca = irf_panel('B');
irf_plot(hca,{gseB1.x,gseB1.y,gseB1.z,gseB1.abs},'comp');
hca.YLabel.String = 'B (nT)';
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.02 0.15])

hca = irf_panel('Vi');
irf_plot(hca,gseVi1);
%hca.YLabel.String = 'v_{i,1} (km/s)'; 
ylabel(hca,'v_{i} (km/s)','interpreter','tex')
irf_legend(hca,{'v_x','v_y','v_z'},[0.02 0.15])

hca = irf_panel('parallel Vi');
irf_plot(hca,gseVi1.dot(gseB1.resample(gseVi1)/gseB1.abs));
%hca.YLabel.String = 'v_{i,1} (km/s)'; 
ylabel(hca,'v_{i,||} (km/s)','interpreter','tex')

hca = irf_panel('ni');
irf_plot(hca,ni1);
%hca.YLabel.String = 'n_i (cm^{-3})';
ylabel(hca,'n_i (cm^{-3})','interpreter','tex')

hca = irf_panel('Ti');
% I want to change unist of temperature from MK to eV. This can be made
% directly on the TSeries object. With matrix, we had to do.
% [Ti1(:,1) Ti1(:,2)*units.MK*units.kB/units.e]
irf_plot(Ti1*units.MK*units.kB/units.e)
%hca.YLabel.String = 'T_i (eV)';
ylabel(hca,'T_i (eV)','interpreter','tex')


% Adjust the figure a little bit
h(1).Title.String = 'Cluster 1 (GSE)';
tintZoom = irf.tint('2008-04-22T18:00:00.00Z/2008-04-22T18:15:00.00Z');
irf_zoom(h,'x',tintZoom); % zoom in figure
irf_zoom(h,'y'); % adjusts y limits so it's good for the new time interval
irf_plot_axis_align % align labels on y axis

%% Other operations
gseVi1+gseB1
gseVi1+gseB1.resample(gseVi1)
gseVi1.cross(gseB1)
gseVi1.cross(gseB1.resample(gseVi1))
gseVi1.cross(gseB1.resample(gseVi1)).tlim(tintZoom).x

%% Loading data with MMS
% Show database structure
% Initialize database
mms.db_init('local_file_db','/Volumes/Samsung/data');

%% For loading, see 
% open mms_2015Oct16.load_srvy_data
% open mms_2015Oct16.load_brst_data
eval('1+1')
c_eval('1+?',2)
c_eval('str?=''hej'';',3); str, whos str
c_eval('!+?',4,5)

%% Show MMS example plots
load /Users/Cecilia/Data/MMS/2015Oct16/workspace2015Nov09.mat
mms_2015Oct16.make_omni_and_perpar

%% Example plot 1
%mms_2015Oct16.plot_psd_projection
% toplot = 4; % knee in gyrotropic perp, 6 pitchangles
% toplot = 3; % 
% toplot = 6; % mva: 3pa, 3 projections

%% Example plot 2
% Define time interval
tint = irf.tint('2015-10-16T10:33:45.00Z',0.04);

% Get mean magnetic field direction
c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint).time).data);',ic); 
hatB0 = double(irf_norm(B0));

% Initialize figure
for ii = 1:4; h(ii) = subplot(2,2,ii); end
isub = 1;

% Plot particle distribution for all direction but a single energy
hca = h(isub); isub = isub + 1;
mms.plot_skymap(hca,desDist1,'tint',tint,'energy',70,'vectors',{hatB0,'B'});

% Plot particle distribution for all direction but a single energy
hca = h(isub); isub = isub + 1;
mms.plot_skymap(hca,desDist1,'tint',tint,'energy',180,'vectors',{hatB0,'B'});

% Plot projection onto a plane perpendicular to B
hca = h(isub); isub = isub + 1;
mms.plot_projection(hca,desDist1,'tint',tint,'xyz',[1 0 0; 0 1 0; hatB0],'elevationlim',20,'vlim',15000);

% Plot particle distribution for pitch angles 0 90 and 180
hca = h(isub); isub = isub + 1;
mms.plot_cross_section_psd(hca,desDist1,dmpaB1brst,'tint',tint,'scPot',P1,'energies',[70 180])

%% Other time interval
tint = irf.tint('2015-10-16T10:33:30.25Z',0.04); 

if 0 % get mva direction for projection plot
  %%
  % Automatic results made from gui above
  ic = 1:4;
  tint_mva = irf.tint('2015-10-16T10:33:13.227751708Z/2015-10-16T10:33:38.076784912Z');
  c_eval('[out?,l?,mva?]=irf_minvar(dmpaB?.tlim(tint_mva));',ic)
  mva = mva3;
  % rotate around x-axis 
  turn_angle = 0; % degrees
  tmpv = mva;
  mva(2,:) = tmpv(2,:)*cosd(turn_angle) + tmpv(3,:)*sind(turn_angle);
  mva(3,:) = tmpv(3,:)*cosd(turn_angle) - tmpv(2,:)*sind(turn_angle);
  % v is 3x3 matrix:
  %v1: [0.0906 0.0896 0.9918] - first row
  %v2: [0.3510 -0.9349 0.0524] - second row
  %v3: [0.9320 0.3434 -0.1161] - third row   
end

% Get mean magnetic field direction
c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint).time).data);',ic); 
hatB0 = double(irf_norm(B0));
energy1 = 50;
energy2 = 180;

% Initialize figure
for ii = 1:4; h(ii) = subplot(2,2,ii); end
isub = 1;

% Plot particle distribution for all direction but a single energy
hca = h(isub); isub = isub + 1;
mms.plot_skymap(hca,desDist1,'tint',tint,'energy',energy1,'vectors',{hatB0,'B'});

% Plot particle distribution for all direction but a single energy
hca = h(isub); isub = isub + 1;
mms.plot_skymap(hca,desDist1,'tint',tint,'energy',energy2,'vectors',{hatB0,'B'});

% Plot projection onto a plane perpendicular to B
hca = h(isub); isub = isub + 1;
mms.plot_projection(hca,desDist1,'tint',tint,'xyz',mva,'elevationlim',20,'vlim',15000,'vectors',{hatB0,'B'},'vlabel',{'v_L','v_M','v_N'});

% Plot particle distribution for pitch angles 0 90 and 180
hca = h(isub); isub = isub + 1;
mms.plot_cross_section_psd(hca,desDist1,dmpaB1brst,'tint',tint,'scPot',P1,'energies',[energy1 energy2])
