&inputfiles
vtk_input_file='div_ansysmesh.vtk',
/
&miscparameters
/
&plotselections
plot_vtk = .true.,
/
&geofilparameters
absorber_objects=2,
beancan_objects=0,
/
&datvtkparameters
description='absorb',
filedata=.true.,
transform_type = 'translate',
end_angle = -1,
start_position = 0,0,0,
finish_position = 0,0,2,
line_divisions = 10
/
&progfiles
rz_input_file='jetbot1.dat',
/
&datvtkparameters
description='absorb',
filedata=.true.,
transform_type = 'translate',
end_angle = 1,
start_position = 0,0,0,
finish_position = 0,0,2,
line_divisions = 10
/
&progfiles
rz_input_file='jetbot2.dat',
/
