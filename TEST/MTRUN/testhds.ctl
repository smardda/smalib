&inputfiles
 vtk_input_file='sduct.vtk',
/
&hdsgenparameters
   geometrical_type=2,
   min_geobj_in_bin=20,
   limit_geobj_in_bin=80,
   margin_type=2,
/
&btreeparameters
   tree_type=3,
   tree_ttalg=2,
   tree_nxyz=8,1,2,
/
&positionparameters
   position_transform=1,
/
&plotselections
   plot_hdsm = .true.,
   plot_hdsq = .true.,
   plot_geobjq = .true.,
   plot_geoptq = .true.,
/
