module namelist

  use const_kind_m

! public types
  integer(ki4), parameter :: MAX_NUMBER_OF_POSITIONS=100 !< maximum number of positions allowed
  integer(ki4), parameter :: MAX_NUMBER_OF_PARAMETERS=10 !< maximum number of parameters allowed
  integer(ki4), parameter :: maxinp=100 !< Maximum number of inputs allowed for each variable
  integer(ki4), parameter  :: maximum_number_of_tracks=100 !< maximum number of tracks allowed
  integer(ki4), parameter  :: maximum_number_of_files=100 !< maximum number of files allowed
  integer(ki4), parameter  :: maximum_number_of_panels=100 !< maximum number of panels allowed
!---------------------------------------------------------------------
!> Entries in namelist btreeparameters.
!! These define the fundamentals of the
!! hierarchical data structure (HDS) used in the
!! fieldline tracing, discussed in the ray-tracing context by the
!! preprint at http://arxiv.org/abs/1403.6750
!! The variables are of two types:
!! 1. Allocates an amount of static storage suitable for a geometry
!! with up to 1 million triangles approximately.
!! 2. Determine the type of HDS to be used, normally multi-octree.
!! - Vector components \f$ (x,y,z) \f$ correspond to those of
!! mapped coordinate system, respectively either
!! \f$ (\psi,\theta,\zeta) \f$ or \f$ (R,Z,\xi) \f$.
  type, public :: btreeparameters
   !> EXPERT  size of binary tree, must be positive,
   !! default 2000000
     integer(ki4):: btree_size !< .
   !> EXPERT  size of binary tree pter, must be positive,
   !! default 3
     integer(ki4):: btree_sizep !< .
   !> EXPERT  size of exten array, must be positive,
   !! default 1000
     integer(ki4):: btree_sizee  !< local variable
   !> EXPERT  size of list array hoc, must be positive,
   !! default 2
     integer(ki4):: btree_sizeh !< .
   !> EXPERT  size of list array, must be positive,
   !! default 2000000
     integer(ki4):: btree_sizel !< .
   !> EXPERT  max depth of tree, must be positive,
   !! default 30
     integer(ki4):: btree_depth !< .
   !> \f$ (n_x,n_y,n_z) \f$. Vector of no. of children at
   !! top of multi-octree, components must be positive,
   !! octree at all levels is default, vector (2,2,2),
   !! change only if ttalg=2
     integer(ki4), dimension(3) :: tree_nxyz !< .
   !> EXPERT  \f$ {\bf h}=(h_x,h_y,h_z) \f$ Use this vector of spacings to set
   !! octree quantisation, components must be positive
   !! default spacings (0.1,0.1,0.1)
     real(kr4), dimension(3) :: tree_hxyz !< .
   !> define quantisation algorithm if multi-octree
   !! 0. automatic
   !! 1. default, use hxyz (spacings) to set top
   !! 2. use nxyz, set split at top directly
     integer(ki4):: tree_ttalg !< .
   !> type of tree, must be positive
   !! 1. BSP  (EXPERT)
   !! 2. octree
   !! 3. multi-octree (octree with special top)
     integer(ki4):: tree_type !< .
  end type btreeparameters

!---------------------------------------------------------------------
!> Entries in namelist datvtkparameters.
!! This is the main parameter control for datvtk and made extensive use
!! of by filgeo, when
!! a ctl file (datvtk -c) is used to generate a surface from file(s)
!! containing \f$ R \f$ and/or \f$ Z \f$ point values.
!! (Note that otherwise datvtk
!! acts as a filter, converting file.dat to file_out.vtk.)
!!
!! Remember that angle \f$ \zeta = -\theta \f$ where \f$ \theta \f$ is
!! the more usual cylindrical polar angle. Despite the additional
!! restriction that \f$ |\zeta| < 360 \f$ (in degrees), it is possible
!! to cover any range of angles by suitable choice of limit angles.
!! As an example, to produce a cylinder with a gap of \f$ 15^o \f$ near
!! the negative \f$ y \f$ axis in \f$ x>0 \f$, set
!! start_angle=-280., finish_angle=65..
!! Order of start and finish matters in that the surface normal is
!! outward from the origin if the finish angle is larger than
!! the start angle, inward otherwise.
!! In the 'translate' case, finish_angle is ignored, translation is
!! the direction given by finish_position-start_position and by a
!! distance equal to the length of finish_position-start_position.
!!
!! Note 1) that if variables start_position and finish_position are
!! used to define a surface, then
!! input variable line_divisions must be explicitly set positive. \ns
!! 2) Presently only two types of surfaces may be produced, depending
!! on whether transform_type is set to 'rotate' or 'translate'. 
!! \arg For surfaces formed by rotation, the range of angles is set by
!! start_angle and finish_angle, unless end_angle is nonzero, when
!! the line in \f$ (R,Z) \f$ is swept through all the angles in the geometry \n
!! \arg For surfaces formed by translation,
!! a planar surface is produced by sweeping in the direction
!! of the vector difference of finish_position and start_position.
!! The plane in which the line in \f$ (R,Z) \f$ is regarded as
!! positioned initially is given by start_angle, unless end_angle is nonzero,
!! when either the smallest or largest angle in the geometry is selected
!! dpending whether end_angle=-1 or +1.
  type, public :: datvtkparameters
   !> surface interaction with particles
   !! - 'absorb' first wall (absorbing boundary)
   !! - 'invisi' transparent geometry, totally ignored
   !! - 'skylit' skylight
   !! - 'escape' escape boundary (loss)
   !! - 'errlos' beancan (error loss)
   !! - 'cutout' cutout (expected loss)
   character(len=80) :: description  !< .
  !> If .TRUE., then data file specifies geometry
  !! which is the default.
  !! If filedata becomes .TRUE. 
  !! there needs to be a matching namelist progfiles following
  !! datvtkparameters in the file,
  !! which refers to .dat-style file(s) to define a surface
  logical :: filedata !< .
   !> angle units, either radian(s or degree(s),
   !! default 'degree'
     character(len=20) :: angle_units  !< .
   !> length units, either mm or me(tres)
   !! default 'mm'
     character(len=20) :: length_units  !< .
   !> transform type,
   !! default  'rotate', rotate about central axis of torus,
   !! alternative is 'translate'
     character(len=80) :: transform_type  !< .
   !> starting  toroidal angle \f$ \zeta \f$,
   !! will be overridden when constructing a skylight object
   !! default  0
     real(kr8) :: start_angle !< .
   !> finishing  toroidal angle \f$ \zeta \f$,
   !! will be overridden when constructing a skylight object
   !! default  90
     real(kr8) :: finish_angle !< .
   !> start position 
   !! default is Cartesian \f$ (1,0,0) \f$
     real(kr8), dimension(3) :: start_position !< .
   !> finish position 
   !! default is Cartesian \f$ (1,1,0) \f$
     real(kr8), dimension(3) :: finish_position !< .
   !> When constructing planes using transform_type='translate',
   !! set to +-1 to use largest/smallest angle in geometry
   integer(ki4) :: end_angle !< +-1, use largest/smallest angle in geometry
   !> number of divisions around axis of rotation
   !! or along line of translation,
   !! default  10
     integer(ki4) :: number_of_divisions !< .
   !> coordinate system positive as posang (q.v.),
   !!default and only valid entry is unity for cylindrical polars
     integer(ki4) :: coordinate_system !< .
   !> number of divisions, default  0
   !! Number of points inserted between each point in .dat file, EXCEPT
   !! If line_divisions > number of positions in the .dat file, then
   !! replace those points with positions uniformly spaced on the line between
   !! the first and last \f$(R,Z)$\f$ points in the file
     integer(ki4) :: line_divisions !< .
  end type datvtkparameters

!!---------------------------------------------------------------------
!> Entries in namelist hdsgenparameters.
!! These control subdivision that generates HDS.
!! A box corresponding to a leaf of the binary tree is referred
!! to as a 'bin'.
  type, public :: hdsgenparameters
   !> Stop dividing if less than this number of objects
   !! in 'bin', default 10.
     integer(ki4) :: limit_geobj_in_bin !< .
   !> INERT numerical parameter
     integer(ki4) :: nbins_child  ns_child!< .
   !> Type of geometry
   !! 1. Points (default)
   !! 2. Triangles
     integer(ki4) :: geometrical_type !< .
   !> INERT Positive tolerance <= max_tolerance
     real(kr4) :: min_tolerance !< .
   !> Positive tolerance, default 0.1.
   !! Expand bounding box by this margin in every direction, units of input geometry
     real(kr4) :: max_tolerance !< .
   !> Number of extra bounding cubes, no more than two, default 0
   !! If delta_inner_length>=0 then add this margin to volume,
   !! and if delta_outer_length>=0 add inner+outer_length margin.
   !! If delta_inner_length<0, then use corner values to define and only
   !! produce one cube. N.B. this will be expanded to include all objects
     integer(ki4) :: no_boundary_cubes !< .
   !> vector defining corner with lesser component values
   !! default origin (absolute coordinates)
   !! Only active if any of the components is less than the smallest
   !! corresponding coordinate in the geometry, and then only that coordinates are affected.
     real(kr4), dimension(3) :: lower_corner !< .
   !> vector defining corner with greater component values (absolute coordinates)
   !! Only active if any of the components is greater than the largest
   !! corresponding coordinate in the geometry, and then only those coordinates are affected.
     real(kr4), dimension(3) :: upper_corner !< .
   !> Inner bounding cube separation from geometry, mm,
   !! default -1
     real(kr4) :: delta_inner_length   !< .
   !> Outer bounding cube separation from inner, mm,
   !! default 0
     real(kr4) :: delta_outer_length   !< .
   !> \f$ N_Q \f$, default 1024
     integer(ki4) :: quantising_number !< .
   !> EXPERT Number of objects in tree, default 0, normally determined automatically
     integer(ki4) :: no_geobj_records !< .
   !> If positive, expand global bounding box bin by
   !! margin determined by other settings, default 0
     integer(ki4) :: margin_type !< .
   !> INERT type of scaling of geobj position values
     integer(ki4) :: type_geobj_coord_scaling !< .
  end type hdsgenparameters

!---------------------------------------------------------------------
!> Entries in namelist inputfiles.
!! There is a namelist inputfiles for many SMARDDA
!! programs including hdsgen and move. The namelist
!! specifies names of key input files in each case.
!!
  type, public :: inputfiles
   !> hdsgen only, legacy vtk format geometry data file, default 'null'.
   !! must be ASCII text and note that only a restricted set of objects
   !! such as triangles and tetrahedra is accepted.
     character(len=80) :: vtk_input_file  !< .
   !> move only, special hds format file output by hdsgen to describe
   !! hierachical data structure HDS to accelerate ray-tracing
     character(len=80) :: hds_input_file  !< .
   !> move only, either special qry format file, which has number of
   !! tracks on first line, or \n
   !! legacy vtk format geometry data file, default 'null'.
   !! must be ASCII text and note that only a restricted set of objects
   !! such as triangles and tetrahedra is accepted.
     character(len=80) :: query_input_file  !< .
  end type inputfiles

!---------------------------------------------------------------------
!> Entries in namelist miscparameters.
!! There is a namelist miscparameters many SMARDDA 
!! programs including vtktfm.
!! The namelist was intended as a catch-all
!! for parameters that did not easily fit into any other category, and
!! primarily intended for expert use. However, for vtktfm
!! miscparameters is (mis)used to set critical user quantities.
!!
!! For vtktfm, it helps to know a little history, in that for early
!! ITER work, each panel was regarded as composed of two separate bodies
!! either side of a central attachment (to avoid self-shadowing). For
!! most use, the option 'panel' should be set so that there is a one-to-one
!! correspondence between isolated geometry bodies and panels.
!!
  type, public :: miscparameters
   !> for defining panels and their transforms,
   !! default  'split_panels', i.e. each panel consists of two disjoint
   !! bodies of geometry, unless set to 'panel' when it is assumed that
   !! there is one-to-one correspondence between bodies and panels.
   !! 'tagged' is for EXPERT use only
     character(len=20) :: option !< .
   !> set .TRUE. to read namelist vtktfmparameters
     logical :: new_controls !< .
   !> number \f$ \geq \f$ number of files for which transform defined,
   !! default  1, needs to be \f$ \leq \f$ maximum_number_of_files
     integer(ki4) :: max_number_of_files !< .
   !> units, either radian(s) or degree(s),
   !! default  'radians'
     character(len=20) :: angle_units  !< .
   !> number \f$ \geq \f$ number of panels in geometry,
   !! default  1, needs to be \f$ \leq \f$ maximum_number_of_panels
     integer(ki4) :: max_number_of_panels !< .
   !> number \f$ \geq \f$ number of transforms defined,
   !! default  1, no direct restriction on number
     integer(ki4) :: max_number_of_transforms !< .
   !> vtktfm only, dimension of bods index array, default 1000
     integer(ki4) :: max_bods_index !< .
   !> vtktfm only, used to generate unique bods numbers over many files, default 100
     integer(ki4) :: max_bods_in_file !< .
   !> vtktfm only, bods remain distinct, default .FALSE.
     logical :: preserve_input_structure !< .
   !> number of panels for which bodies defined,
   !! default  0, needs to be \f$ \leq \f$ maximum_number_of_panels
     integer(ki4) :: number_of_panels !< .
   !> number of transforms defined,
   !! default  0, no direct restriction on number
     integer(ki4) :: number_of_transforms !< .
  end type miscparameters

!---------------------------------------------------------------------
!> Entries in namelist numericalparameters
  type, public :: numericalparameters
  integer(ki4) :: geometrical_type !< type of geometry, default  1
  real(kr4) :: min_tolerance !< numerical parameter, default  1.0e-3
  real(kr4) :: max_tolerance !< numerical parameter, default  1.0e-1
  integer(ki4) :: quantising_number !< local variable, default 1024
  integer(ki4) :: no_geobj_records !< local variable, default  0
  integer(ki4) :: margin_type !< local variable, default  0
  !> type of scaling of geobj x value, default 2,
  !! allows for offset
  integer(ki4) :: type_geobj_coord_scaling !< .
  !> content of query file, which need not have full vtk
  !! header, only number of tracks on first line
  !! 0. Default, alternate start and end points
  !! 1. Alternate start points and velocity values
  !! 2. Velocity field is separate VECTORS field_name in .vtk file
  integer(ki4) :: query_file_content ! < .
  !> nominal timestep if velocity in vtk test file, default 1
  real(kr8) :: nominal_dt !< .
  !> actual name of velocity vector field in vtk test file,
  !! default 'Vel'
  character(len=80) :: field_name !< .
  end type numericalparameters

!---------------------------------------------------------------------
!> Entries in namelist panelarrayparameters.
!! This assigns transforms to bodies, strictly objects defined by
!! the numbers defined by process_by_name in namelist
!! vtktfmparameters (if used). It is easiest to understand
!! if every body appears in the
!! panel_bodies list and has a corresponding transform,
!! unless the 'split' option is set when each transform corresponds
!! to two adjacent bodies in the list. However, the software will by
!! default enter all bodies into a list of length up to
!! maximum_number_of_panels and assign the
!! identity transformation to each.
!!
  type, public :: panelarrayparameters
   !> bodies defining the geometry,
   !! if no bodies are listed, then all are considered
     integer(ki4), dimension(2*maximum_number_of_panels) :: panel_bodies !< .
   !> number of transform to apply,
   !! default -1 gives error so must be changed
     integer(ki4), dimension(2*maximum_number_of_panels) :: panel_transform !< .
  end type panelarrayparameters

!---------------------------------------------------------------------
!> Entries in namelist plotselections.
!! There is a namelist plotselections for each of the major component
!! programs, datvtk, hdsgen and move (although required for datvtk
!! if a .ctl file is used, the namelist controls have no effect). The plots that
!! can be produced vary with each program, but the following conventions
!! are observed:
!! 1. All variables are logicals, and must be set .TRUE. to produce a
!! data file for plotting, or be omitted, or set .FALSE. (otherwise an error will
!! terminate execution).
!! 2. All variable names begin with 'plot_', then followed by an identifier
!! for the particular program if a vtk file is wanted, thus
!!     - hds or geo for hdsgen
!!     - mov or geo for move
!!     - if gnu appears, a file for gnuplot is produced
!! 3. All vtk files (for visualisation using ParaView) contain only surface
!! information, unless the string
!! 'vol' appears following the program identifier, or they show the HDS
!! which is inherently a 3-D wireframe.
!! 4. The coordinate system used is indicated by the last letter before
!! the .vtk filetype, viz:
!!     - x or m for Cartesians (mm)
!!     - q for quantised coordinates, the system in which lines are actually
!! tracked, typically ranging over 0 to 1024.
!! 5. There are also names which are a legacy of previous developments,
!! if they are set .TRUE., a warning will be issued.
!! 6. Corresponding output files to variable plot_xxx have names root_xxx.vtk
!! or root_xxx.gnu or root_gnu(_opt).gnu, where root.ctl is the
!! control file used to generate them
!! 7. If it is not sensible to produce a plot file, its output will be
!! suppressed
  type, public :: plotselections
   !> hdsgen and move only, default .false., above conventions apply.
     logical :: plot_hdsm !< .
   !> hdsgen and move only, default .false., above conventions apply.
     logical :: plot_hdsq !< .
   !> hdsgen and move only, default .false., plot geometry
   !! in quantised coordinates
     logical :: plot_geobjq !< .
   !> hdsgen and move only, default .false., plot points of geometry
     logical :: plot_geoptq !< .
   !> move only, default  .true., plot points of geometry
     logical :: plot_geopt !< .
   !> move only, default  .true., plot tracks including 
   !! intersection points with geometry
     logical :: plot_moves !< vtk plot selector, default  .true.
   !> move only, default  .false., plot tracks including 
   !! intersection points with geometry
   !! in quantised coordinates
     logical :: plot_movesq !< .
  end type plotselections

!---------------------------------------------------------------------
!> Entries in namelist positionparameters.
!! These describe a transform of position in a 3-D coordinate system.
!!
  type, public :: positionparameters
   !> transformation 3x3 matrix,
   !! default is identity matrix
     real(kr4), dimension(3,3):: position_matrix  !< .
   !> transformation scale 3-vector,
   !! default \f$ (1.,  1.,  1. ) \f$
     real(kr4), dimension(3):: position_scale  !< .
   !> transformation offset 3-vector,
   !! default \f$ (0.,  0.,  0. ) \f$
     real(kr4), dimension(3):: position_offset  !< .
   !> transformation type,
   !! default 0 (identity),
   !! - 1. scale - scale in Cartesian coordinates
   !! - 2. offsca - scale and offset  in Cartesians
   !! - 6. poltil - poloidal tilt, rotate about an axis parallel to the
   !! toroidal direction through the centre of the body
   !! - 7. tortil - toroidal tilt, rotate about an axis in the
   !! vertical direction through the centre of the body
   !! - 12. torrot - rotate in toroidal direction
   !! - 22. polrot - rotate in poloidal direction
   !! - 42. radial - displace towards a point \f$ (R,Z) \f$, which will
   !! normally be the centre of the torus (i.e. displace in minor radius),
   !! negative displacements imply outwards translation. The first entry is
   !! the displacement, and the next two give the point coordinates
     integer(ki4):: position_transform  !< .
   !> describes transformation type,
   !! default 'cartesian_scale' (inert)
     character(len=20) :: transform_desc !< .
   !> identifies transformation,
   !! default 'unit' (inert)
     character(len=20) :: transform_id !< .
  end type positionparameters

!---------------------------------------------------------------------
!> Entries in namelist progfiles.
!! These set the file(s) used by datvtk -c.
!! Either r_input_file and z_input_file  should contain
!! equal numbers of \f$ R \f$ and \f$ Z \f$ coordinates respectively, one per line, with
!! a null value for rz_input_file, or rz_input_file  should
!! contain the \f$ (R,Z) \f$ positions, one pair of numbers per line, and the other two
!! filenames should be 'null'.
!!
!! The file format used here is referred to as the 'gnu dat' format, to
!! distinguish it from the NASTRAN .dat format, as it is compatible with the
!! gnuplot plotting package. It allows for a single header line beginning '#'
!! with a closely defined layout to describe the content.
!!
  type, public :: progfiles
   !> position coordinate 1 \f$ (R) \f$ data input file,
   !! default 'null'
     character(len=80) :: r_input_file  !< .
   !> position coordinate 2 \f$ (Z) \f$ data input file,
   !! default 'null'
     character(len=80) :: z_input_file !< .
   !> position coordinate 1 and 2 \f$ (R,Z) \f$ data input file,
   !! default 'null'
     character(len=80) :: rz_input_file  !< .
  end type progfiles

!---------------------------------------------------------------------
!> Entries in namelist vtkfiles.
!! This sets the names of the vtk format files to be used by vtktfm,
!! and how many copies in total are needed of each file.
!! Note that if a file is to be copied, it is assumed that it contains
!! only one body (use datvtk -xxm to remove body identifiers from dat files
!! when converting to vtk).
!! The copies of the first file in the vtk_input_file list are numbered
!! 101, 102, 103, ... etc, the second are numbered 201, 202, 203, ...etc.
!! (so that if 4 occurrences of the geometry in file 2 are needed,
!! set the number of copies to 4, and reference them as 201-204).
!!
  type, public :: vtkfiles
   !> vtk format geometry data input files,
   !! default 'null' gives error so must be changed
     character(len=80), dimension(maximum_number_of_files) :: vtk_input_file !< .
   !> number of copies needed of each corresponding input file,
   !! default 1 for each entry in vtk_input_file list
     integer(ki4), dimension(maximum_number_of_files) :: number_of_copies !< .
   !> labels vtk format geometry data input file (not used),
   !! default 'null'
     character(len=80), dimension(maximum_number_of_files) :: vtk_label  !< .
   !> numeric label for input file (not used),
   !! default 0
     integer(ki4), dimension(maximum_number_of_files) :: vtk_label_number !< .
   !> vtk format geometry data output file,
   !! defaults to root of ctl file (plus vtk suffix)
     character(len=80) :: vtk_output_file  !< .
  end type vtkfiles

!---------------------------------------------------------------------
!> Entries in namelist vtktfmparameters.
!!
!! It helps to know a little history, in that for early
!! ITER work, each panel was regarded as composed of two separate bodies
!! either side of a central attachment (to avoid self-shadowing). For
!! more recent use, the option 'panel' should be set so that there is a one-to-one
!! correspondence between isolated geometry bodies and panels. Further,
!! it became possible for the code to distinguish any integer-labelled attribute
!! in the input, e.g. using set process_by_name=Surf for surfaces, hence the
!! term 'bods' is now preferred.
!!
!! There are important differences in behaviour depending on whether
!! the main geometry input is a single file or not. In the former case
!! it is possible to split the file by a named attribute (typically 'Body')
!! and/or to displace each body by a numbered transform. Provided
!! there are fewer than max_bods_in_file=100 bodies, then it is not necessary to specify all
!! the bodies in panel_bodies, those omitted will be assigned the
!! identity transform. Note the variable make_same must be explicitly set
!! .FALSE., also then geometry which is not assigned to a body
!! or surface (Body "0" or zero) cannot be transformed.
!!
!! If there is more than one input file,
!! then by default each file is treated as a single object with a body with
!! number greater than 100. The numbering convention should be obvious
!! but might take some unravelling, especially as there are bugs in vtktfm (TO DO).
!! The first file is numbered 101, the second 201, the third 301, etc.
!! The number of copies may be specified for each file, default 1
!! (i.e. a single occurrence of the geometry). Other copies of
!! file 1 are numbered 102, 103, 104, etc, of file 2 are numbered 202,
!! 203,204, etc. max_bods_in_file may be freely changed on input to
!! accommodate more complicated geometries. If preserve_input_structure
!! is set .TRUE., then the bods in the input files are (re)numbered
!! consecutively starting with unity, but there must not be more than
!! max_bods_index bods. \n
!! WARNING (1) unless the bodies are numbered consecutively setting
!! preserve_input_structure=.TRUE. will usually cause a crash error.\n
!! (2) Most
!! combinations of preserve_input_structure, extract, split_file and
!! make_same are untested and may well also cause vtktfm to crash.
!!
  type, public :: vtktfmparameters
   !> Option for defining panels and their transforms
   !! 'split' default
   !! 'panel' bodies in pairs to define panel
   !! 'tagged' bodies tagged with their transform number
     character(len=20) :: option !< .
   !> number \f$ \geq \f$ number of files for which transform defined,
   !! default  1, needs to be \f$ \leq \f$ maximum_number_of_files
     integer(ki4) :: max_number_of_files !< .
   !> units, either radian(s) or degree(s),
   !! default  'radians'
     character(len=20) :: angle_units  !< .
   !> number \f$ \geq \f$ number of panels in geometry,
   !! default  1, needs to be \f$ \leq \f$ maximum_number_of_panels
     integer(ki4) :: max_number_of_panels !< .
   !> number \f$ \geq \f$ number of transforms defined,
   !! default  1, no direct restriction on number
     integer(ki4) :: max_number_of_transforms !< .
   !> number of panels for which bodies defined,
   !! default  0, needs to be \f$ \leq \f$ maximum_number_of_panels
     integer(ki4) :: number_of_panels !< .
   !> number of transforms defined,
   !! default  0, no direct restriction on number
     integer(ki4) :: number_of_transforms !< .
   !> set .TRUE. to split file by attribute name, default .FALSE.
     logical :: split_file !< .
   !> set .TRUE. to make all attribute names the same, default .TRUE.
     logical :: make_same !< .
   !>  value for attribute to take (e.g. body number), default 1
     integer(ki4) :: same_value !< .
   !>  name of attribute to be split / homogenised, default Body
     character(len=80) :: process_by_name !< .
     logical :: paneltfm !< apply transform only if .TRUE. (default)
     logical :: extract !< extract objects according to key and limits, default .FALSE.
     character(len=80) :: extract_key !< key for extraction, default 'toroidal'
   !> dimension of bods index array, default 1000
     integer(ki4) :: max_bods_index !< .
   !> used to generate unique bods numbers over many files, default 100
     integer(ki4) :: max_bods_in_file !< .
   !> bods remain distinct, default .FALSE.
     logical :: preserve_internal_structure !< .
     real(kr8), dimension(2) :: plasma_centre !< centre of discharge in \f$ (R,Z) \f$
     real(kr8) :: minimum_angle !< minimum angle for extraction, default -22.5
     real(kr8) :: maximum_angle !< maximum angle for extraction, default 22.5


  end type vtktfmparameters

end module namelist
