!>\mainpage SMARDDA-LIB Ray tracing User Manual
!! 
!! 
!! \section sec1 What is SMARDDA-LIB?
!! 
!! SMARDDA-LIB (SMALIB for short) is a set of object-oriented 
!! Fortran-95 modules for ray, particle and fieldline tracing 
!! using the SMARDDA ray-tracing algorithm. SMARDDA is a published
!! algorithm for triangulated surface geometries, designed to handle
!! short tracks as well as 'long' rays efficiently. The central
!! geometry format is the legacy vtk format, with the smint tool
!! provided for basic geometry manipulations including conversions
!! from other mesh formats.
!!
!! Three basic mesh manipulation codes also appear in SMARDDA-LIB, namely
!! datvtk, geofil and vtktfm, together with code to produce a data structure
!! needed to accelerate ray-tracing (hdsgen) and a simple test harness (move).
!! Codes vtktfm, hdsgen and move need ParaView for visualisation.
!!
!! SMARDDA-LIB is offered under the LGPL 3 licence to
!! encourage developments using SMARDDA modules both
!! inside and outside the fusion community. (Within the
!! community a range of magnetic fusion physics codes have
!! been built, notably SMARDDA-PFC as the kernel of the ITER
!! power deposition package SMITER, see \ref sec1resource)
!! 
!! Pronunciation SMARDDA is pronounced smarter.
!! 
!! SMARDDA-LIB technical features
!! 
!! \arg Hybrid algorithm for long and short rays \n
!!   -Short rays optimal for discretising curved particle tracks \n
!!   -Long rays optimal for long neutral particle tracks \n
!!   -Published in detail and extensively used in fusion physics
!! \arg Uses surface geometry described by legacy vtk format \n
!!   -Primarily triangles \n
!!   -Utility (datvtk) to convert from dat format output by CATIA-TM \n
!!   -Other utilities based on ParaView to convert STL and other meshes \n
!!   -Mesh manipulation by vtktfm \n
!!   -Capability to construct elementary surfaces using geofil
!! \arg Toroidally periodic feature \n
!!   -Only necessary to specify repeating fraction of geometry \n
!! \arg Written in published Fortran-95 style \n
!! \arg On-going development programme \n
!!   -MPI parallelism \n
!!   -ANSYS integration  \n
!!   -GPU implementation  \n
!!   -Light-weight provenance capture \n
!!   -Developer documentation
!! 
!>
!! 
!! \subsection sec1coord Coordinate Systems
!! So Many Accurate Real Descriptions of DAta
!! \arg Accurate first wall geometry originates in 
!! Computer-Aided-Design (CAD) geometry database,
!! coordinates (X,Y,Z) in millimetres (mm).
!! \arg Conversion to cylindrical polars/toroidal
!! coordinates. (TO DO)
!! 
!! 
!! 
!>
!! \image html xy.png
!! \image latex xy.eps "Cartesian XYZ coordinates" width=5cm
!! Cartesian \f$ XYZ \f$  coordinates
!! 
!>
!! \image html rz.png
!! \image latex rz.eps "Cylindrical polars" width=5cm
!! Cylindrical polars \f$ RZ\zeta \f$, N.B. not \f$ R\phi Z \f$. Distances in
!! metres (m), angles in radians. Also \f$ RZ\xi \f$.
!>
!! \image html circ.png
!! \image latex circ.eps "Toroidal coordinates" width=5cm
!! Toroidal coordinates \f$ r\theta\zeta \f$, angles in radians. Angle
!! measured from vertical through X-point. Also \f$ r\theta\xi \f$.
!!
!>
!! Angle \f$ \xi \f$ is \f$ \zeta \f$ scaled so that \f$ \xi \f$ runs from 
!! \f$ -\pi \f$  to \f$ +\pi \f$ in a toroidal segment, which is a
!! fraction \f$ 1/N_{\zeta p} \f$ of
!! the device geometry (often \f$ N_{\zeta p} = 18 \f$).
!! 
!! 
!! \subsection sec1surf Defining Surfaces Accurately
!! So Many AccuRate Definitions of surface DAta
!! \arg CAD databases contain curves and surfaces defined
!! as NURBS, a special sort of rational polynomial that
!! enables exact representation of for example the conic sections.
!! \arg Intersections of NURBS surfaces are only
!! approximately represented by NURBS curves,
!! implying small gaps.
!! \arg Small gaps avoided by use of triangular faceting of
!! every surface. \n
!! Lots of triangles to test for intersection with rays or curve
!! segments, unless you get SMARDDA.
!! 
!! \subsection sec1alg SMARDDA Algorithm
!! Reduces number of triangles to test by:
!! \arg Recursive division of 3-D space into cuboidal boxes, aka octree.
!! Recursion stops when number of triangles in a box aka bin
!! falls below a preset value. (Always know which bin your fieldline
!! segment is in.)
!! \arg Track segment from box to box by binary tricks, thanks to use of
!! quantised coordinates. The quantum lengths \f$ h_i \f$ are chosen so
!! that \f$ L_i=1024h_i \f$ or similar power of two. (The biggest box has
!! size \f$ L_1L_2L_3 \f$. There is a single precision restriction \f$ 3\times 10<32 \f$)
!! \arg Jargon -  an octree is an example of an Hierarchical Data
!! Structure (HDS), designed to Help Do Surfaces quickly.
!! 
!! \subsection sec1morecoord Coordinates Summary
!! SuMmary of AccuRate Descriptions of coorDinAtes
!! 1. Cartesian \f$ XYZ \f$ mm from CAD
!! 2. Cylindrical polars \f$ RZ\zeta \f$ metres and radians (\f$ RZ\xi \f$) (TO DO)
!! 3. Toroidal polars \f$ r\theta\zeta \f$ metres and radians (\f$ r\theta\xi \f$) (TO DO)
!! 4. Quantised coordinates \f$ xyz \f$ or \f$ 123 \f$, versions of 2 and 3 above
!! 
!! \subsection sec1work Workflow
!>
!! \image html flomm.png
!! \image latex flomm.eps "SMARDDA-MOVE flow chart" width=10cm
!!
!! 
!! Steps MAking ReaDy for geometry DAta
!! 1. CAD database input -> CAD defeaturing of part body -> CAD
!! "healing" -> Meshing -> Mesh output of part body. \n
!! (.CATpart file
!! -> .dat file, preferably produced by specially trained Design
!! Office staff)
!! 2. Conversion of .dat format (proprietary MSC Nastran or
!! a project specific 'gnu dat' format \ref sec2formats) to
!! legacy .vtk format (OpenSource) by SMALIB datvtk.
!! 3. Combination of part bodies to form PFC geometry under study, by
!! SMARDDA-PFC vtktfm.
!! 4. Conversion of STL and obj files to legacy .vtk format (OpenSource) by
!! ParaView, scripts stl2vtk and obj2vtk provided.
!! 5. Use of ParaView "Sources" menu  to produce simple objects in legacy .vtk,
!! format, other objects by meshing constant surfaces of 3-D scalar data.
!! 6. Check using ParaView that surface normals point outward, and if
!! necessary correct (a body at a time) using datvtk.
!!
!! Aside: the UKAEA approach: .CATpart .x_t or .step file imported by
!! CADfix-TM on Windows for user-controlled defeaturing, healing,
!! orient normals and meshing. cfmcn(p) converts CADfix proprietary format to .vtk
!! 
!! 
!! \subsection sec1cli Command line interface
!!
!! (If SMARDDA-LIB has been installed in directory ~/smardda/smardda-lib, the SMALIB
!! codes may be added to the user's path with the shell commands: \n
!! export SMALIB_DIR=~/smardda/smardda-lib;source $SMALIB_DIR/Extras/setup
!! )
!!
!! \arg The standard way to execute codes is \n
!! \Emph{code .ctl file}, \n
!! e.g. \Emph{hdsgen wall}, where file wall.ctl contains
!! inputs controlling hdsgen run, except
!! \arg datvtk was designed as a filter, so instead \n
!! datvtk -[opt] ".dat file", or datvtk -c ".ctl file" \n
!! In the latter case, the .ctl file
!! contains instructions to generate surfaces by rotating
!! straight lines defined by their end-points, and the points
!! may be defined by files in the gnu dat format (\ref sec2formats).
!! 
!! Aside: \Emph{opt=ABC}, where \Emph{A=v} -> .vtk input; \Emph{B=o} -> orient
!! triangles relative to first in file, \Emph{B=s}->shell tetrahedra (BROKEN, use ParaView
!! instead);\Emph{B=f} -> flip triangle orientation; \Emph{C=d} -> divide into
!! separate block body files, \Emph{C=m} ->suppress block
!! information, e.g. datvtk -xom FW34.dat (x as spacer)
!! 
!! The .ctl files
!! contain Fortran namelists which must appear in a
!! strict order although the variables within a namelist may be
!! (re)defined in any order.
!! Variables are documented using doxygen 1.8, which produces
!! .html files for input to firefox or more generally xdg-open,
!! which are linked to the \Bold{smint} script \ref sec2smint.
!! 
!! Example .ctl file input to geofil \n
!! &inputfiles \n
!! vtk_input_file='sduct.vtk', \n
!! / \n
!! &miscparameters \n
!! / \n
!! &plotselections \n
!! plot_vtk = .true., \n
!! / \n
!! &geofilparameters \n
!! absorber_objects=2, \n
!! beancan_objects=0, \n
!! / \n
!! &datvtkparameters \n
!! description='absorb', \n
!! filedata=.true., \n
!! transform_type = 'translate', \n
!! end_angle = -1, \n
!! start_position = 0,0,0, \n
!! finish_position = 0,0,2, \n
!! line_divisions = 10 \n
!! / \n
!! &progfiles \n
!! rz_input_file='jetbot1.dat', \n
!! / \n
!! &datvtkparameters \n
!! description='absorb', \n
!! filedata=.true., \n
!! transform_type = 'translate', \n
!! end_angle = 1, \n
!! start_position = 0,0,0, \n
!! finish_position = 0,0,2, \n
!! line_divisions = 10 \n
!! / \n
!! &progfiles \n
!! rz_input_file='jetbot2.dat', \n
!! / \n
!! 
!! Namelist &inputfiles  obvious \n
!! Namelist &miscparameters 
!! for future/expert use etc. \n
!! Namelist &plotselections 
!! almost obvious, see index of outputs in \ref sec1output \n
!! Namelist & geofilparameters 
!! go to doxygen documentation headed "Data Types List", 
!! also accessible directly using
!! firefox or xdg-open ./smardda/smardda-lib/doc/srcdoc/html/classes.html
!! (index under "Data Fields" sub-heading) \n
!! Namelist &datvtkparameters as geofilparameters etc. \n
!! In this example, two absorbing objects are to be defined by
!! translation of lines defined by pairs of points in the files
!! jetbot1.dat and jetbot2.dat. Nonzero end_angle parameter  denotes
!! that the angles at which the absorber appears are determined by the
!! angular extent of the object.
!! 
!! \subsection sec1output SMALIB output files
!!
!! All the codes produce a log file, which contains typically time-and-date stamps,
!! information about the files used, a selection of key parameters, execution times and
!! at end a summary of errors and warnings, terminating with "END OF LOG FILE".
!! Output to the terminal (Fortran print command) consists of the code
!! header, errors associated with log-file writing and a small sub-set of the log-file output,
!! notably serious errors, and the execution times if the code
!! completes normally.
!!
!! The suite is set up so as to terminate codes when a serious error
!! occurs, other errors are classified as warnings, and execution continues.
!! Most also produce a .out file, which is primarily used for inter-code
!! communication, and also acts as a "lock-file", in that the same .ctl
!! file cannot be re-run until the corresponding .out file is deleted.
!! (Good practice is to change the .ctl file-name if the contents change.)
!!
!! Index of suffix strings of output files, with
!! terminating .vtk suffix implies may be visualised using ParaView.
!!
!! The .vtk files begin with the geometry definition, then the quantities in angled brackets.
!! The coordinate system used is indicated by the last letter before
!! the .vtk filetype, viz:
!!     - x for Cartesians (mm)
!!     - m for mapped coordinates (in the current context also Cartesians, possibly in metres)
!!     - q for quantised coordinates, the system in which lines are
!! actually tracked, typically ranging over 0 to 1024.
!!
!! Mapped geometry usually refers to Cartesian coordinates.
!! In the \Bold{smint} single directory,
!! hdsgen and move outputs acquire initial distinguishing "_hds" and "_mov" suffices respectively.
!! \arg _geobj < surface > hdsgen and move output in Cartesian geometry
!! \arg _geobjq < surface > hdsgen and move output in quantised geometry
!! \arg _geopt  points only, hdsgen and move output in Cartesian geometry
!! \arg _geoptq  points only, hdsgen and move output in quantised geometry
!! \arg _hdsm < HDS > hdsgen and move output in mapped geometry, note
!! the move version stores the data less efficiently, both have many points repeated
!! \arg _hdsq < HDS > hdsgen and move output in quantised geometry, note
!! the move version stores the data less efficiently, both have many points repeated
!! \arg _movq  tracks as 3 points move output on quantised geometry
!! \arg _movx tracks as 3 points move output on Cartesian geometry
!! 
!!
!! \subsection sec1resource SMALIB Resources
!! 
!! \arg ./smardda/smardda-lib/doc/srcdoc/html as above, includes links to SMARDDA papers on arXiv,
!! preprint for ray-tracing at http://arxiv.org/abs/1403.6750 ,
!! describing fieldline following at http://arxiv.org/abs/1403.7142 .
!! \arg ./smardda/smardda-lib/Extras \n
!!   Utilities for stl to vtk conversion \n
!! The test-deck should be regarded as a valuable resource, providing
!! examples that can quickly be modified to address new problems.
!! \arg ParaView Guide - a printed version for ParaView 4.3
!! may be purchased via the web-site http://www.paraview.org, or an electronic version
!! of the latest release (5.4 as of January 2019) is a free download.
!! For hints for using ParaView with SMALIB see \ref sec2paraview
!! \arg Repositories based at ITER (IMAS) and CCFE (gitlab) contain codes for use
!! by the fusion community, viz
!! -SMARDDA-PFC, constituent codes geoq and powcal for calculating power deposition
!! on tokamak plasma facing components using simple model of energy flow along 
!! magnetic fieldines 
!! -FLDIFF to model relatively small diffusion of particles normal to fieldlines
!! -SMARDDA-NUCODE for tracking neutrals and reioinised particles in ducts, electrons in
!! bucket sources &c
!! -MC2AT used in paper by Arter and Loughlin,
!! published at http://dx.doi.org/10.1016/j.fusengdes.2008.11.051
!! 
!! \subsection sec1install SMALIB Installation
!! 
!! \arg Clone or pull from IMAS git repository
!! http://git.iter.org/projects/BND/repos/smardda-lib/browse
!! and follow instructions in README, noting that
!! \arg produce executables and documentation by typing in the smalib top directory \n
!! source ./exec/build intel \n
!! gfortran and the debug options  intel_dbg and gfortran_dbg are also available
!! \arg Ensure correct modules are invoked as described in documentation
!! \arg ParaView needs a one-off re-setting of defaults (for
!! better \Bold{smint} integration)
!! \arg SMALIB has a test-deck,
!! see \ref sec3 for more details, which should be used to check the installation.
!!
!! It may be helpful to note also that
!! \arg the commands \n
!! ./exec/decompile and ./exec/compile \n
!! are provided for use after modifying source and/or when changing compilers
!! \arg after the initial build, the compiler may be changed by going to
!! the ./config directory and creating a new config.inc file (link), thus \n
!! cd ./config;ln -sf config_new-compiler.inc config.inc
!! 
!!
!! \section sec2 Description of Constituent Codes and Scripts
!!
!! \ref sec1work  indicates how the codes are combined to perform a
!! ray-tracing calculation.  Each of the codes
!! is briefly described below. As described in the introduction, they are
!! driven ultimately by namelists in .ctl files. The namelists expected in the
!! respective .ctl files are listed below. Precise descriptions of 
!! the usage of the variables in each
!! namelist may be referenced in the online documentation, as data types
!! at ./smardda/smardda-lib/doc/srcdoc/html/classes.html.
!!	
!! \subsection sec2datvtk datvtk code
!! datvtk was designed as a filter, so unlike the other main executables,
!! it accepts options immediately after the name, e.g. \n
!!  datvtk -[opt]  .dat file \n
!!  datvtk -c .ctl file \n
!! In the latter case, the .ctl file
!! contains instructions to generate surfaces by rotating
!! straight lines defined by their end-points, optionally these
!! points may be defined in separate gnu dat format files (\ref sec2formats).
!! Options: \Emph{opt=ABC}, where \n
!! \Emph{A=v} -> .vtk input\n
!! \Emph{B=o} -> orient triangles\n
!! \Emph{B=s} -> shell tetrahedra (BROKEN, use ParaView  instead)\n
!! \Emph{B=f} -> flip triangle orientation\n
!! \Emph{C=d} -> divide into separate block body files\n
!! \Emph{C=m} -> suppress block information \n
!! e.g. datvtk -xom FW34.dat (x as spacer).
!!
!! datvtk can also be as a filter to convert to STL format: \n
!! \Emph{C=b} INERT, output in STL format (divide into separate block body files) \n
!! \Emph{C=t} output as triangles in STL format (suppress block information) \n
!! e.g. datvtk -vxt FW34.vtk (x as spacer).
!!
!! Note that datvtk accepts directly only the .dat Small and Large Field Formats,
!! and not the Free Field Format, as described
!! in the MSC Nastran 2013 Quick Reference Guide. (Unfortunately this
!! document may be obtained only after
!! registering with the MSC Nastran web-site.) There is a further restriction to
!! the small set of "bulk data entries" output by CATIA mainly to define
!! simple meshes. Beware that CATIA is also capable of producing a more
!! compact, non-conforming .dat file output, which SMALIB will not read.
!!
!! Namelist order in .ctl is as follows
!! \arg progfiles 
!! \arg datvtkparameters 
!! \arg plotselections 
!! 
!! \subsection sec2geofil geofil code
!! This code is designed to add special classes of geometry to .vtk
!! files, for example special termination planes, beancans, cutouts, and
!! periodic boundaries, the last so that tracks can interact with repeats
!! of a basic geometry, for example a 20 degree segment of a torus.
!! geofil has a similar structure to the other main executables,
!!
!! Namelist order in .ctl is as follows
!! \arg inputfiles
!! \arg miscparameters
!! \arg plotselections
!! \arg geofilparameters
!!
!! Then following in the same file, there should be namelist(s) corresponding to
!! each of the objects enumerated in geofilparameters. For each object
!! there must appear
!! \arg datvtkparameters
!!
!! and if filedata=.TRUE. in this namelist, there must be a second namelist
!! so for such objects there must appear
!! \arg datvtkparameters
!! \arg progfiles
!! 
!! \subsection sec2hdsgen hdsgen code
!! This computes the multi-octree HDS introduced in \ref sec1alg,
!! which is designed to accelerate the computation of track/ray-triangle intersection.
!! The user should not normally need to be concerned with details of
!! hdsgen operation, except for the need to increase the parameter limit_geobj_in_bin for
!! geometries with a large number of triangles, i.e. greater than
!! approximately 100000 triangles, depending on their degree of clustering.
!!
!! Namelist order in .ctl as follows
!! \arg inputfiles 
!! \arg hdsgenparameters 
!! \arg btreeparameters 
!! \arg positionparameters 
!! \arg plotselections 
!! 
!! \subsection sec2move move code
!! This is the test harness for particle straight-line moves, equivalently rays that
!! uses a hybrid data structure
!! linked to a set of triangles from CAD. 
!! (The HDS is created using hdsgen from the set of
!! triangle coordinates from the CAD model)
!!
!! Namelist order in .ctl as follows
!! \arg inputfiles 
!! \arg numericalparameters 
!! \arg plotselections 
!! 
!! \subsection sec2smint smint script
!! \Bold{smint} is the main GUI for SMALIB, designed to help the beginning user:
!! \arg Assemble geometry from meshed parts or bodies
!! \arg Run and analyse simple problems interactively
!! and start:
!! \arg Constructing scripts for investigations of problem sequences
!! \arg Editing .ctl files to access more advanced featuresa
!!
!! THIS SCRIPT STILL NEEDS A LOT OF WORK TO DO
!! 
!! Primarily \Bold{smint} is a bash 3.0 script using Linux dialog, with features
!! \arg Self-contained as an introduction to SMALIB
!! \arg Point-and-click is optional
!! \arg Help on dialog usage from top-level menu
!! \arg Firefox start-up with links to web-based documentation is also available from within script
!!
!! (bash 3.0 is described in the book Learning the bash shell 3rd Edition, by
!! C.Newham and B.Rosenblatt, O'Reilly 2005)
!! 
!! \Bold{smint} start
!! \arg Create directory and preferably copy geometry
!! files into it (although you can
!! copy later): Linux  \Emph{mkdir},  \Emph{cp} and  \Emph{cd} commands.
!! \arg Attach (\Emph{cd}) to directory, type \Bold{smint}
!! \arg First time, click on Help button \n
!! Note: Ctrl/C (undocumented) can get you out of
!! problems on menus below top menu
!! 
!! \Bold{smint} finish
!! \arg \Bold{smint} as it finishes tells you
!! 1. Where it logged the codes it ran, this file is potential starting-point
!! for a control script (contains !* where errors occurred)
!! 2. Which scratch directory it used (usually ./ismtemp),
!! contains files logging \Bold{smint} script input data for
!! reuse and potentially provenancing.
!! \arg Directory now full of input and output files.
!! .vtk and .png files may be visualised at any time
!! using ParaView, or use eog for png.
!! \arg To visualise results, the recommended approach is to use ParaView to
!! select the geometry, see \ref sec2paraview for details.
!!
!! In the example below, the shadow geometry was assembled from files 
!! partshad1_R_modifier.dat \n
!! partshad1_Z_modifier.dat \n
!! partshad2_R.dat \n
!! partshad2_Z.dat \n
!! beancan_RZ.dat \n
!! results.dat, plus transformations
!! 
!! and the target or results geometry from results.dat.
!! 
!! The beancan is a simple geometrical device to trap any anomalous fieldline
!! integrations within a "beancan", which should lie outside and preferably
!! enclose the volume of interest.
!! 
!! The resulting files in the \Bold{smint} directory were: \n
!! beancan.ctl control file for datvtk execution \n
!! beancan.log file logging output of datvtk \n
!! beancan_RZ.dat point data used to define a surface by datvtk \n
!! beancan.vtk ASCII file defining geometry as set of oriented triangles \n
!! cmdwrap\Emph{nnn}.log start-point for script construction \n
!! ismfiles.sav records principal geometry file(s) used by smint \n
!! ismform.sav records principal parameters used by smint \n
!! ismtemp scratch directory \n
!! partshad1_modifier.ctl control file for vtktfm execution \n
!! partshad1_modifier.log file logging output of vtktfm \n
!! partshad1_modifier.vtk ASCII file defining geometry as set of oriented triangles \n
!! partshad1_R_modifier.dat point data used to define a surface by datvtk \n
!! partshad1_Z_modifier.dat point data used to define a surface by datvtk \n
!! partshad2.ctl control file for datvtk execution \n
!! partshad2.log file logging output of datvtk \n
!! partshad2_R.dat point data used to define a surface by datvtk \n
!! partshad2.vtk ASCII file defining geometry as set of oriented triangles \n
!! partshad2_Z.dat point data used to define a surface by datvtk \n
!! results.dat mesh data converted to .vtk surface format by datvtk \n
!! results.ctl control file for move execution \n
!! results.log file logging output of move \n
!! results_move.out output file from move \n
!! results_mov.vtk output of tracks in Cartesian geometry \n
!! results_movq.vtk output of tracks in quantised geometry \n
!! results_tfm1.ctl control file for vtktfm execution \n
!! results_tfm1.ctlin control file for ctlgen execution \n
!! results_tfm1.log file logging output of vtktfm \n
!! results_tfm1.vtk vtktfm output of geometry after transformation in Cartesians \n
!! results_tfm2.ctl control file for vtktfm execution \n
!! results_tfm2.ctlin control file for ctlgen execution \n
!! results_tfm2.log file logging output of vtktfm \n
!! results_tfm2.vtk vtktfm output of geometry after transformation in Cartesians \n
!! results.vtk ASCII file defining geometry as set of oriented triangles \n
!! shadow.ctl control file for geofil execution \n
!! shadow.ctlin control file for ctlgen execution \n
!! shadow.vtk geofil output of geometry \n
!! shadow_geofil.out output file from geofil \n
!! shadow_gnu_cart.gnu file for input to gnuplot with all 3 Cartesian components of B-field \n
!! shadow_hds.ctl control file for hdsgen execution \n
!! shadow_hds_geobjq.vtk hdsgen output of geometry in quantised space \n
!! shadow_hds_geoptq.vtk hdsgen output of geometry as points only  in quantised space \n
!! shadow_hds_hdsgen.out output file from hdsgen \n
!! shadow_hds_hds.hds file defining HDS produced by hdsgen \n
!! shadow_hds_hdsm.vtk hdsgen geometric output of HDS in mapped coordinates \n
!! shadow_hds_hdsq.vtk hdsgen geometric output of HDS in quantised space \n
!! shadow_hds.log file logging output of hdsgen \n
!! shadow.log file logging output of geofil \n
!! shadow.vtk ASCII file defining geometry as set of oriented triangles \n
!! 
!!
!! \subsection sec2vtktfm vtktfm code
!! The code transforms and/or combines .vtk format surface geometry files (ASCII format).
!! In the \Bold{smint} script, vtktfm is used either to apply a solid body
!! transformation or to combine files, but in fact 
!! the code can simultaneously transform and combine, although no more
!! than one transformation on any body in one code usage.
!! It can also or instead  extract geometry according to angle constraints,
!! useful for cutting down geometry to  a periodic cell.
!!
!! Namelist order in .ctl as follows
!! \arg miscparameters 
!! \arg vtktfmparameters (only if new_controls=.TRUE. in miscparameters)
!! \arg vtkfiles 
!! \arg panelarrayparameters 
!! \arg positionparameters 
!! 
!! 
!! \subsection sec2paraview ParaView and SMALIB
!! The initial chapters of the official guide explain basic concepts such as the 
!! graphics pipeline and should be at least skim-read before reading these
!! SMARDDA-LIB-specific hints. (As of ParaView 5.2, there is an initial pop-up window
!! with information for the beginner.) This preliminary should help the user  acquire
!! the necessary familiarity with basic ParaView GUI 
!! layout, knowledge of when to use the "Apply" button, and thereby avoid
!! the commonest mistakes and sources of confusion.
!!
!! General points worth noting are that
!! \arg ParaView is continually being developed, and the hints here assume
!! version 4.4 or later.
!! \arg SMARDDA-LIB in normal use only generates surface information, so many ParaView
!! functions which require volume data are not available or needed.
!! \arg ParaView starts counting from zero,
!! whereas SMARDDA-LIB is more conventional, so that for example, the first track
!! starts from triangle zero.
!! \arg Toggling the gear-wheel symbol gives access to a wider and often useful
!! range of extra functionality.
!! 
!! The functionality most often needed for SMARDDA-LIB concerns the making and saving
!! of selections of the data in the .vtk files.
!! The selection process is activated by toggling the "Select Cells on" icon
!! immediately above the main "render" window ("Select Cells With Polygon"
!! offers finer control). A selection box is defined by the mouse and the
!! selected cells (triangles) are highlighted. The selection may be inspected
!! in more detail, by use of the Selection Display Inspector
!! (View -> Selection Display Inspector). If satisfactory, the "Extract
!! Selection" icon which is on the menu bar just above the Select Cells icons,
!! should be clicked upon. File->Save data can then be used
!! to save the selection to a Legacy vtk (.vtk) format file, using the 
!! options to save all timesteps to the Ascii file type. The new .vtk file
!! thereby produced can be used as an input to a new SMARDDA-LIB run.
!! datvtk preserves block information from each .dat file when it
!! is converted to .vtk (unfortunately it is discarded in subsequent
!! processing). If block numbers are known, or identified from
!! ParaView visulation, their geometry can be selected by first
!! using Edit -> "Find data", where "is one of" is the most flexible
!! selection option. Having created the selection use the "Extract 
!! Selection" icon as above, clicking on "Copy active selection".
!! 
!! Most other functionality is accessed by the Filters menu. Of the myriad
!! options, the ones most relevant to SMARDDA-LIB are
!! \arg Integrate variables, for total area.
!! \arg Calculator, operates on velocity fields to give V.n etc.
!! \arg Point data -> Cell data, and vice versa, since Calculator only operates
!! on variables of the same data type.
!! \arg Cell centres, when producing spread-sheet .csv output.
!! \arg Generate Ids, to help locate elements and points on the surface.
!! \arg Append attributes, enables data from different .vtk files with same
!! geometry to be combined.
!! \arg Plot on intersection curves, for line plots of data as functions
!! of arclength on surface, see "Line plots" below.
!! 
!! Views in the render window  may be saved to .png file for plotting
!! using File->Save Screenshot.
!! The same view may be applied to a separate session by using File->Save State,
!! then File->Load State in the new ParaView session (same release required).
!!
!! Spread-sheets may be produced using File->Save Data, for which .csv is
!! the default option. For point data this is straighforward, as the
!! position is saved with the data (consistent with ParaView usage, the 
!! first vector component is zero). Cell centre points need to be
!! explicitly defined by the user with the Cell centre filter before using
!! File->Save Data with cell data such as Q.
!!
!! Line plots using ParaView
!!
!! These may require interpolation from cell-centred values to
!! point values, which takes place in a manner than does not seem to
!! documented officially. Checks indicate that the point value is an average over
!! cell-centred values for all triangles meeting at the point (probably
!! area-weighted), hence edges get blurred where typically only 2 cells
!! out of 6 or 8 have non-zero field values.  The intersection curve is likely to
!! be generated by the intersection of a plane and the surface, and
!! when generically the intersection between the plotting plane does not
!! align with the surface triangulation points, there is a further interpolation
!! to points on the intersection. The intersection curve is not sampled uniformly
!! but at points which seem to correspond to the projections of the vertices
!! of intersected triangles. \n
!! Hints for line plots\n
!! It is often easier to reduce the size of the geometry in the window
!! by using  "select cells through" icon just above followed by  \n
!! Extract selection, (on bar just above) then \n
!! Cell data-> point data \n
!! The filter is "Plot on intersection curves", when
!! the menu brought up is comprehensive if the gear-wheel is clicked.
!! However it will often be found helpful to save the data to a .csv file
!! for graphing using gnuplot (comment out first line with leading "#", replace
!! commas with spaces) or spreadsheet software.
!!
!!
!! \subsection sec2formats Format of gnu dat and .vtk file headers
!!
!! Special header lines have been introduced for gnu dat and .vtk files 
!! to encode information about the geometrical data contained, as to
!! which coordinate system and units are used therein. In a few cases, it
!! may help the user to know the format of these headers, hence a
!! description follows. There are examples in the test deck
!! 
!! For gnu dat files, the header should begin with '#', then \f$29\f$ characters
!! of free format, non-blank description terminated by space. After the
!! space is the string 'Number_Parameters=', followed by the number of
!! integer parameters in the rest of the line. Currently, this number is
!! six in Fortran I3 integer format, followed by numbers in 1X,I4 format, meanings as follows
!! 1. default 1
!! 2. descriptive code for geometry, see geobj for full list, 
!! eg. 0. geometry, 4. skylight, 5. terminating geometry, 6. beancan, 7. cutout.
!! 3. units (consistent with posang, -3 for mm, 0 for m)
!! 4. coordinate system, after posang, namely \n
!!   0 - position in Cartesians \f$ (x,y,z) \f$ \n
!!   1 - position in cylindrical polars \f$ (R,Z,\zeta), |\zeta| \leq \pi \f$ \n
!!   2 - position in flux coordinates  \f$ (\psi,\theta,\zeta) \f$ \n
!! 5. default 0
!! 6. default 0
!!
!! For .vtk files, the header is actually the second line, starting with
!! \f$30\f$ characters of free format, non-blank description terminated by space. After the
!! space is the string 'Number_Parameters=', followed by the number of
!! integer parameters in the rest of the line. Currently, this number is
!! six in Fortran I3 integer format, followed by numbers in 1X,I4 format, meanings as follows
!! 1. infilelevel (level of refinement in file), default 1
!! 2. indicates triangles geometry coded as to type
!! 3. units (consistent with posang, -3 for mm, 0 for m)
!! 4. coordinate system, after posang, namely \n
!!   0 - position in Cartesians \f$ (x,y,z) \f$ \n
!!   1 - position in cylindrical polars \f$ (R,Z,\zeta), |\zeta| \leq \pi \f$ \n
!!   2 - position in flux coordinates  \f$ (\psi,\theta,\zeta) \f$ \n
!! 5. quantised (if 1, else 0)
!! 6. (number of transforms applied) - (number of inverse transforms applied)
!!
!!
!! \section sec3 Test-Deck
!! 
!! \subsection sec3intro Introduction
!! 
!! The test-deck, including the validation test cases,  is designed to
!! demonstrate correct functioning of all the commonly used parts of
!! SMARDDA-LIB. Data for the test-deck is in the TEST sub-directory,
!! together with files "_notes.txt" on how to produce geometry for the validation
!! test cases.
!!
!! \subsection sec3cases Test Cases
!! A brief description of each case follows so that the
!! user may determine which is closest to the problem of current
!! interest and therefore might be most easily modified to solve it. 
!! Test cases from the SMARDDA-LIB creation are performed
!! using data in subdirectories of TEST.
!!
!! As more cases are added,
!! the convention will be to create new directories
!! of which the names start with "Test-". The remaining components of 
!! the names in order will be 
!! an identifier of the shadowing geometry and
!! the source geometry respectively. 
!! Thus the cases may be easily tied to the geometries among the .vtk files
!! in ./smardda/smardda-lib/TEST and visualised using ParaView, although note that
!! some .vtk files  have to be generated from the
!! .dat files.
!!
!! \subsubsection sec3T1 GFRUN
!! \image html absb2.png
!! \image latex absb2.eps "Augmented geometry" width=10cm
!! This test-case is to add vertical bounding planes to a simplified
!! JET divertor model confined to one octant. The original divertor
!! geometry corresponds to Body 1 in the figure and the bounding planes
!! to Body 2. The planes are described as absorbers, so that physically
!! they behave like the rest of the geometry.
!! 
!! To execute the test case, copy $SMALIB_DIR/TEST/GFRUN to a work
!! directory, removing the comparison output files, run geofil, then 
!! compare results, ie: \n
!! cp -r $SMALIB_DIR/TEST/GFRUN workdir1 \n
!! cd workdir1 \n
!! rm -f absb_geofil.out absb.vtk absb.log \n
!! geofil absb \n
!! diff -b absb.vtk $SMALIB_DIR/TEST/GFRUN/absb.vtk
!!
!! The figure should be reproducible using ParaView
!!
!! \subsubsection sec3T2 MTRUN
!! \image html hds.png
!! \image latex hds.eps "HDS generation for ITER duct geometry" width=10cm
!! \image html move.png
!! \image latex move.eps "Ray intersections with ITER duct geometry" width=10cm
!! This test-case geometry is used in the first major SMARDDA publication,
!! preprint at http://arxiv.org/abs/1403.6750 ,
!! where the concept of hierarchical data structure (HDS) is
!! discussed in the ray-tracing context.
!! 
!! The tracks in file exps00000200.qry may be produced as described in the README file in
!! directory $SMALIB_DIR/TEST/RAYS .
!! A file with this name already exists in directory
!! $SMALIB_DIR/TEST/RAYS and is used by program move in the test case
!! below.
!!
!! To execute the test case, copy $SMALIB_DIR/TEST/MTRUN to a work
!! directory, removing the comparison output files, run hdsgen and move, then 
!! compare results, ie: \n
!! cp -r $SMALIB_DIR/TEST/MTRUN workdir2 \n
!! cd workdir2 \n
!! rm -f *_* *.log \n
!! hdsgen testhds \n
!! diff -b testhds_hds.hds $SMALIB_DIR/TEST/MTRUN/testhds_hds.hds \n
!! move test \n
!! diff -b test_movx.vtk $SMALIB_DIR/TEST/MTRUN/test_movx.vtk
!!
!! The figures should be reproducible using ParaView
!!
