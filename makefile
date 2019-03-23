# YouJunHu's makefile
PROGRAM_NAME:=GTAW
COMPILER:=	gfortran
BUILDER:=	$(COMPILER)
#OPTION:= -fbounds-check -fopenmp
OPTION:= -fbounds-check 
COMPILE:=	$(COMPILER) $(OPTION) -c 
BUILD:=		$(BUILDER) $(OPTION) -o $(PROGRAM_NAME) 
#lapack_location:=/home/yj/installed/lapack-3.1.1-gfortran/
lapack_location:=/usr/lib/lapack/liblapack.so.3
blas_location:=/usr/lib/libblas/libblas.so.3


f90sources:= modules.f90  main.f90 read_gfile.f90 arrange_lcfs.f90 calculate_contours.f90 poloidal_magnetic_flux.f90 calculate_poloidal_flux_partial_derivatives.f90  calculate_fourier_integral.f90 calculate_weight_functions.f90 calculate_matrix_elements.f90 calculate_continous_spectrum.f90 local_magnetic_shear.f90 mnewt.f90 calculate_jacobian.f90 calculate_bsq.f90 curvature.f90 filter_kappas.f90 arc_length.f90 normalized_parallel_current_density.f90 cylindrical_continua.f90 math.f90 calculate_radial_matrix.f90 invert_matrix.f90 shoot.f90 set_density.f90 find_global_mode.f90 calculate_kernel.f90 interpolate_flux_function.f90 calculate_phase.f90 poloidal_displacement_and_divergence_of_displacement.f90 midplane.f90  major_radius_on_lfs_midplane.f90 interpolate.f90 find_q_equal_one_surface.f90 contour.f90 diagnostic.f90 set_radial_range.f90 choose_radial_coordinate.f90 fast_ions_pressure.f90 artifial_mode.f90
f77sources:=splie2.for splin2.for  spline.for  splint.for rtbis.for dftcor.for  dftint.for  four1.for  polint.for  realft.for newt_complex.for  fmin.for lnsrch.for lubksb.for ludcmp.for broydn.for qrdcmp.for qrupdt.for rsolv.for rotate.for 

f90objs:= $(f90sources:.f90=.o)
f77objs:= $(f77sources:.for=.o)


$(PROGRAM_NAME):	$(f90objs) $(f77objs)
	$(BUILD) $(f90objs) $(f77objs) $(lapack_location) $(blas_location)


modules.o: modules.f90
	$(COMPILE) $< -o $@
artifial_mode.o: artifial_mode.f90 modules.f90
	$(COMPILE) $< -o $@
choose_radial_coordinate.o: choose_radial_coordinate.f90 modules.f90
	$(COMPILE) $< -o $@
fast_ions_pressure.o: fast_ions_pressure.f90 modules.f90
	$(COMPILE) $< -o $@
set_radial_range.o: set_radial_range.f90 modules.f90
	$(COMPILE) $< -o $@
poloidal_magnetic_flux.o: poloidal_magnetic_flux.f90
	$(COMPILE) $< -o $@
contour.o: contour.f90 modules.f90
	 $(COMPILE) $< -o $@
diagnostic.o: diagnostic.f90 modules.f90
	$(COMPILE) $< -o $@
find_q_equal_one_surface.o: find_q_equal_one_surface.f90 modules.f90
	 $(COMPILE) $< -o $@
main.o: main.f90 modules.f90
	$(COMPILE) $< -o $@
filter_kappas.o: filter_kappas.f90 modules.f90
	$(COMPILE) $< -o $@
arrange_lcfs.o: arrange_lcfs.f90 modules.f90
	$(COMPILE) $< -o $@
midplane.o: midplane.f90 modules.f90
	$(COMPILE) $< -o $@
calculate_phase.o: calculate_phase.f90 modules.f90
	$(COMPILE) $< -o $@
find_global_mode.o: find_global_mode.f90 modules.f90
	$(COMPILE) $< -o $@
poloidal_displacement_and_divergence_of_displacement.o: poloidal_displacement_and_divergence_of_displacement.f90 modules.f90
	$(COMPILE) $< -o $@
calculate_kernel.o: calculate_kernel.f90 modules.f90
	$(COMPILE) $< -o $@
interpolate_flux_function.o: interpolate_flux_function.f90 modules.f90
	$(COMPILE) $< -o $@
read_gfile.o: read_gfile.f90 modules.f90
	$(COMPILE) $< -o $@
set_density.o: set_density.f90 modules.f90
	$(COMPILE) $< -o $@
cylindrical_continua.o: cylindrical_continua.f90 modules.f90
	$(COMPILE) $< -o $@
math.o: math.f90 modules.f90
	$(COMPILE) $< -o $@

calculate_radial_matrix.o: calculate_radial_matrix.f90 modules.f90
	$(COMPILE) $< -o $@

shoot.o: shoot.f90 modules.f90
	$(COMPILE) $< -o $@
calculate_poloidal_flux_partial_derivatives.o: calculate_poloidal_flux_partial_derivatives.f90 modules.f90
	$(COMPILE) $< -o $@
curvature.o: curvature.f90 modules.f90
	$(COMPILE) $< -o $@
calculate_bsq.o: calculate_bsq.f90 modules.f90
	$(COMPILE) $< -o $@
local_magnetic_shear.o: local_magnetic_shear.f90 modules.f90
	$(COMPILE) $< -o $@
normalized_parallel_current_density.o: normalized_parallel_current_density.f90 modules.f90
	$(COMPILE) $< -o $@
calculate_contours.o: calculate_contours.f90 modules.f90
	$(COMPILE) $< -o $@
arc_length.o: arc_length.f90 modules.f90
	$(COMPILE) $< -o $@
calculate_weight_functions.o: calculate_weight_functions.f90 modules.f90
	$(COMPILE) $< -o $@

calculate_fourier_integral.o: calculate_fourier_integral.f90 modules.f90
	$(COMPILE) $< -o $@
calculate_matrix_elements.o: calculate_matrix_elements.f90 modules.f90
	$(COMPILE) $< -o $@
calculate_continous_spectrum.o: calculate_continous_spectrum.f90 modules.f90
	$(COMPILE) $< -o $@
calculate_jacobian.o: calculate_jacobian.f90 modules.f90
	$(COMPILE) $< -o $@
major_radius_on_lfs_midplane.o:major_radius_on_lfs_midplane.f90 modules.f90
	$(COMPILE) $< -o $@
interpolate.o: interpolate.f90
	$(COMPILE) $< -o $@
mnewt.o: mnewt.f90 modules.f90
	$(COMPILE) $< -o $@
broydn.o: broydn.for modules.f90
	$(COMPILE) $< -o $@

qrdcmp.o: qrdcmp.for modules.f90
	$(COMPILE) $< -o $@
qrupdt.o: qrupdt.for modules.f90
	$(COMPILE) $< -o $@
rsolv.o: rsolv.for modules.f90
	$(COMPILE) $< -o $@
rotate.o: rotate.for modules.f90
	$(COMPILE) $< -o $@
splie2.o: splie2.for modules.f90
	$(COMPILE) $< -o $@
invert_matrix.o: invert_matrix.f90 modules.f90
	$(COMPILE) $< -o $@
splin2.o: splin2.for modules.f90
	$(COMPILE) $< -o $@

spline.o: spline.for modules.f90
	$(COMPILE) $< -o $@
splint.o: splint.for modules.f90
	$(COMPILE) $< -o $@
rtbis.o: rtbis.for modules.f90
	$(COMPILE) $< -o $@

dftcor.o: dftcor.for modules.f90
	$(COMPILE) $< -o $@

dftint.o: dftint.for modules.f90
	$(COMPILE) $< -o $@
four1.o: four1.for modules.f90
	$(COMPILE) $< -o $@

polint.o: polint.for modules.f90
	$(COMPILE) $< -o $@

realft.o: realft.for modules.f90
	$(COMPILE) $< -o $@
newt_complex.o: newt_complex.for modules.f90
	$(COMPILE) $< -o $@
fmin.o: fmin.for modules.f90
	$(COMPILE) $< -o $@
lnsrch.o: lnsrch.for modules.f90
	$(COMPILE) $< -o $@
lubksb.o: lubksb.for modules.f90
	$(COMPILE) $< -o $@
ludcmp.o: ludcmp.for modules.f90
	$(COMPILE) $< -o $@

.PHONY : run clean tarfile merge publish

run: $(PROGRAM_NAME)
	./$(PROGRAM_NAME)
clean :
	 rm -f $(PROGRAM_NAME)  $(f90objs) $(f77objs) *.mod  
tarfile:
	mkdir $(PROGRAM_NAME)_version`date +"%Y-%m-%d"` && cp -r $(f90sources) $(f77sources) makefile gtaw.in gfile $(PROGRAM_NAME)_version`date +"%Y-%m-%d"` && tar -cf $(PROGRAM_NAME)_version`date +"%Y-%m-%d"`.tar $(PROGRAM_NAME)_version`date +"%Y-%m-%d"` && rm -r $(PROGRAM_NAME)_version`date +"%Y-%m-%d"` 

merge:
	cat $(f90sources) $(f77sources) > source_merged.f90

publish: tarfile
	rsync -avz $(PROGRAM_NAME)_version`date +"%Y-%m-%d"`.tar yj@shenma.ipp.ac.cn:/home/yj/public_html/codes/GTAW_src.tar
sededit:
	for i in $(f90sources);do sed -i.yjyjt 's/poloidal_flux_normalized/pfn/gI' $$i;done
	 #for i in $(f90sources);do sed 's/poloidal_flux_normalized/pfn/gI' $$i>/tmp/GTAW/$$i;done;cp $(f77sources) makefile gtaw.in g0*.* gfile_solovev electron_number_density_g038300.03900.dat electron_number_density_uniform.dat /tmp/GTAW
