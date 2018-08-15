test:
	gfortran -fno-second-underscore -DBL_SPACEDIM=2 -DBL_FORT_USE_UNDERSCORE -g tetrahedron_keast_rule.F90 -c -o tetrahedron_keast_rule.o
	gfortran -fno-second-underscore -DBL_SPACEDIM=2 -DBL_FORT_USE_UNDERSCORE -g triangle_fekete_rule.F90 -c -o triangle_fekete_rule.o  
	gfortran -fno-second-underscore -DBL_SPACEDIM=2 -DBL_FORT_USE_UNDERSCORE -g PROBCOMMON.F90 -c -o PROBCOMMON.o  
	gfortran -fno-second-underscore -DBL_SPACEDIM=2 -DBL_FORT_USE_UNDERSCORE -g GLOBALUTIL.F90 -c -o GLOBALUTIL.o  
	gfortran -fno-second-underscore -DBL_SPACEDIM=2 -DBL_FORT_USE_UNDERSCORE -g MOF.F90 -c -o MOF.o  
	gfortran -fno-second-underscore -DBL_SPACEDIM=2 -DBL_FORT_USE_UNDERSCORE -g vof_cisl.F90 -c -o vof_cisl.o
	gfortran -fno-second-underscore -DBL_SPACEDIM=2 -DBL_FORT_USE_UNDERSCORE -g vfrac_pair.F90 -c -o vfrac_pair.o
	gfortran -fno-second-underscore -DBL_SPACEDIM=2 -DBL_FORT_USE_UNDERSCORE -g multimat_FVM.F90 -c -o multimat_FVM.o
	gfortran -fno-second-underscore -DBL_SPACEDIM=2 -DBL_FORT_USE_UNDERSCORE -g BICGSTAB_Yang_MULTI.F90 -c -o BICGSTAB_Yang_MULTI.o 
	gfortran -fno-second-underscore -DBL_SPACEDIM=2 -DBL_FORT_USE_UNDERSCORE -g interp_code.F90 -c -o interp_code.o 
	gfortran -fno-second-underscore -DBL_SPACEDIM=2 -DBL_FORT_USE_UNDERSCORE -g main.F90 -c -o main.o 
	gfortran -g -o test tetrahedron_keast_rule.o triangle_fekete_rule.o PROBCOMMON.o GLOBALUTIL.o MOF.o BICGSTAB_Yang_MULTI.o vof_cisl.o vfrac_pair.o multimat_FVM.o interp_code.o main.o

run:
	./test
clean:
	rm *.mod core.* test *.jpg *.o *.tec wavedata*
clear:	
	rm *.dat
