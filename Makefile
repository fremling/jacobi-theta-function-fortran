
.PHONY.:test
test:run_test_suite_jacobitheta
	./run_test_suite_jacobitheta

run_test_suite_jacobitheta: test_suite_jacobitheta.o
	gfortran test_suite_jacobitheta.o jacobitheta.o typedef.o test_utilities.o dedekind_eta_function.o misc_random.o -o run_test_suite_jacobitheta

%.o:%.f90
	gfortran -c  $<


test_suite_jacobitheta.o:jacobitheta.o typedef.o test_utilities.o dedekind_eta_function.o misc_random.o
test_utilities.o:typedef.o
dedekind_eta_function.o:typedef.o


