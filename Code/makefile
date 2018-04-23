CC=gcc
#-----------------CFLAGS-------------------------------------------------------
OBJS_geom=./lib/my_geom.o ./lib/my_memory.o ./lib/messages.o ./lib/geom_prop.o ./lib/histogram.o ./lib/quaternions.o
#OBJS_psi=./psi_move/theta1.o ./psi_move/theta3.o ./psi_move/theta5.o ./psi_move/theta7.o ./psi_move/theta2.o ./psi_move/theta4.o ./psi_move/theta6.o ./psi_move/T7.o ./psi_move/robo_precomputed.o 
OBJSCAT=./lib/Caterpillar.o ./lib/Caterpillar_energies.o ./lib/my_geom.o ./lib/CAT_moves.o ./generic_move/CR_DerT7.o  ./generic_move/CR_JacT7.o  ./generic_move/CR_precomp.o ./generic_move/CR_T7.o
OBJSOUT=./lib/ABM_pdbio.o ./lib/Caterpillar_IO.o
#CFLAGS=-Wall -g -std=c99 -O0 -fstack-protector
CFLAGS=-Wall -g -std=c99 -O3 -p -fno-omit-frame-pointer -flto -fstrict-aliasing -DCBETA -Werror=implicit-function-declaration -Werror=incompatible-pointer-types -funsafe-loop-optimizations -march=native
#CFLAGS=-Wall -g -std=c99 -p -fno-omit-frame-pointer -flto -fstrict-aliasing -DCBETA -Werror=implicit-function-declaration -Werror=incompatible-pointer-types -funsafe-loop-optimizations -march=native
#CFLAGS=-static -Wall -g -std=c99 -p -fno-omit-frame-pointer -flto -fstrict-aliasing -I/home/mirekj/gsl2.4/include -L/home/mirekj/gsl2.4/lib -DCBETA -Werror=implicit-function-declaration -Werror=incompatible-pointer-types -funsafe-loop-optimizations -march=native
#CFLAGS=-static -Wall -g -std=c99 -p -fno-omit-frame-pointer -flto -fstrict-aliasing -I/software/gsl-1.16/gcc/include/ -L/software/gsl-1.16/gcc/lib -DCBETA -Werror=implicit-function-declaration -Werror=incompatible-pointer-types -funsafe-loop-optimizations -march=native

minimal.x: minimal.c ./lib/messages.o ./generic_move/CR_DerT7.o ./generic_move/CR_JacT7.o ./generic_move/CR_precomp.o ./generic_move/CR_T7.o
	$(CC) $(CFLAGS)  $^ -lm -lgsl  -lgslcblas -o $@

concerted_rot.x: concerted_rot.c $(OBJS_geom) $(OBJSCAT) $(OBJSOUT) 
	$(CC) $(CFLAGS)  $^ -lm -lgsl  -lgslcblas -o $@

test_CR_speed.x: ./tests/CR_speed/test_CR_speed.c $(OBJS_geom) $(OBJSCAT) $(OBJSOUT)
	$(CC) $(CFLAGS)  $^ -lm -lgsl  -lgslcblas -o $@

test_quaternions.x: ./tests/quaternions/quaternions.c $(OBJS_geom) $(OBJSCAT) $(OBJSOUT)
	$(CC) $(CFLAGS) $^ -lm -lgsl -lgslcblas -o $@

test_phi_psi_flat_histo.x: ./tests/phi_psi_flat_histo/phi_psi_flat_histo.c $(OBJS_geom) $(OBJSCAT) $(OBJSOUT)
	$(CC) $(CFLAGS) $^ -lm -lgsl -lgslcblas -o $@

test_phi_psi_excludedVolume.x: ./tests/phi_psi_excludedVolume/phi_psi_excludedVolume.c $(OBJS_geom) $(OBJSCAT) $(OBJSOUT)
	$(CC) $(CFLAGS) $^ -lm -lgsl -lgslcblas -o $@

test_phi_psi_excludedVolumeInteraction.x: ./tests/phi_psi_excludedVolumeInteraction/phi_psi_excludedVolumeInteraction.c $(OBJS_geom) $(OBJSCAT) $(OBJSOUT)
	$(CC) $(CFLAGS) $^ -lm -lgsl -lgslcblas -o $@


./lib/%.o:./lib/%.c
	$(CC) $(CFLAGS) $< -o $@ -c

./tests/CR_speed/%.o:./tests/CR_speed/%.c
	$(CC) $(CFLAGS) $< -o $@ -c

./tests/quaternions/%.o:./tests/quaternions/%.c
	$(CC) $(CFLAGS) $< -o $@ -c

./generic_move/%.o:./lib/%.c
	$(CC) $(CFLAGS) $< -o $@ -c
