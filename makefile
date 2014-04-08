#MAKEFILE for the BEMRI project

#Variable Declarations
FILES=cosmo.c file.c kb.c load_params.c main.c matrix.c BEMRI_functions.c kepler_eq.c random.c quadrupole.c ODE.c newtonian_n_body.c chain.c leapfrog.c args.c
CFLAGS=-c -std=c99 
SRCDIR=./code/
OBJDIR=./obj/
SOURCES=$(FILES:%.c=$(SRCDIR)%.o)
OBJECTS=$(FILES:%.c=$(OBJDIR)%.o)

#all and clean targets
all: BEMRI

#make all and run
run: BEMRI
	 ./BEMRI

clean:
	rm -rf $(OBJECTS) BEMRI

#Executable target
BEMRI: make_dir $(OBJECTS)
	gcc -lm $(OBJECTS) -o BEMRI

#Object File targets
#note: $@ is replaced by the name of the target
# and $< is the name of the first dependency
$(OBJDIR)cosmo.o: $(SRCDIR)cosmo.c
	gcc -lm $(CFLAGS) $< -o $@
	
$(OBJDIR)file.o: $(SRCDIR)file.c
	gcc $(CFLAGS) $< -o $@
	
$(OBJDIR)kb.o: $(SRCDIR)kb.c
	gcc $(CFLAGS) $< -o $@
	
$(OBJDIR)load_params.o: $(SRCDIR)load_params.c
	gcc $(CFLAGS) $< -o $@

$(OBJDIR)main.o: $(SRCDIR)main.c
	gcc $(CFLAGS) $< -o $@

$(OBJDIR)matrix.o: $(SRCDIR)matrix.c
	gcc $(CFLAGS) $< -o $@
	
$(OBJDIR)BEMRI_functions.o: $(SRCDIR)BEMRI_functions.c
	gcc $(CFLAGS) $< -o $@
	
$(OBJDIR)newtonian_n_body.o: $(SRCDIR)newtonian_n_body.c
	gcc $(CFLAGS) $< -o $@	
	
$(OBJDIR)kepler_eq.o: $(SRCDIR)kepler_eq.c
	gcc $(CFLAGS) $< -o $@		

$(OBJDIR)random.o: $(SRCDIR)random.c
	gcc $(CFLAGS) $< -o $@	

$(OBJDIR)quadrupole.o: $(SRCDIR)quadrupole.c
	gcc $(CFLAGS) $< -o $@	

$(OBJDIR)ODE.o: $(SRCDIR)ODE.c
	gcc $(CFLAGS) $< -o $@	

$(OBJDIR)chain.o: $(SRCDIR)chain.c
	gcc $(CFLAGS) $< -o $@	

$(OBJDIR)leapfrog.o: $(SRCDIR)leapfrog.c
	gcc $(CFLAGS) $< -o $@	

$(OBJDIR)args.o: $(SRCDIR)args.c
	gcc $(CFLAGS) $< -o $@	
		
make_dir:
	if test -d ./obj; then echo; else mkdir ./obj; fi		
	
