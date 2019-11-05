################################### BEWARE OF SPACE INSIDE $()!!!! 
# Macros
EXE = thermocline
F95 = gfortran
#FLAGS = -fdefault-real-8 -g -fimplicit-none -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid -O0
FLAGS = -fdefault-real-8 -O2

# Source and object
SRCF95 = InOut.f95 SemiDiscret.f95 CodeVerif.f95 main.f95

# Pattern rules
%.o: %.f95
	$(F95) $(FLAGS) -c $< -o $@

OBJF95 = $(SRCF95:.f95=.o)

# Target
$(EXE): $(OBJF95)
	$(F95) $(OBJF95) -o $@

.PHONY: run
run: $(EXE) plot.py
	./$(EXE)
	python3 plot.py

clean:
	rm $(OBJF95) *.mod $(EXE)

# Dependencies
main.o: InOut.o SemiDiscret.o CodeVerif.o
CodeVerif.o: InOut.o SemiDiscret.o

