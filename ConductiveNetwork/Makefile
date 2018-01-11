#  Declaration
CXXFLAGS = -Wall -v -O3 -std=c++0x
PARFLAGS = -fopenmp
CC	=	gcc
C++	=	g++
F77    	=	f77

.SUFFIXES:

.SUFFIXES: .f .F .cpp .c .o .h

%.o: %.f
	${F77} ${FFLAGS} $< -c
	
%.o: %.F
	${F77} ${FFLAGS} $< -c

%.o: %.cpp %.h
	${C++} ${CXXFLAGS} $< -c

%.o: %.c %.h
	${CC} ${CFLAGS} $< -c


objects =  App_Network_3D.o Backbone_Network.o Background_vectors.o Clusters_fractions.o Contact_grid.o Cutoff_Wins.o Direct_Electrifying.o \
           Electrical_analysis.o Fem_3D.o GenNetwork.o Geometry_3D.o Hns.o Hoshen_Kopelman.o Input_Reader.o \
           MathMatrix.o Percolation.o Printer.o Tecplot_Export.o Triangulation.o MainPro.o \
	                   
necn : $(objects)        
	${C++} ${PARFLAGS} -o necn $(objects)
	              
MainPro.o : MainPro.cpp
	${C++} ${CXXFLAGS} MainPro.cpp -c

# clean
.PHONY : clean
clean :
	-rm necn $(objects)
