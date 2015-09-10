#  Declaration
CXXFLAGS = -Wall -v -O3 
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


objects	=  App_Network_3D.o Backbone_Network.o Cutoff_Wins.o Direct_Electrifying.o \
	   Fem_3D.o GenNetwork.o Geometry_3D.o Hns.o Input_Reader.o \
	   Hoshen_Kopelman.o MathMatrix.o Tecplot_Export.o MainPro.o \
	                   
necn : $(objects)        
	${C++} ${PARFLAGS} -o necn $(objects)
	              
MainPro.o : MainPro.cpp
	${C++} ${CXXFLAGS} MainPro.cpp -c

# 伪目标文件
.PHONY : clean
clean :
	-rm necn $(objects)
