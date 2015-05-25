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


objects	= Algorithm_FEM.o Fem_3D.o Gauss.o Geometry_3D.o \
	  GeoNano.o Global_Stiff_Matrix.o Hns.o \
	  MatBase.o MathMatrix.o MatPro.o Mesher.o Neca.o Postprocessor.o \
	  Preprocessor.o SolveEqu.o RNetwork.o MainPro.o \
	                   
neca : $(objects)        
	${C++} ${PARFLAGS} -o neca $(objects)
	              
MainPro.o : MainPro.cpp Neca.h
	${C++} ${CXXFLAGS} MainPro.cpp -c

Global_Stiff_Matrix.o : Global_Stiff_Matrix.cpp
	${C++} ${PARFLAGS} ${CXXFLAGS} $< -c 

Algorithm_FEM.o : Algorithm_FEM.cpp
	${C++} ${PARFLAGS} ${CXXFLAGS} $< -c 
		
# 伪目标文件
.PHONY : clean
clean :
	-rm neca $(objects)
