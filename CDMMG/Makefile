HEADERS =-I./ -I/usr/local/include/ -I/home/dbrg/APPS/boost_1_64_0/include
LIBS =-L./ -L/usr/local/lib -L/home/dbrg/APPS/boost_1_64_0/stage/lib

facilityAccess:  facilityaccess.cpp
	mpicxx $(HEADERS) $(LIBS) -std=c++11 -fopenmp -Wall -O3 -o facilityAccess facilityaccess.cpp -lgdal -lmpi
clean:
	rm  facilityAccess
