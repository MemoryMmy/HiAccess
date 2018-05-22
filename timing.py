# -*- coding: UTF-8
# timing program
import time
import os

cmd ="mpirun -np 4 ./CDM/facilityAccess  --road ./data/dataset1/roads.shp --facility ./data/dataset1/education.shp --resolution 100 --tolarence 10 --maxBound 3000 --rate 2 --output ./results/dataset1_100m_CDM.tif"
N=10

time1 = time.time()
for i in range(N):
	os.system(cmd)

time2 = time.time()
print str((time2-time1)/N) +'s('+ str(N)+' times average)'
