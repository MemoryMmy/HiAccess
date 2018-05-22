# -*- coding: UTF-8
import time
import os

cmd ="./facilityAccess  --road ./data/龙山路网/龙山路网.shp --facility ./data/龙山县卫生局二级单位/龙山卫生局二级单位.shp --resolution 1.3234 "
N=10

time1 = time.time()
for i in range(N):
	#print i
	os.system(cmd)

time2 = time.time()
print str((time2-time1)/N) +'s('+ str(N)+' times average)'
