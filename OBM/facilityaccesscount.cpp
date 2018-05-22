/*
* date: 2017-08-10
* author: Mengyu Ma@National University of Defense Technology
* e-mail: mamengyu10@nudt.edu.cn
* description: Apply HiAccess to the opportunity-based measurement.
* run: 
* 	 $ mpirun -np 4 ./facilityAccesscount  --road ../data/Huayuan/Roads/CountyRoad.shp,../data/Huayuan/Roads/RuralRoad.shp,../data/Huayuan/Roads/NationalHighway.shp --facility ../data/Huayuan/Jobs/FactoryEnterprise.shp --resolution 0.001 --tolarence 0.0001 --maxBound 0.03 --rate 2 --output ./job.tif
*/



#include "ogrsf_frmts.h"
#include "ogr_p.h"
#include "cpl_conv.h"
#include "cpl_string.h"
#include <stdio.h>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <dirent.h>
#include <omp.h>
#include <sys/time.h>
#include <regex.h>
#include <math.h>
#include "mpi.h"

#include <queue>
#include <utility>

#include <boost/unordered_map.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#define MAX_NODE_SIZE 8
#define MAX_SPLIT_NUM 4 
#define MAX_NEIGHBOR_NUM 8
#define MAX_DOUBLE 1000000000
#define MAX_ROAD_FILE 256
#define MAX_POINT_NUM 100 

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
namespace bgm = boost::geometry::model;


void GetFilelist(char* argv, char* result[], char* flag, int& count)
{
	char* string = strdup(argv); 
	char* p;
	int i = 0;
	while((p = strsep(&string, flag)) != NULL)
	{
		result[i] = p;
		i++;
	}
	result[i] = string;
	count = i;
}

inline double distanceP2P(double x1,double y1,double x2,double y2)
{
	double dx = x2 - x1;
	double dy = y2 - y1;
	return sqrt(dx*dx + dy*dy);
}

inline double distanceP2PSquare(double x1,double y1,double x2,double y2)
{
	double dx = x2 - x1;
	double dy = y2 - y1;
	return dx*dx + dy*dy;
}

bool iscommonLinenear(const int s0,const double d0,const int s1,const double d1,double maxBound,
const double* seglengthList,double &distance)
{
	distance=d0+d1;
	if(!(distance<maxBound))
		return false;
	if(s0<s1)
	{
		for(int i=s0;i<s1;i++)
		{
			
			distance+=seglengthList[i];
			if(!(distance<maxBound))
				return false;
		}	
	}else
	{
		for(int i=s1;i<s0;i++)
		{
			
			distance+=seglengthList[i];
			if(!(distance<maxBound))
				return false;
		}
	}
	return true;
}


bool getDistance(const double x0,const double y0,const double x1,const double y1,const double px,const double py, 
const double* seglengthList,const int segmentIndex,const int seglen,double rate, double maxBound,
double &segFDistance,int &segFIndex,double &headDistance,double &endDistance)
{
	
	double x,y;
	if(x0==x1)
	{
		x=x0;
		y=py;
	}
	else if(y0==y1)
	{
		y=y0;
		x=px;
	}
	else
	{
		double k = (y1-y0)/(x1-x0);
		x=(k*k*x0 + k*(py-y0) + px)/(k*k+1);
		y=k*(x-x0)+y0;
	}
	double dSquare0=distanceP2PSquare(px,py,x0,y0);
	double dSquare1=distanceP2PSquare(px,py,x1,y1);
	double dSquare=distanceP2PSquare(px,py,x,y);
	if(dSquare<dSquare0&&dSquare<dSquare1)
	{
		double d=sqrt(dSquare)*rate;
		headDistance=d+distanceP2P(x,y,x0,y0);
		endDistance=d+distanceP2P(x,y,x1,y1);
		if(headDistance<endDistance)
		{
			segFDistance=headDistance;
			segFIndex=segmentIndex;
		}else
		{
			segFDistance=endDistance;
			segFIndex=segmentIndex+1;
		}
		if(!(segFDistance<maxBound))
			return false;
		for(int i=0;i<segmentIndex;i++)
			headDistance+=seglengthList[i];
		for(int i=segmentIndex+1;i<seglen;i++)
			endDistance+=seglengthList[i];	
	}
	else if(dSquare0>dSquare1)
	{
		double d=sqrt(dSquare1)*rate;
		headDistance=d;
		endDistance=d;
		segFDistance=d;
		segFIndex=segmentIndex+1;
		if(!(segFDistance<maxBound))
			return false;
		for(int i=0;i<segmentIndex+1;i++)
			headDistance+=seglengthList[i];
		for(int i=segmentIndex+1;i<seglen;i++)
			endDistance+=seglengthList[i];	
	}
	else
	{
		double d=sqrt(dSquare0)*rate;
		headDistance=d;
		endDistance=d;
		segFDistance=d;
		segFIndex=segmentIndex;
		if(!(segFDistance<maxBound))
			return false;
		for(int i=0;i<segmentIndex;i++)
			headDistance+=seglengthList[i];
		for(int i=segmentIndex;i<seglen;i++)
			endDistance+=seglengthList[i];							
	}
	return true;
}

//Dijkstra
class cmp
{
	public:
		bool operator() (std::pair<int,double> a,std::pair<int,double> b) 
		{
			return a.second>b.second;
		}
	
};
			
void shortestDistance(const int n,const int v0,double** adjacencydistList1,int **adjacencyList1,
int * neighborcList,double *dist) 
{ 
	std::priority_queue< std::pair< int,double >, std::vector<std::pair<int,double> >, cmp > priorQ;
	double (*adjacencydistList)[MAX_NEIGHBOR_NUM]  =(double(*)[MAX_NEIGHBOR_NUM])adjacencydistList1;
	int (*adjacencyList)[MAX_NEIGHBOR_NUM]  =(int(*)[MAX_NEIGHBOR_NUM])adjacencyList1;
    bool *visited=(bool *)malloc(sizeof(bool)*n);
	memset(dist,0x43,n*sizeof(double));
	memset(visited,false,n);
	int u;
	double disttmp;
    dist[v0]=0;
	    
	for(int i=0;i<n;i++) 
		priorQ.push(std::make_pair(i,dist[i]));
    
    while(!priorQ.empty())
    {
		u = priorQ.top().first;
		priorQ.pop();
		if (dist[u]>MAX_DOUBLE){
		    break;
		}
		if(visited[u])
		{
			continue;
		}
		visited[u]=true;
        for(int i=0;i<neighborcList[u];i++) 
        {
			disttmp=dist[u]+adjacencydistList[u][i];
            if(!visited[adjacencyList[u][i]])
            {
				disttmp=dist[u]+adjacencydistList[u][i];
				if(disttmp<dist[adjacencyList[u][i]])
				{
					dist[adjacencyList[u][i]]=disttmp;
					priorQ.push(std::make_pair(adjacencyList[u][i],disttmp));
			    }
			}    
        } 
  
    }
}


void dfsGraph(unsigned short* subgraphFlag,const int currentI,const int graphCount,
const int* neighborcList,int ** adjacencyList1, int &graphvertexCount)
{
	int (*adjacencyList)[MAX_NEIGHBOR_NUM]  =(int(*)[MAX_NEIGHBOR_NUM])adjacencyList1;
	for(int i=0;i<neighborcList[currentI];i++)
		if(subgraphFlag[adjacencyList[currentI][i]]<graphCount)
		{
			//~ printf("----i:%d \n",adjacencyList[currentI][i]);
			subgraphFlag[adjacencyList[currentI][i]]=graphCount;
			graphvertexCount++;
			dfsGraph(subgraphFlag,adjacencyList[currentI][i],graphCount,neighborcList,(int **)adjacencyList,graphvertexCount);
		}
}

void Usage()
{
    printf( "Usage:           [--road:       道路文件(shp)       ]\n"
	    "                     [--facility:   设施文件(shp)       ]\n"
	    "                     [--resolution: 分辨率              ]\n"
	    "                     [--tolarence:  生成路网容差         ]\n"
	    "                     [--maxBound:   最大搜索范围         ]\n"
	    "                     [--rate:       无道路代价比         ]\n"
	    "                     [--output:     输出文件(tif)       ]\n"	    	    
	    );
}



int main( int nArgc, char ** papszArgv )
{
	typedef bgm::d2::point_xy<double> point;
	typedef bgm::box<point> box;
	typedef bgm::segment<point> segment; 
	//~ typedef bgm::linestring<point> linestring;  
    typedef std::pair<point, unsigned> valueP;
    //~ typedef std::pair<box, unsigned> valueb;  
    typedef boost::tuple<segment,unsigned,unsigned> valueS;
    typedef boost::unordered_multimap<unsigned,int*> umultimap;

	int myId, numProcs;
	MPI_Init(&nArgc,&papszArgv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myId);
	MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
	double t1=MPI_Wtime();

	const char* roadShpList = NULL;     
	const char* facilityShp = NULL; 
	const char* rasterFile=NULL; 
	double resolution = -1;	  
	double tolarence = 0.00001;
	double maxBound = 0.1;
	double rate =2;
	const char* roadAppidListS = NULL; 
	const char* facilityAppid =NULL;
	
	char shpPath[256] = "/cluster/higis/data/shape/";
	nArgc = OGRGeneralCmdLineProcessor( nArgc, &papszArgv, 0 );

	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	CPLSetConfigOption("SHAPE_ENCODING", "UTF-8");
   
	for( int iArg = 1; iArg < nArgc; iArg++ )
	{
		if( EQUAL(papszArgv[iArg], "--road") )
	    {
			roadShpList = papszArgv[iArg+1];
			
		}
		else if( EQUAL(papszArgv[iArg], "--facility") )
	    {
			facilityShp = papszArgv[iArg+1];
			
		}
		else if( EQUAL(papszArgv[iArg], "--resolution") )
	    {
			resolution = atof(papszArgv[iArg+1]);
			
		}
		else if( EQUAL(papszArgv[iArg], "--tolarence") )
	    {
			tolarence = atof(papszArgv[iArg+1]);
			
		}
		else if( EQUAL(papszArgv[iArg], "--maxBound") )
	    {
			maxBound = atof(papszArgv[iArg+1]);
			
		}
		else if( EQUAL(papszArgv[iArg], "--rate") )
	    {
			rate = atof(papszArgv[iArg+1]);
			
		}
		else if( EQUAL(papszArgv[iArg], "--output") )
	    {
			rasterFile = papszArgv[iArg+1];
		}
		else if( EQUAL(papszArgv[iArg], "--roadAppid") )
	    {
			roadAppidListS = papszArgv[iArg+1];
		}
		else if( EQUAL(papszArgv[iArg], "--facilityAppid") )
	    {
			facilityAppid = papszArgv[iArg+1];
		}	
	}
	
	if (roadShpList == NULL || facilityShp == NULL || resolution == -1 ||rasterFile==NULL)
	{
		
		if(myId==0)
			Usage();
		return -3;
	}
	
	if(myId==0)
	{
		printf("--road: %s\n", roadShpList);
		printf("--facility: %s\n", facilityShp);
		printf("--resolution: %lf\n", resolution);
		printf("--tolarence: %lf\n", tolarence);
		printf("--maxBound: %lf\n", maxBound);
		printf("--rate: %lf\n", rate);
		printf("--output: %s\n", rasterFile);
		
	}
	
	
	bgi::rtree< valueP, bgi::quadratic<MAX_NODE_SIZE> > rtreeVertex;
	std::vector<int> edgeheadListpre;
	std::vector<int> edgeendListpre;
	int vertexIndex=0;
	double minXOut=MAX_DOUBLE,minYOut=MAX_DOUBLE,maxXout=-1*MAX_DOUBLE,maxYOut=-1*MAX_DOUBLE;
	
	int roadCount = 0;
	char* roadList[MAX_ROAD_FILE];
	char* roadAppidList[MAX_ROAD_FILE];
	bool isroadinPG = false;
	GetFilelist( (char*)roadShpList, roadList, (char*)",", roadCount);
	if(roadAppidListS!=NULL)
	{
		isroadinPG=true;
		GetFilelist( (char*)roadAppidListS, roadAppidList, (char*)",", roadCount);
		
	}
	if(facilityAppid!=NULL)
	{
		char *tmp=new char[256];
		char s1[2]="/";
		char s2[5]=".shp";
		sprintf(tmp,"%s%s%s%s%s",shpPath,facilityAppid,s1,facilityShp,s2);
		facilityShp=tmp;	
	}
	OGRRegisterAll();
	GDALDataset* facilityDS =(GDALDataset*)GDALOpenEx(facilityShp, GDAL_OF_VECTOR,NULL, NULL, NULL );
	if(facilityDS==NULL)
	{
		printf("[ERROR] Open facilityShp failed.\n");
		exit(1);
	}
	OGRLayer* facilityLayer=facilityDS->GetLayer(0);
	OGRFeature *facilityFeature;
	OGRSpatialReference *pOgrSRS = NULL;
	pOgrSRS = facilityLayer->GetSpatialRef();
	
	OGRLayer** roadLayer= new OGRLayer*[roadCount];
	OGRFeature *roadFeature;
	for(int i=0;i<roadCount;i++)
	{
		char *roadShp;
		if(isroadinPG)
		{
			char *tmp=new char[256];
			char s1[2]="/";
			char s2[5]=".shp";
			sprintf(tmp,"%s%s%s%s%s",shpPath,roadAppidList[i],s1,roadList[i],s2);
			roadShp=tmp;
		}else
			roadShp=roadList[i];
		GDALDataset* roadDS =(GDALDataset*)GDALOpenEx(roadShp, GDAL_OF_VECTOR,NULL, NULL, NULL );
		if(roadDS==NULL)
		{
			printf("[ERROR] Open roadShp failed.\n");
			exit(1);
		}
		roadLayer[i]=roadDS->GetLayer(0);
		OGREnvelope env;
		roadLayer[i]->GetExtent(&env);
		minXOut=minXOut<env.MinX? minXOut:env.MinX;
		minYOut=minYOut<env.MinY? minYOut:env.MinY;
		maxXout=maxXout>env.MaxX? maxXout:env.MaxX;
		maxYOut=maxYOut>env.MaxY? maxYOut:env.MaxY;

		roadLayer[i]->ResetReading(); 
		while((roadFeature=roadLayer[i]-> GetNextFeature())!= NULL)
		{
			OGRGeometry *poGeometry=roadFeature->GetGeometryRef();
			int eType = wkbFlatten(poGeometry->getGeometryType());
			
			if(eType == wkbLineString)
			{
				
				OGRLineString* pOGRLineString=(OGRLineString*) poGeometry;
				int pointCount = pOGRLineString->getNumPoints();
				
				double hx=pOGRLineString->getX(0);
				double hy=pOGRLineString->getY(0);
				double ex=pOGRLineString->getX(pointCount-1);
				double ey=pOGRLineString->getY(pointCount-1);
				
				box headBuffer;
				box endBuffer;
				bg::buffer(bg::return_envelope<box>(point(hx,hy)),headBuffer,tolarence);
				bg::buffer(bg::return_envelope<box>(point(ex,ey)),endBuffer,tolarence);
	
				std::vector<valueP> headResult;
				std::vector<valueP> endResult;
				
				rtreeVertex.query(bgi::intersects(headBuffer), std::back_inserter(headResult));
				if(headResult.size()<1)
				{
					rtreeVertex.insert(std::make_pair(point(hx,hy), vertexIndex));
					edgeheadListpre.push_back(vertexIndex);
					vertexIndex++;
				}else{
					edgeheadListpre.push_back(headResult.front().second);
				}
				
				
				rtreeVertex.query(bgi::intersects(endBuffer), std::back_inserter(endResult));
				if(endResult.size()<1)
				{
					rtreeVertex.insert(std::make_pair(point(ex,ey), vertexIndex));
					edgeendListpre.push_back(vertexIndex);
					vertexIndex++;
				}else{
					edgeendListpre.push_back(endResult.front().second);
				}
			}
			else if( eType == wkbMultiLineString)
			{
				//~ printf("wkbMultiLineString");
				OGRMultiLineString* pOGRMultiLineString=(OGRMultiLineString*) poGeometry;
				int lineCount =pOGRMultiLineString->getNumGeometries();
				for(int j=0;j<lineCount;j++)
				{
					OGRLineString* pOGRLineString=(OGRLineString*)pOGRMultiLineString->getGeometryRef(j);
					int pointCount = pOGRLineString->getNumPoints();
				
					double hx=pOGRLineString->getX(0);
					double hy=pOGRLineString->getY(0);
					double ex=pOGRLineString->getX(pointCount-1);
					double ey=pOGRLineString->getY(pointCount-1);
					
					box headBuffer;
					box endBuffer;
					bg::buffer(bg::return_envelope<box>(point(hx,hy)),headBuffer,tolarence);
					bg::buffer(bg::return_envelope<box>(point(ex,ey)),endBuffer,tolarence);
		
					std::vector<valueP> headResult;
					std::vector<valueP> endResult;
					
					rtreeVertex.query(bgi::intersects(headBuffer), std::back_inserter(headResult));
					if(headResult.size()<1)
					{
						rtreeVertex.insert(std::make_pair(point(hx,hy), vertexIndex));
						edgeheadListpre.push_back(vertexIndex);
						vertexIndex++;
					}else{
						edgeheadListpre.push_back(headResult.front().second);
					}
					
					
					rtreeVertex.query(bgi::intersects(endBuffer), std::back_inserter(endResult));
					if(endResult.size()<1)
					{
						rtreeVertex.insert(std::make_pair(point(ex,ey), vertexIndex));
						edgeendListpre.push_back(vertexIndex);
						vertexIndex++;
					}else{
						edgeendListpre.push_back(endResult.front().second);
					}
							
				}
			}
		}	
		
	}
		


	double (*adjacencydistList)[MAX_NEIGHBOR_NUM] = (double(*)[MAX_NEIGHBOR_NUM])malloc(MAX_NEIGHBOR_NUM*vertexIndex*sizeof(double));
	int (*adjacencyList)[MAX_NEIGHBOR_NUM] = (int(*)[MAX_NEIGHBOR_NUM])malloc(MAX_NEIGHBOR_NUM*vertexIndex*sizeof(int));
	int* neighborcList = (int*)malloc(sizeof(int)*vertexIndex);
	memset(neighborcList,0,sizeof(int)*vertexIndex);

	int* edgeheadList=(int*)malloc(sizeof(int)*(vertexIndex*MAX_SPLIT_NUM));
	int* edgeendList=(int*)malloc(sizeof(int)*(vertexIndex*MAX_SPLIT_NUM));
	double** edgexList=new double*[vertexIndex*MAX_SPLIT_NUM];
	double** edgeyList=new double*[vertexIndex*MAX_SPLIT_NUM];
	
	double** edgelengthList=new double*[vertexIndex*MAX_SPLIT_NUM];
	int* edgelenList=(int*)malloc(sizeof(int)*vertexIndex*MAX_SPLIT_NUM);
	
	int edgeLen=0; 
	int lineIndex=0;
	for(int i=0;i<roadCount;i++)
	{

	
	
		roadLayer[i]->ResetReading(); 
		
		while((roadFeature=roadLayer[i]-> GetNextFeature())!= NULL)
		{
			
			OGRGeometry *poGeometry=roadFeature->GetGeometryRef();
			int eType = wkbFlatten(poGeometry->getGeometryType());
			if(eType == wkbLineString)
			{	
	
				OGRLineString* pOGRLineString=(OGRLineString*) poGeometry;
				int pointCount = pOGRLineString->getNumPoints();
				edgeheadList[edgeLen]=edgeheadListpre[lineIndex];
				
				edgexList[edgeLen]=(double*)malloc(sizeof(double)*(pointCount));
				edgeyList[edgeLen]=(double*)malloc(sizeof(double)*(pointCount));
				
				edgelengthList[edgeLen]=(double*)malloc(sizeof(double)*(pointCount-1));
				
				int len=1;
				double edgeLength=0;
				edgexList[edgeLen][0]=pOGRLineString->getX(0);
				edgeyList[edgeLen][0]=pOGRLineString->getY(0);
				for(int i=1;i<pointCount-1;i++)
				{
	
					double px=pOGRLineString->getX(i);
					double py=pOGRLineString->getY(i);
					box pointBuffer;
					bg::buffer(bg::return_envelope<box>(point(px,py)),pointBuffer,tolarence);
					std::vector<valueP> pointResult;
					rtreeVertex.query(bgi::intersects(pointBuffer), std::back_inserter(pointResult));
	
					edgexList[edgeLen][len]=px;
					edgeyList[edgeLen][len]=py;
					edgelengthList[edgeLen][len-1]=distanceP2P(edgexList[edgeLen][len-1],edgeyList[edgeLen][len-1],px,py);
					edgeLength+=edgelengthList[edgeLen][len-1];		
				
					if(pointResult.size()<1)
					{
						len++;			
					}else
					{
						edgeendList[edgeLen]=pointResult.front().second;
						edgelenList[edgeLen]=len;
	
						if(edgeheadList[edgeLen]!=edgeendList[edgeLen])
						{
							adjacencydistList[edgeheadList[edgeLen]][neighborcList[edgeheadList[edgeLen]]]=edgeLength;
							adjacencyList[edgeheadList[edgeLen]][neighborcList[edgeheadList[edgeLen]]]= edgeendList[edgeLen];
							neighborcList[edgeheadList[edgeLen]]++;
							adjacencydistList[edgeendList[edgeLen]][neighborcList[edgeendList[edgeLen]]]=edgeLength;
							adjacencyList[edgeendList[edgeLen]][neighborcList[edgeendList[edgeLen]]]=edgeheadList[edgeLen];
							neighborcList[edgeendList[edgeLen]]++;
						}				
	
						edgeLen++;
						edgeheadList[edgeLen]=pointResult.front().second;
						
						edgexList[edgeLen]=(double*)malloc(sizeof(double)*(pointCount-len));
						edgeyList[edgeLen]=(double*)malloc(sizeof(double)*(pointCount-len));
						edgelengthList[edgeLen]=(double*)malloc(sizeof(double)*(pointCount-len-1));
						len=1;
						edgeLength=0;
						edgexList[edgeLen][0]=px;
						edgeyList[edgeLen][0]=py;
					}
	
				}
				edgeendList[edgeLen]=edgeendListpre[lineIndex];
				edgelenList[edgeLen]=len;
				edgexList[edgeLen][len]=pOGRLineString->getX(pointCount-1);
				edgeyList[edgeLen][len]=pOGRLineString->getY(pointCount-1);
				edgelengthList[edgeLen][len-1]=distanceP2P(edgexList[edgeLen][len-1],edgeyList[edgeLen][len-1],edgexList[edgeLen][len],edgeyList[edgeLen][len]);
				edgeLength+=edgelengthList[edgeLen][len-1];	
	
				if(edgeheadList[edgeLen]!=edgeendList[edgeLen])
				{
					adjacencydistList[edgeheadList[edgeLen]][neighborcList[edgeheadList[edgeLen]]]=edgeLength;
					adjacencyList[edgeheadList[edgeLen]][neighborcList[edgeheadList[edgeLen]]]= edgeendList[edgeLen];
					neighborcList[edgeheadList[edgeLen]]++;
					adjacencydistList[edgeendList[edgeLen]][neighborcList[edgeendList[edgeLen]]]=edgeLength;
					adjacencyList[edgeendList[edgeLen]][neighborcList[edgeendList[edgeLen]]]=edgeheadList[edgeLen];
					neighborcList[edgeendList[edgeLen]]++;
				}			
	
				edgeLen++;
				lineIndex++;
			}
			else if( eType == wkbMultiLineString)
			{
				//~ printf("wkbMultiLineString");
				OGRMultiLineString* pOGRMultiLineString=(OGRMultiLineString*) poGeometry;
				int lineCount =pOGRMultiLineString->getNumGeometries();
				for(int j=0;j<lineCount;j++)
				{
					OGRLineString* pOGRLineString=(OGRLineString*)pOGRMultiLineString->getGeometryRef(j);
					int pointCount = pOGRLineString->getNumPoints();
					edgeheadList[edgeLen]=edgeheadListpre[lineIndex];
					
					edgexList[edgeLen]=(double*)malloc(sizeof(double)*(pointCount));
					edgeyList[edgeLen]=(double*)malloc(sizeof(double)*(pointCount));
					
					edgelengthList[edgeLen]=(double*)malloc(sizeof(double)*(pointCount-1));
					
					int len=1;
					double edgeLength=0;
					edgexList[edgeLen][0]=pOGRLineString->getX(0);
					edgeyList[edgeLen][0]=pOGRLineString->getY(0);
					for(int i=1;i<pointCount-1;i++)
					{
		
						double px=pOGRLineString->getX(i);
						double py=pOGRLineString->getY(i);
						box pointBuffer;
						bg::buffer(bg::return_envelope<box>(point(px,py)),pointBuffer,tolarence);
						std::vector<valueP> pointResult;
						rtreeVertex.query(bgi::intersects(pointBuffer), std::back_inserter(pointResult));
		
						edgexList[edgeLen][len]=px;
						edgeyList[edgeLen][len]=py;
						edgelengthList[edgeLen][len-1]=distanceP2P(edgexList[edgeLen][len-1],edgeyList[edgeLen][len-1],px,py);
						edgeLength+=edgelengthList[edgeLen][len-1];		
					
						if(pointResult.size()<1)
						{
							len++;			
						}else
						{
							edgeendList[edgeLen]=pointResult.front().second;
							edgelenList[edgeLen]=len;
		
							if(edgeheadList[edgeLen]!=edgeendList[edgeLen])
							{
								adjacencydistList[edgeheadList[edgeLen]][neighborcList[edgeheadList[edgeLen]]]=edgeLength;
								adjacencyList[edgeheadList[edgeLen]][neighborcList[edgeheadList[edgeLen]]]= edgeendList[edgeLen];
								neighborcList[edgeheadList[edgeLen]]++;
								adjacencydistList[edgeendList[edgeLen]][neighborcList[edgeendList[edgeLen]]]=edgeLength;
								adjacencyList[edgeendList[edgeLen]][neighborcList[edgeendList[edgeLen]]]=edgeheadList[edgeLen];
								neighborcList[edgeendList[edgeLen]]++;
							}				
		
							edgeLen++;
							edgeheadList[edgeLen]=pointResult.front().second;
							
							edgexList[edgeLen]=(double*)malloc(sizeof(double)*(pointCount-len));
							edgeyList[edgeLen]=(double*)malloc(sizeof(double)*(pointCount-len));
							edgelengthList[edgeLen]=(double*)malloc(sizeof(double)*(pointCount-len-1));
							len=1;
							edgeLength=0;
							edgexList[edgeLen][0]=px;
							edgeyList[edgeLen][0]=py;
						}
		
					}
					edgeendList[edgeLen]=edgeendListpre[lineIndex];
					edgelenList[edgeLen]=len;
					edgexList[edgeLen][len]=pOGRLineString->getX(pointCount-1);
					edgeyList[edgeLen][len]=pOGRLineString->getY(pointCount-1);
					edgelengthList[edgeLen][len-1]=distanceP2P(edgexList[edgeLen][len-1],edgeyList[edgeLen][len-1],edgexList[edgeLen][len],edgeyList[edgeLen][len]);
					edgeLength+=edgelengthList[edgeLen][len-1];	
		
					if(edgeheadList[edgeLen]!=edgeendList[edgeLen])
					{
						adjacencydistList[edgeheadList[edgeLen]][neighborcList[edgeheadList[edgeLen]]]=edgeLength;
						adjacencyList[edgeheadList[edgeLen]][neighborcList[edgeheadList[edgeLen]]]= edgeendList[edgeLen];
						neighborcList[edgeheadList[edgeLen]]++;
						adjacencydistList[edgeendList[edgeLen]][neighborcList[edgeendList[edgeLen]]]=edgeLength;
						adjacencyList[edgeendList[edgeLen]][neighborcList[edgeendList[edgeLen]]]=edgeheadList[edgeLen];
						neighborcList[edgeendList[edgeLen]]++;
					}			
		
					edgeLen++;
					lineIndex++;
				
				}
			
			}
		}
	
	
	}
	
	if(myId==0)
	printf("vertex count:%d edge count:%d\n",vertexIndex,edgeLen);
	
	//~ double t2= MPI_Wtime();
	//~ printf("[network build] Successfully. cores:%d coreid:%d cost:<total>%f</total> s\n",numProcs,myId, t2-t1);
	//create rtreeEdge
	//max connected subgraph NUM 255
	unsigned short* subgraphFlag= (unsigned short*)malloc(vertexIndex*2);
	memset(subgraphFlag,0,vertexIndex*2);
	int graphvertexCount=0;
	int graphCount=0;
	while(graphvertexCount<vertexIndex)
	{
		int i;
		for(i=0;i<vertexIndex;i++)
		{
			
			if(subgraphFlag[i]<1)
				break;
		}
		graphCount++;
		graphvertexCount++;
		subgraphFlag[i]=graphCount;
		dfsGraph(subgraphFlag,i,graphCount,neighborcList,(int **)adjacencyList,graphvertexCount);
	}
	if(myId==0)
	printf("graphCount:%d\n",graphCount);

	bgi::rtree< valueS, bgi::quadratic<MAX_NODE_SIZE> > rtreeList[graphCount];

	for(int i=0;i<edgeLen;i++)
	{
		int graphIndex =subgraphFlag[edgeheadList[i]]-1;
		for(int j=0;j<edgelenList[i];j++)
		{	
			rtreeList[graphIndex].insert(boost::make_tuple(segment(point(edgexList[i][j],edgeyList[i][j]),
			point(edgexList[i][j+1],edgeyList[i][j+1])), i,j));
		}
	}
	

	int facilityCount=facilityLayer->GetFeatureCount();
	if(myId==0)
	printf("facilityCount :%d \n",facilityCount);
	
	
	double* facilitypxPolygon = (double*)malloc(facilityCount*MAX_POINT_NUM*sizeof(double));
	double* facilitypyPolygon = (double*)malloc(facilityCount*MAX_POINT_NUM*sizeof(double));
	int facilityPolygonCount=0;
	for(int i=0;i<facilityCount;i++)
	{
		facilityFeature=facilityLayer->GetFeature(i);
		OGRGeometry *poGeometry=facilityFeature->GetGeometryRef();
		int eType = wkbFlatten(poGeometry->getGeometryType());
	
		if( eType == wkbMultiPoint)
		{
			//~ printf("wkbMultiPoint\n");
			OGRMultiPoint* pOGRMultiPoint=(OGRMultiPoint*) poGeometry;
			int pointCount=pOGRMultiPoint->getNumGeometries();
			for(int k=0;k<pointCount;k++)
			{
				OGRPoint* pOGRPoint=(OGRPoint*) pOGRMultiPoint->getGeometryRef(k);
				double px = pOGRPoint->getX();
				double py = pOGRPoint->getY();
				facilitypxPolygon[facilityPolygonCount]=px;
				facilitypyPolygon[facilityPolygonCount]=py;
				facilityPolygonCount++;
			
			}		
		}
		else if( eType ==wkbPolygon)
		{
			//~ printf("wkbPolygon\n");
			OGRPolygon* pOGRPolygon=(OGRPolygon*) poGeometry;
			OGRLinearRing *pLinearRing = pOGRPolygon->getExteriorRing();
			int pointCount=pLinearRing->getNumPoints();
			double px0,py0,px1,py1;
			px0 =pLinearRing->getX(0);
			py0 =pLinearRing->getY(0);
			facilitypxPolygon[facilityPolygonCount]=px0;
			facilitypyPolygon[facilityPolygonCount]=py0;
			facilityPolygonCount++;
			for (int j=1;j<pointCount-1;j++)
			{
				px1 =pLinearRing->getX(j);
				py1 =pLinearRing->getY(j);
				if(fabs(px1-px0)>resolution||fabs(py1-py0)>resolution)
				{
					px0=px1;
					py0=py1;
					facilitypxPolygon[facilityPolygonCount]=px0;
					facilitypyPolygon[facilityPolygonCount]=py0;
					facilityPolygonCount++;
				}
			}
		}
		else if( eType ==wkbMultiPolygon)
		{
			//~ printf("wkbMultiPolygon\n");
			OGRMultiPolygon* pOGRMultiPolygon=(OGRMultiPolygon*) poGeometry;
			int polygonCount=pOGRMultiPolygon->getNumGeometries();
			for(int k=0;polygonCount;k++)
			{
				OGRPolygon* pOGRPolygon=(OGRPolygon*)pOGRMultiPolygon->getGeometryRef(k);
				OGRLinearRing *pLinearRing = pOGRPolygon->getExteriorRing();
				int pointCount=pLinearRing->getNumPoints();
				double px0,py0,px1,py1;
				px0 =pLinearRing->getX(0);
				py0 =pLinearRing->getY(0);
				facilitypxPolygon[facilityPolygonCount]=px0;
				facilitypyPolygon[facilityPolygonCount]=py0;
				facilityPolygonCount++;
				for (int j=1;j<pointCount-1;j++)
				{
					px1 =pLinearRing->getX(j);
					py1 =pLinearRing->getY(j);
					if(fabs(px1-px0)>resolution||fabs(py1-py0)>resolution)
					{
						px0=px1;
						py0=py1;
						facilitypxPolygon[facilityPolygonCount]=px0;
						facilitypyPolygon[facilityPolygonCount]=py0;
						facilityPolygonCount++;
					}
				}
			}
			
		}
	}
	//~ if(myId==0)
	//~ printf("facilityPolygonCount:%d\n",facilityPolygonCount);
	
	
	
	int totalCount=facilityPolygonCount/numProcs+facilityCount/numProcs+2;
	int distenceListsize=sizeof(int)+vertexIndex*totalCount*sizeof(double)+
	sizeof(int)+(sizeof(int)*3+sizeof(double))*totalCount*graphCount;
	int* distenceList=(int *)malloc(distenceListsize);
	distenceList[0]=0;
	double* distenceListBegin=(double *)(distenceList+ 1);
	memset(distenceListBegin,0x43,vertexIndex*totalCount*sizeof(double));
	
	
	
	int* fdList=(int*)(distenceList+1+vertexIndex*totalCount*2);
	fdList[0]=0;
	int* fdListBegin=fdList+1;
	
	void* distenceListFinal=(void *)malloc(distenceListsize*numProcs);
	int totalCountExact=0;
	for(int graphIndex=0;graphIndex<graphCount;graphIndex++)
	{	
		bgi::rtree< valueS, bgi::quadratic<MAX_NODE_SIZE> >* rtreesubgraphEdge;
		rtreesubgraphEdge=&rtreeList[graphIndex];
		for(int i=myId;i<facilityPolygonCount;i+=numProcs)
		{
			if(graphIndex<1)
			totalCountExact++;
			
			if(distenceList[0]==totalCountExact)
			distenceList[0]=0;
			
			double px = facilitypxPolygon[i];
			double py = facilitypyPolygon[i];
			
			box pointBuffer;
			bg::buffer(bg::return_envelope<box>(point(px,py)),pointBuffer,maxBound);
			std::vector<valueS> segmentResult;
			(*rtreesubgraphEdge).query(bgi::intersects(pointBuffer)&&bgi::nearest(point(px,py),1), std::back_inserter(segmentResult));
			if(segmentResult.size()>0)
			{
				
				//~ printf("-----%ld\n",segmentResult.size());
				segment pSegment= boost::get<0>(segmentResult.front());
				double x0= bg::get<0,0>(pSegment);
				double y0= bg::get<0,1>(pSegment);
				double x1= bg::get<1,0>(pSegment);
				double y1= bg::get<1,1>(pSegment);
				int edgeIndex= boost::get<1>(segmentResult.front());
				int segmentIndex=boost::get<2>(segmentResult.front());
				
				int headId = edgeheadList[edgeIndex];
				int endId = edgeendList[edgeIndex];
				
				//caculate headDistance, endDistance;
				double headDistance, endDistance;
				double segFDistance;
				int segFIndex;
				bool s=getDistance(x0,y0,x1,y1,px,py,edgelengthList[edgeIndex],segmentIndex,edgelenList[edgeIndex],rate,maxBound,
				segFDistance,segFIndex,headDistance,endDistance);
				if(s)
				{	
					fdListBegin[5*fdList[0]]=edgeIndex;
					fdListBegin[5*fdList[0]+1]=segFIndex;
					fdListBegin[5*fdList[0]+2]=distenceList[0]+myId*totalCount;
					*((double *)(&fdListBegin[5*fdList[0]+3]))=segFDistance;
					fdList[0]++;
					
					double* distenceheadList=(double*)malloc(vertexIndex*sizeof(double));
					double* distenceendList=(double*)malloc(vertexIndex*sizeof(double));

					shortestDistance(vertexIndex,headId, (double **)adjacencydistList,(int **)adjacencyList,neighborcList,distenceheadList);
					shortestDistance(vertexIndex,endId, (double **)adjacencydistList,(int **)adjacencyList,neighborcList,distenceendList);
					double* currentBegin=(double*)(distenceListBegin+distenceList[0]);
					for(int j=0;j<vertexIndex;j++)
					{
						double hd=distenceheadList[j]+headDistance;
						double ed=distenceendList[j]+endDistance;
						double mind= hd<ed? hd:ed;

						if(maxBound>mind)
							currentBegin[j*totalCount]=mind;	
					}
				}

			}
			distenceList[0]++;
		}
			
		for(int i=myId;i<facilityCount;i+=numProcs)
		{
			facilityFeature=facilityLayer->GetFeature(i);
			OGRGeometry *poGeometry=facilityFeature->GetGeometryRef();
			int eType = wkbFlatten(poGeometry->getGeometryType());
			//~ printf("!distenceList[0]:%d facility:%d\n",distenceList[0],i);
			if(eType == wkbPoint)
			{
				//~ printf("wkbPoint\n");
				if(graphIndex<1)
					totalCountExact++;
				if(distenceList[0]==totalCountExact)
					distenceList[0]=0;
				OGRPoint* pOGRPoint=(OGRPoint*) poGeometry;
				double px = pOGRPoint->getX();
				double py = pOGRPoint->getY();
				
				box pointBuffer;
				bg::buffer(bg::return_envelope<box>(point(px,py)),pointBuffer,maxBound);
				std::vector<valueS> segmentResult;
				(*rtreesubgraphEdge).query(bgi::intersects(pointBuffer)&&bgi::nearest(point(px,py),1), std::back_inserter(segmentResult));
				if(segmentResult.size()>0)
				{
					
					//~ printf("-----%ld\n",segmentResult.size());
					segment pSegment= boost::get<0>(segmentResult.front());
					double x0= bg::get<0,0>(pSegment);
					double y0= bg::get<0,1>(pSegment);
					double x1= bg::get<1,0>(pSegment);
					double y1= bg::get<1,1>(pSegment);
					int edgeIndex= boost::get<1>(segmentResult.front());
					int segmentIndex=boost::get<2>(segmentResult.front());
					//~ printf("-----%d %d\n",edgeIndex,segmentIndex);
					//~ printf("%lf   %lf \n",distanceP2P(x0,y0,x1,y1),edgelengthList[edgeIndex][segmentIndex-1]);
					int headId = edgeheadList[edgeIndex];
					int endId = edgeendList[edgeIndex];
					
					//caculate headDistance, endDistance;
					double headDistance, endDistance;
					double segFDistance;
					int segFIndex;
					bool s=getDistance(x0,y0,x1,y1,px,py,edgelengthList[edgeIndex],segmentIndex,edgelenList[edgeIndex],rate,maxBound,
					segFDistance,segFIndex,headDistance,endDistance);
					if(s)
					{
						fdListBegin[5*fdList[0]]=edgeIndex;
						fdListBegin[5*fdList[0]+1]=segFIndex;
						fdListBegin[5*fdList[0]+2]=distenceList[0]+myId*totalCount;
						*((double *)(&fdListBegin[5*fdList[0]+3]))=segFDistance;
						fdList[0]++;
						
						double* distenceheadList=(double*)malloc(vertexIndex*sizeof(double));
						double* distenceendList=(double*)malloc(vertexIndex*sizeof(double));
						shortestDistance(vertexIndex,headId, (double **)adjacencydistList,(int **)adjacencyList,neighborcList,distenceheadList);
						shortestDistance(vertexIndex,endId, (double **)adjacencydistList,(int **)adjacencyList,neighborcList,distenceendList);
						
						double* currentBegin=(double*)(distenceListBegin+distenceList[0]);
						for(int j=0;j<vertexIndex;j++)
						{
							double hd=distenceheadList[j]+headDistance;
							double ed=distenceendList[j]+endDistance;
							double mind= hd<ed? hd:ed;
							if(maxBound>mind)
								currentBegin[j*totalCount]=mind;
								
						}
					}
					
					
				}
				//~ printf("distenceList[0]:%d facility:%d\n",distenceList[0],i);
				distenceList[0]++;
				if(distenceList[0]>totalCountExact)
				distenceList[0]=distenceList[0]%totalCountExact;
				//~ printf("totalCountExact %d distenceList[0]:%d\n",totalCountExact,distenceList[0]);
			}

		}
	}

	MPI_Allgather(distenceList,distenceListsize,MPI_BYTE,
	distenceListFinal,distenceListsize,MPI_BYTE,MPI_COMM_WORLD);

	umultimap btreeFlist;
	int beginGap=sizeof(int)+vertexIndex*totalCount*sizeof(double);
	for(int i=0;i<numProcs;i++)
	{
		int * beginListI=(int *)((char *)distenceListFinal + distenceListsize*i+beginGap);
		int lenI=beginListI[0];
		beginListI=beginListI+1;
		
		for(int j=0;j<lenI;j++)

			btreeFlist.insert(std::make_pair(beginListI[j*5],beginListI+j*5));
	}
	
	int widthOut,heightOut;
	minXOut-=maxBound;
	minYOut-=maxBound;
	maxXout+=maxBound;
	maxYOut+=maxBound;
	widthOut=(maxXout-minXOut)/resolution;
	heightOut=(maxYOut-minYOut)/resolution;
	//~ printf("%d %d \n",widthOut,heightOut);
	
	int heightP =heightOut/numProcs+1;
	unsigned short* rasterOutFinal=(unsigned short*)malloc(widthOut*heightP*numProcs*sizeof(unsigned short));
	unsigned short* rasterP=(unsigned short*)malloc(widthOut*heightP*sizeof(unsigned short));
	memset(rasterP,0,widthOut*heightP*sizeof(unsigned short));
	
	for(int itmp=heightP*myId;(itmp<heightP*(myId+1))&&(itmp<heightOut);itmp++)
	{
		int i=itmp%heightP;
		double rpx,rpy;
		rpy=maxYOut-itmp*resolution;
		for(int j =0;j<widthOut;j++)
		{
			rpx=minXOut+j*resolution;
			bool * facilityFlag=(bool *)malloc(totalCount*numProcs);
			memset(facilityFlag,true,totalCount*numProcs);	
			for(int graphIndex=0;graphIndex<graphCount;graphIndex++)
			{
				//~ printf("graphIndex:%d \n",graphIndex);
				bgi::rtree< valueS, bgi::quadratic<MAX_NODE_SIZE> >* rtreesubgraphEdge;
				rtreesubgraphEdge=&rtreeList[graphIndex];

				box pointBuffer;
				bg::buffer(bg::return_envelope<box>(point(rpx,rpy)),pointBuffer,maxBound);
				std::vector<valueS> segmentResult;
				(*rtreesubgraphEdge).query(bgi::intersects(pointBuffer)&&bgi::nearest(point(rpx,rpy),1), std::back_inserter(segmentResult));
				if(segmentResult.size()>0&&(bg::distance(point(rpx,rpy),boost::get<0>(segmentResult.front()))<maxBound))
				{
					segment pSegment= boost::get<0>(segmentResult.front());
					double x0= bg::get<0,0>(pSegment);
					double y0= bg::get<0,1>(pSegment);
					double x1= bg::get<1,0>(pSegment);
					double y1= bg::get<1,1>(pSegment);
					int edgeIndex= boost::get<1>(segmentResult.front());
					int segmentIndex=boost::get<2>(segmentResult.front());
					int headId = edgeheadList[edgeIndex];
					int endId = edgeendList[edgeIndex];
					
					double headDistance, endDistance;
					double segFDistance;
					int segFIndex;
					bool s=getDistance(x0,y0,x1,y1,rpx,rpy,edgelengthList[edgeIndex],segmentIndex,edgelenList[edgeIndex],rate,maxBound,
					segFDistance,segFIndex,headDistance,endDistance);					
					if(s)
					{	
						if(headDistance<maxBound||endDistance<maxBound)
						{
							for(int k=0;k<numProcs;k++)
							{
								int * distenceListtmp=(int *)(distenceListFinal)+distenceListsize/4*k;
								double *distenceListHead=(double*)((double *)(distenceListtmp+1)+headId*totalCount);
								double *distenceListEnd=(double*)((double *)(distenceListtmp+1)+endId*totalCount);
								for(int kk=0;kk<distenceListtmp[0];kk++)
								{
									if(facilityFlag[k*totalCount+kk]&&
									(headDistance+distenceListHead[kk]<maxBound||endDistance+distenceListEnd[kk]<maxBound))
									{	
										rasterP[i*widthOut+j]++;
										facilityFlag[k*totalCount+kk]=false;
									}
	
								}
								
							}
						}								
					}
									
					auto rangetmp=btreeFlist.equal_range(edgeIndex);
					for(auto k=rangetmp.first;k!=rangetmp.second;k++)
					{
						int* kbegin = (*k).second;
						double distmp;
						if(facilityFlag[kbegin[2]]&&
						iscommonLinenear(segFIndex,segFDistance,kbegin[1],*((double *)&kbegin[3]),
						maxBound,edgelengthList[edgeIndex],distmp))
						{
							rasterP[i*widthOut+j]++;
							facilityFlag[kbegin[2]]=false;		
						}
					}
				}
		
			}
			free(facilityFlag);	
		}
	}

	MPI_Gather(rasterP,widthOut*heightP,MPI_UNSIGNED_SHORT,
	rasterOutFinal,widthOut*heightP,MPI_UNSIGNED_SHORT,0,MPI_COMM_WORLD);

	if(myId==0)
	{
		
		double adfGeoTransform[6];  
		adfGeoTransform[0] = minXOut;  
		adfGeoTransform[1] = resolution;  
		adfGeoTransform[2] = 0;  
		adfGeoTransform[3] = maxYOut;  
		adfGeoTransform[4] = 0;  
		adfGeoTransform[5] = -resolution;                
		
		const char *pszFormat = "GTiff"; 
		GDALDriver *poDriver = NULL;  
		poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat); 
		GDALDataset *poNewDS = poDriver->Create(rasterFile,widthOut,heightOut,1,GDT_UInt16,NULL); 
		
		GDALSetGeoTransform( poNewDS, adfGeoTransform );
		if (pOgrSRS != NULL)
		{
			char *pPrj = NULL;
			pOgrSRS->exportToWkt(&pPrj);
			poNewDS->SetProjection(pPrj);
		} 
		GDALRasterBand *outBand = poNewDS->GetRasterBand(1);
		outBand->RasterIO( GF_Write, 0, 0, widthOut, heightOut, rasterOutFinal, widthOut, heightOut, GDT_UInt16, 0, 0);
		//~ outBand->SetNoDataValue(0);
		GDALClose(poNewDS);
		double t5 = MPI_Wtime();
		printf("ourputsize: width: %d height: %d \n",widthOut, heightOut);
		printf("[RESULT] Successfully. cores:%d cost:<total>%f</total> s\n",numProcs, t5-t1);
	}
	
	GDALClose(facilityDS);
	MPI_Finalize();
	
}
