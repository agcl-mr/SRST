#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/Convex_hull_d_traits_3.h>
#include <CGAL/Convex_hull_d_to_polyhedron_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/centroid.h>
#include <fstream>

#include <stdlib.h>
#include <queue>
#include <iostream>
#include <sys/time.h>
#include <boost/bind.hpp>
#include <nanoflann.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K> Vb;
typedef CGAL::Triangulation_data_structure_3<Vb> Tds;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay;
typedef CGAL::Triangle_3<K> Triangle;
typedef CGAL::Point_3<K> Point_3;
typedef CGAL::Tetrahedron_3<K> Tetrahedron;
typedef Delaunay::Point Point;
typedef Delaunay::Edge Edge;
typedef Delaunay::Facet Facet;
typedef CGAL::Segment_3<K> Segment;
typedef CGAL::Convex_hull_d_traits_3<K> Hull_traits_3;
typedef CGAL::Convex_hull_d< Hull_traits_3 > Convex_hull_3;
typedef std::list<Delaunay::Cell_handle> CellList;
typedef unsigned long long timestamp_t;
typedef std::list<int> FaceNoList;
typedef FaceNoList::iterator FaceNoIterator;

using namespace nanoflann;

Delaunay T;


/*
The struct PointCloud is used to create the kd-tree for the point cloud.
A point cloud of the centroid of the cells(tetrahedrons) is created.
Using the centroid the index of the neighboring cell is obtained.
*/
struct PointCloud
{
  struct Point
  {
    float  x, y, z;
  };

  std::vector<Point>  pts;

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const { return pts.size(); }

  inline float kdtree_distance(const float *p1, const size_t idx_p2, size_t /*size*/) const
  {
    const float d0 = p1[0] - pts[idx_p2].x;
    const float d1 = p1[1] - pts[idx_p2].y;
    const float d2 = p1[2] - pts[idx_p2].z;
    return d0*d0 + d1*d1 + d2*d2;
  }

  inline float kdtree_get_pt(const size_t idx, int dim) const
  {
    if (dim == 0) return pts[idx].x;
    else if (dim == 1) return pts[idx].y;
    else return pts[idx].z;
  }

  template <class BBOX>
  bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }

};

typedef KDTreeSingleIndexAdaptor< L2_Simple_Adaptor< float, PointCloud >, PointCloud, 3 >  KDTree;

/*
get_timestamp ()
The function is used to get the current time of the day.
This in tern is used to measure the run time of the program
*/
static timestamp_t get_timestamp ()
{
  struct timeval now;
  gettimeofday (&now, NULL);
  return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

/*
TetrahedronCell class corresponds to each tetrahedron in the Delaunay Triangulation.
Cell_handle is the cell handle for that partiualar cell.
centroid is the centroid point of that cell.
triangle[] has the list of triangles in that cell.
visited[] is used to mark whether that particular triangle is visited or not.
deleted[] is used to mark whether that particular triangle is deleted or not.
neighborTriangleIndex[] consists of the index of that triangle with respect the neighboring cell.
neighborCellIndex[] consists of the index of the neighboring cell.
neighbor[] has the index of the neighbor triangles in both the tetrahedrons(current and neighbor cell).
*/
class TetrahedronCell
{
	public:
		Delaunay::Cell_handle cell_handle;
		Point_3 centroid;

		Triangle triangle[4];
		bool retained[4];
		int faceNo[4];
	  float triangleCircumradius[4];
		bool visitedTriangle[4];
		int neighborTriangleIndex[4];
		int neighborCellIndex[4];
		int neighbor[6];
};

/*
edgeIndex - unique id for each edgeIndex
faceNoList - list of face id of the faces having that edge
*/
class EdgeIdentifier
{
	public:
  	int edgeIndex;
  	FaceNoList faceNoList;
};

/*
index - unique index for each faceNoList
circumRadius - circumRadius of the facet
facet - facet
deleted - whether the facet is deleted or not
cellId - index of the cell to which the facet belongs to
triangleId - index of the triangle in the cell
*/
class FaceId
{
	public:
	  FaceId()
	  {
	    deleted = false;
	  }
	  FaceId(int i, float cr, Facet f)
	  {
	    index = i;
	    circumRadius = cr;
	    facet = f;
	    deleted = false;
	  }
	  int index;
	  float circumRadius;
	  Facet facet;
	  bool deleted;
	  int cellId;
	  int triangleId;
};

struct FaceNode
{
  int index;
  float circumRadius;
  Point_3 centroid;
  Facet facet;
};

struct CompareFace
{
  bool operator()(const FaceNode* lhs, const  FaceNode* rhs) const
  {
    return lhs->circumRadius < rhs->circumRadius;
  }
};

std::priority_queue < FaceNode*, std::vector<FaceNode*>, CompareFace > faceQueue;

void insertFace(int index, float circumRadius, Point_3 a, Facet f)
{
  FaceNode *p;
  p = (struct FaceNode*)malloc(sizeof(struct FaceNode));
  p->index = index;
  p->circumRadius = circumRadius;
  p->centroid = a;
  p->facet = f;
  faceQueue.push(p);
}

/*
Used to check whether a cell is an infinte cell or not
*/
bool is_infinite(Delaunay::Cell_handle ch)
{
	if (ch->vertex(0) == T.infinite_vertex() || ch->vertex(1) == T.infinite_vertex() || ch->vertex(2) == T.infinite_vertex() || ch->vertex(3) == T.infinite_vertex())
		return true;
	return false;
}

/*Computes the centroid of the cell */
Point_3 centroidOfCell(Delaunay::Cell_handle cell)
{
  Point_3 p1 = cell->vertex(0)->point();
  Point_3 p2 = cell ->vertex(1)->point();
  Point_3 p3 = cell->vertex(2)->point();
  Point_3 p4 = cell->vertex(3)->point();

  float c1 = (p1.x() + p2.x() + p3.x() + p4.x()) / 4;
  float c2 = (p1.y() + p2.y() + p3.y() + p4.y()) / 4;
  float c3 = (p1.z() + p2.z() + p3.z() + p4.z()) / 4;

  Point_3 centroid(c1, c2, c3);

  return centroid;
}

Point_3 centroidOfTriangle(Triangle t)
{
  Point_3 p1 = t.vertex(0);
  Point_3 p2 = t.vertex(1);
  Point_3 p3 = t.vertex(2);

  float c1 = (p1.x() + p2.x() + p3.x()) / 3;
  float c2 = (p1.y() + p2.y() + p3.y()) / 3;
  float c3 = (p1.z() + p2.z() + p3.z()) / 3;

  Point_3 centroid(c1, c2, c3);

  return centroid;
}

/*computes the distance between two points*/
double distance(Point a, Point b)
{
	return sqrt(((a.x() - b.x())*(a.x() - b.x())) + ((a.y() - b.y())*(a.y() - b.y())) + ((a.z() - b.z())*(a.z() - b.z())));
}

double calculateCircumradius(Triangle t)
{
  return distance(CGAL::circumcenter(t), t.vertex(0));
}

/*checks whether two triangles are equal or not */
bool triEqual(Triangle t1, Triangle t2)
{
	Point_3 c10 = t1.vertex(0);
	Point_3 c11 = t1.vertex(1);
	Point_3 c12 = t1.vertex(2);

	Point_3 c20 = t2.vertex(0);
	Point_3 c21 = t2.vertex(1);
	Point_3 c22 = t2.vertex(2);

	int number = 0;
	if (c10 == c20 || c10 == c21 || c10 == c22)
		number++;
	if (c11 == c20 || c11 == c21 || c11 == c22)
		number++;
	if (c12 == c20 || c12 == c21 || c12 == c22)
		number++;

	if (number == 3)
		return true;
	return false;
}

int getFaceid(TetrahedronCell *tetrahedronCells, int neighborCellIndex, Triangle t)
{
	for(int j = 0; j < 4; j++)
	{
		if(triEqual(tetrahedronCells[neighborCellIndex].triangle[j], t))
		{
			return j;
		}

	}
}

/*Computes the midPoint of the two points*/
Point_3 calculateMidPoint(Point_3 a, Point_3 b)
{
  float x = (a.x() + b.x())/2;
  float y = (a.y() + b.y())/2;
  float z = (a.z() + b.z())/2;

  Point_3 mp(x, y, z);
  return mp;
}

/*computes the circumRadius of a triangle formed by the points a,b,c */
double calculate_area_triangle(Point a, Point b, Point c)
{
	return distance(CGAL::circumcenter(a, b, c), a);
}

/*checks whether the cell c2 has the triangle c1 or not */
bool cellsEqual(Triangle c1, Delaunay::Cell_handle c2)
{
	Point_3 c10 = c1.vertex(0);
	Point_3 c11 = c1.vertex(1);
	Point_3 c12 = c1.vertex(2);

	Point_3 c20 = Point_3(c2->vertex(0)->point());
	Point_3 c21 = Point_3(c2->vertex(1)->point());
	Point_3 c23 = Point_3(c2->vertex(3)->point());
	Point_3 c22 = Point_3(c2->vertex(2)->point());
	int number = 0;
	if (c10 == c20 || c10 == c21 || c10 == c22 || c10 == c23)
		number++;
	if (c11 == c20 || c11 == c21 || c11 == c22 || c11 == c23)
		number++;
	if (c12 == c20 || c12 == c21 || c12 == c22 || c12 == c23)
		number++;

	if (number == 3)
		return true;
	return false;
}

/*writes the triangle to the STL file*/
void writeTriangle(Triangle t, std::ofstream &outFile)
{
  outFile << " facet normal 0.0 0.0 0.0" << std::endl;
  outFile << "  outer loop" << std::endl;

  outFile << "   vertex " << t.vertex(0) << std::endl;
  outFile << "   vertex " << t.vertex(1) << std::endl;
  outFile << "   vertex " << t.vertex(2) << std::endl;

  outFile << "  endloop" << std::endl;
  outFile << " endfacet" << std::endl;
}


int main(int argv, char** argc)
{
  float para;
  std::cout << "Enter the parameter:" ;
  std::cin >> para;
  std::cout << std::endl;


	/*Input and output files*/
  std::ifstream infile(argc[1]);
  std::ofstream outFile(argc[2]);

	/*Reading the point cloud from the input file */
  float a, b, c;
  std::cout << "Reading input..."<< std::endl;
  while (infile >> a >> b >> c)
  {
    T.insert(Point(a, b, c));
  }
  infile.close();

	/*Finding the number of cells in the Delaunay Triangulation */
  int numberOfCells = T.number_of_finite_cells();
  std::cout << "Number of Cells: " << numberOfCells << std::endl;

	/*
	num_results - Number of neighbors to get using the kNN
	ret_index - index of the points returned by kNN
	out_dist_sqr - square of the distance between the query point and the returend points
	*/
	const size_t num_results = 2;
  std::vector<size_t>   ret_index(num_results);
  std::vector<float> out_dist_sqr(num_results);

	double smallest_area_triangle = 999999999.0;
	double area_triangle;

	/*
	edgeCloud is used to store the mid point of the edges as kd-tree
	The mid point is later used to query the edge index to link the facets
	*/
	PointCloud edgeCloud;
  int numberOfEdges = T.number_of_finite_edges();
  edgeCloud.pts.resize(numberOfEdges);
  std::cout << "Number of Edges: " << numberOfEdges << std::endl;

  EdgeIdentifier* edgeIdentifier = new EdgeIdentifier[numberOfEdges];
  timestamp_t t0 = get_timestamp(); //starting the timer
  int edgeCount = 0;
  for(Delaunay::Finite_edges_iterator fei = T.finite_edges_begin(); fei != T.finite_edges_end(); fei++)
  {
    Edge edge = *fei;
    Segment seg = T.segment(edge);
    Point_3 a = seg.source();
    Point_3 b = seg.target();
    Point_3 midPoint = calculateMidPoint(a, b);
    float len = distance(a, b);

    edgeCloud.pts[edgeCount].x = midPoint.x();
    edgeCloud.pts[edgeCount].y = midPoint.y();
    edgeCloud.pts[edgeCount].z = midPoint.z();

    edgeIdentifier[edgeCount].edgeIndex = edgeCount;

    edgeCount++;
  }

  KDTree edgeIndex(3 /*dim*/, edgeCloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
  edgeIndex.buildIndex();

  std::cout << "Edges Indexed..." << std::endl;

	int faceCount = 0;
  int numberOfFacets = T.number_of_finite_facets();
  PointCloud faceCloud;
  faceCloud.pts.resize(numberOfFacets);
  std::cout << "Number of Facets: " << numberOfFacets << std::endl;

  FaceId* faceId = new FaceId[numberOfFacets];

	for(Delaunay::Finite_facets_iterator ffi = T.finite_facets_begin(); ffi != T.finite_facets_end(); ffi++)
  {
    Facet facet = *ffi;
    Triangle tri = T.triangle(facet);
    Point_3 centroid = centroidOfTriangle(tri);
    float circumRadius = calculateCircumradius(tri);

    insertFace(faceCount, circumRadius, centroid, facet);

    faceId[faceCount].index = faceCount;
    faceId[faceCount].circumRadius = circumRadius;
    faceId[faceCount].facet = facet;

    faceCloud.pts[faceCount].x = centroid.x();
    faceCloud.pts[faceCount].y = centroid.y();
    faceCloud.pts[faceCount].z = centroid.z();

    Point_3 midPoint0 = calculateMidPoint(tri.vertex(1), tri.vertex(2));
    float queryPoint1[3] = {midPoint0.x(), midPoint0.y(), midPoint0.z()};
    edgeIndex.knnSearch(&queryPoint1[0], num_results, &ret_index[0], &out_dist_sqr[0]);
    int edgeId = ret_index[0];
    edgeIdentifier[edgeId].faceNoList.push_back(faceCount);

    Point_3 midPoint1 = calculateMidPoint(tri.vertex(0), tri.vertex(2));
    float queryPoint2[3] = {midPoint1.x(), midPoint1.y(), midPoint1.z()};
    edgeIndex.knnSearch(&queryPoint2[0], num_results, &ret_index[0], &out_dist_sqr[0]);
    edgeId = ret_index[0];
    edgeIdentifier[edgeId].faceNoList.push_back(faceCount);

    Point_3 midPoint2 = calculateMidPoint(tri.vertex(0), tri.vertex(1));
    float queryPoint3[3] = {midPoint2.x(), midPoint2.y(), midPoint2.z()};
    edgeIndex.knnSearch(&queryPoint3[0], num_results, &ret_index[0], &out_dist_sqr[0]);
    edgeId = ret_index[0];
    edgeIdentifier[edgeId].faceNoList.push_back(faceCount);

    faceCount++;
  }

  KDTree faceIndex(3 /*dim*/, faceCloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
  faceIndex.buildIndex();

  std::cout << "Faces Indexed..." << std::endl;

	PointCloud cloud; //stores the centroid of the cells
  cloud.pts.resize(numberOfCells);

  TetrahedronCell* tetrahedronCells = new TetrahedronCell[numberOfCells];

  int currentCellNumber = 0;

  for(Delaunay::Finite_cells_iterator fci = T.finite_cells_begin(); fci != T.finite_cells_end(); fci++)
  {
    tetrahedronCells[currentCellNumber].cell_handle = fci;
    tetrahedronCells[currentCellNumber].centroid  = centroidOfCell(fci);

    cloud.pts[currentCellNumber].x = tetrahedronCells[currentCellNumber].centroid.x();
    cloud.pts[currentCellNumber].y = tetrahedronCells[currentCellNumber].centroid.y();
    cloud.pts[currentCellNumber].z = tetrahedronCells[currentCellNumber].centroid.z();

    tetrahedronCells[currentCellNumber].triangle[0] = Triangle(fci->vertex(1)->point(), fci->vertex(2)->point(), fci->vertex(3)->point());
    tetrahedronCells[currentCellNumber].triangle[1] = Triangle(fci->vertex(0)->point(), fci->vertex(2)->point(), fci->vertex(3)->point());
    tetrahedronCells[currentCellNumber].triangle[2] = Triangle(fci->vertex(0)->point(), fci->vertex(1)->point(), fci->vertex(3)->point());
    tetrahedronCells[currentCellNumber].triangle[3] = Triangle(fci->vertex(0)->point(), fci->vertex(1)->point(), fci->vertex(2)->point());

    Point_3 centroid0 = centroidOfTriangle(tetrahedronCells[currentCellNumber].triangle[0]);
    Point_3 centroid1 = centroidOfTriangle(tetrahedronCells[currentCellNumber].triangle[1]);
    Point_3 centroid2 = centroidOfTriangle(tetrahedronCells[currentCellNumber].triangle[2]);
    Point_3 centroid3 = centroidOfTriangle(tetrahedronCells[currentCellNumber].triangle[3]);

    float queryPoint0[3] = {centroid0.x(), centroid0.y(), centroid0.z()};
    float queryPoint1[3] = {centroid1.x(), centroid1.y(), centroid1.z()};
    float queryPoint2[3] = {centroid2.x(), centroid2.y(), centroid2.z()};
    float queryPoint3[3] = {centroid3.x(), centroid3.y(), centroid3.z()};

    faceIndex.knnSearch(&queryPoint0[0], num_results, &ret_index[0], &out_dist_sqr[0]);
    int faceNo = ret_index[0];
    tetrahedronCells[currentCellNumber].faceNo[0] = faceNo;

    faceIndex.knnSearch(&queryPoint1[0], num_results, &ret_index[0], &out_dist_sqr[0]);
    faceNo = ret_index[0];
    tetrahedronCells[currentCellNumber].faceNo[1] = faceNo;

    faceIndex.knnSearch(&queryPoint2[0], num_results, &ret_index[0], &out_dist_sqr[0]);
    faceNo = ret_index[0];
    tetrahedronCells[currentCellNumber].faceNo[2] = faceNo;

    faceIndex.knnSearch(&queryPoint3[0], num_results, &ret_index[0], &out_dist_sqr[0]);
    faceNo = ret_index[0];
    tetrahedronCells[currentCellNumber].faceNo[3] = faceNo;

    tetrahedronCells[currentCellNumber].triangleCircumradius[0] = calculateCircumradius(tetrahedronCells[currentCellNumber].triangle[0]);
    tetrahedronCells[currentCellNumber].triangleCircumradius[1] = calculateCircumradius(tetrahedronCells[currentCellNumber].triangle[1]);
    tetrahedronCells[currentCellNumber].triangleCircumradius[2] = calculateCircumradius(tetrahedronCells[currentCellNumber].triangle[2]);
    tetrahedronCells[currentCellNumber].triangleCircumradius[3] = calculateCircumradius(tetrahedronCells[currentCellNumber].triangle[3]);

    tetrahedronCells[currentCellNumber].retained[0] = true;
    tetrahedronCells[currentCellNumber].retained[1] = true;
    tetrahedronCells[currentCellNumber].retained[2] = true;
    tetrahedronCells[currentCellNumber].retained[3] = true;

    tetrahedronCells[currentCellNumber].visitedTriangle[0] = false;
    tetrahedronCells[currentCellNumber].visitedTriangle[1] = false;
    tetrahedronCells[currentCellNumber].visitedTriangle[2] = false;
    tetrahedronCells[currentCellNumber].visitedTriangle[3] = false;

    tetrahedronCells[currentCellNumber].neighborCellIndex[0] = -2;
    tetrahedronCells[currentCellNumber].neighborCellIndex[1] = -2;
    tetrahedronCells[currentCellNumber].neighborCellIndex[2] = -2;
    tetrahedronCells[currentCellNumber].neighborCellIndex[3] = -2;

    tetrahedronCells[currentCellNumber].neighborTriangleIndex[0] = -2;
    tetrahedronCells[currentCellNumber].neighborTriangleIndex[1] = -2;
    tetrahedronCells[currentCellNumber].neighborTriangleIndex[2] = -2;
    tetrahedronCells[currentCellNumber].neighborTriangleIndex[3] = -2;

		if (is_infinite(fci->neighbor(0)))
		{
			area_triangle = calculate_area_triangle(fci->vertex(1)->point(), fci->vertex(2)->point(), fci->vertex(3)->point());
			if (area_triangle<smallest_area_triangle)
			{
				smallest_area_triangle = area_triangle;
			}
		}

		if (is_infinite(fci->neighbor(1)))
		{
			area_triangle = calculate_area_triangle(fci->vertex(0)->point(), fci->vertex(2)->point(), fci->vertex(3)->point());

			if (area_triangle<smallest_area_triangle)
			{
				smallest_area_triangle = area_triangle;
			}
		}

		if (is_infinite(fci->neighbor(2)))
		{
			area_triangle = calculate_area_triangle(fci->vertex(0)->point(), fci->vertex(1)->point(), fci->vertex(3)->point());

			if (area_triangle<smallest_area_triangle)
			{
				smallest_area_triangle = area_triangle;
			}
		}

		if (is_infinite(fci->neighbor(3)))
		{
			area_triangle = calculate_area_triangle(Point(fci->vertex(1)->point()), Point(fci->vertex(2)->point()), Point(fci->vertex(0)->point()));
			if (area_triangle<smallest_area_triangle)
			{
				smallest_area_triangle = area_triangle;
			}
		}
    currentCellNumber++;
  }

  KDTree index(3 /*dim*/, cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
  index.buildIndex();

  std::cout << "KD Tree constructed.." << std::endl;

  for (int n = 0; n < numberOfCells; n++) // To visit all the cells
  {
    for (int j = 0; j < 4; j++) //To visit all the triangles of the cell
    {
      if (!tetrahedronCells[n].visitedTriangle[j]) // Check whether the triangle is already visited or not
      {
        if (!is_infinite(tetrahedronCells[n].cell_handle->neighbor(j)))  // if not infinite find the index of neighbor cell
        {
          Delaunay::Cell_handle neighborCell = tetrahedronCells[n].cell_handle->neighbor(j);  // get the neighbor of the triangle j
          Point_3 cen = centroidOfCell(neighborCell);
          float queryPoint[3] = { cen.x(), cen.y(), cen.z() };
          index.knnSearch(&queryPoint[0], num_results, &ret_index[0], &out_dist_sqr[0]);
          for (int x = 0; x < num_results; x++)
          {
            if (tetrahedronCells[ret_index[x]].cell_handle == neighborCell)
            {
              tetrahedronCells[n].visitedTriangle[j] = true;
              int neighbourIndex = ret_index[x];
              tetrahedronCells[n].neighborCellIndex[j] = ret_index[x];
              int faceid = getFaceid(tetrahedronCells, tetrahedronCells[n].neighborCellIndex[j], tetrahedronCells[n].triangle[j]);
              tetrahedronCells[n].neighborTriangleIndex[j] = faceid;
              tetrahedronCells[ret_index[x]].visitedTriangle[faceid] = true;
              tetrahedronCells[ret_index[x]].neighborCellIndex[faceid] = n;
              tetrahedronCells[ret_index[x]].neighborTriangleIndex[faceid] = j;
              int fno = tetrahedronCells[n].faceNo[j];

              faceId[fno].cellId =  n;
              faceId[fno].triangleId =  j;

              break;
            }
          }
        }
        else
        {
          tetrahedronCells[n].visitedTriangle[j] = true;
          tetrahedronCells[n].neighborCellIndex[j] = -1;
          tetrahedronCells[n].neighborTriangleIndex[j] = -1;

          faceId[tetrahedronCells[n].faceNo[j]].cellId =  n;
          faceId[tetrahedronCells[n].faceNo[j]].triangleId =  j;
        }
      }
    }
  }

  std::cout << "Neighbor Details Updated" << std::endl;

  while (!faceQueue.empty())
	{
		FaceNode *n;
		n = faceQueue.top();
		faceQueue.pop();
		int neighbor[6];

    int faceN = n->index;
    int triIndex = faceId[faceN].triangleId;
    int cIndex = faceId[faceN].cellId;
    int nTriIndex = tetrahedronCells[cIndex].neighborTriangleIndex[triIndex];
    int nCellIndex = tetrahedronCells[cIndex].neighborCellIndex[triIndex];

		if (triIndex == 0)
		{
			neighbor[0] = 1;
			neighbor[1] = 2;
			neighbor[2] = 3;
		}
		else if (triIndex == 1)
		{
			neighbor[0] = 0;
			neighbor[1] = 2;
			neighbor[2] = 3;
		}
		else if (triIndex == 2)
		{
			neighbor[0] = 0;
			neighbor[1] = 1;
			neighbor[2] = 3;
		}
		else if (triIndex == 3)
		{
			neighbor[0] = 0;
			neighbor[1] = 1;
			neighbor[2] = 2;
		}

		if (nTriIndex == 0)
		{
			neighbor[3] = 1;
			neighbor[4] = 2;
			neighbor[5] = 3;
		}
		else if (nTriIndex == 1)
		{
			neighbor[3] = 0;
			neighbor[4] = 1;
			neighbor[5] = 3;
		}
		else if (nTriIndex == 2)
		{
			neighbor[3] = 0;
			neighbor[4] = 1;
			neighbor[5] = 3;
		}
		else if (nTriIndex == 3)
		{
			neighbor[3] = 0;
			neighbor[4] = 1;
			neighbor[5] = 2;
		}

		if ((tetrahedronCells[cIndex].neighborCellIndex[triIndex] != -1 && !tetrahedronCells[cIndex].retained[neighbor[0]] && !tetrahedronCells[cIndex].retained[neighbor[1]] && !tetrahedronCells[cIndex].retained[neighbor[2]] && !tetrahedronCells[nCellIndex].retained[neighbor[3]] && !tetrahedronCells[nCellIndex].retained[neighbor[4]] && !tetrahedronCells[nCellIndex].retained[neighbor[5]]) || (tetrahedronCells[cIndex].neighborCellIndex[triIndex] == -1 && !tetrahedronCells[cIndex].retained[neighbor[0]] && !tetrahedronCells[cIndex].retained[neighbor[1]] && !tetrahedronCells[cIndex].retained[neighbor[2]]) || ( n->circumRadius < para*smallest_area_triangle) )
		{
			if (!(tetrahedronCells[cIndex].neighborCellIndex[triIndex] != -1 && !tetrahedronCells[cIndex].retained[neighbor[0]] && !tetrahedronCells[cIndex].retained[neighbor[1]] && !tetrahedronCells[cIndex].retained[neighbor[2]] && !tetrahedronCells[nCellIndex].retained[neighbor[3]] && !tetrahedronCells[nCellIndex].retained[neighbor[4]] && !tetrahedronCells[nCellIndex].retained[neighbor[5]]) && !(tetrahedronCells[cIndex].neighborCellIndex[triIndex] == -1 && !tetrahedronCells[cIndex].retained[neighbor[0]] && !tetrahedronCells[cIndex].retained[neighbor[1]] && !tetrahedronCells[cIndex].retained[neighbor[2]]))
			{
				if (tetrahedronCells[cIndex].retained[neighbor[0]] && tetrahedronCells[cIndex].retained[neighbor[1]] && tetrahedronCells[cIndex].retained[neighbor[2]])
				{
					float circumRad[4];
					circumRad[0] = distance(CGAL::circumcenter(Point(tetrahedronCells[cIndex].triangle[0].vertex(0)), Point(tetrahedronCells[cIndex].triangle[0].vertex(1)), Point(tetrahedronCells[cIndex].triangle[0].vertex(2))), Point(tetrahedronCells[cIndex].triangle[0].vertex(0)));
					circumRad[1] = distance(CGAL::circumcenter(Point(tetrahedronCells[cIndex].triangle[1].vertex(0)), Point(tetrahedronCells[cIndex].triangle[1].vertex(1)), Point(tetrahedronCells[cIndex].triangle[1].vertex(2))), Point(tetrahedronCells[cIndex].triangle[1].vertex(0)));
					circumRad[2] = distance(CGAL::circumcenter(Point(tetrahedronCells[cIndex].triangle[2].vertex(0)), Point(tetrahedronCells[cIndex].triangle[2].vertex(1)), Point(tetrahedronCells[cIndex].triangle[2].vertex(2))), Point(tetrahedronCells[cIndex].triangle[2].vertex(0)));
					circumRad[3] = distance(CGAL::circumcenter(Point(tetrahedronCells[cIndex].triangle[3].vertex(0)), Point(tetrahedronCells[cIndex].triangle[3].vertex(1)), Point(tetrahedronCells[cIndex].triangle[3].vertex(2))), Point(tetrahedronCells[cIndex].triangle[3].vertex(0)));
					float maximum = circumRad[0];
					int location = 0;
					for (int c = 1; c < 4; c++)
					{
						if (circumRad[c] > maximum)
						{
							maximum = circumRad[c];
							location = c;
						}
					}
					tetrahedronCells[cIndex].retained[location] = false;
          faceId[faceN].deleted = true;

					if(nCellIndex != -1)
					{
						int nti = -10;
						for(int p = 0; p<4; p++)
						{
							if((triEqual(tetrahedronCells[cIndex].triangle[location], tetrahedronCells[nCellIndex].triangle[p])))
							{
								nti = p;
								break;
							}
						}
						nti = tetrahedronCells[cIndex].neighborTriangleIndex[location];
						tetrahedronCells[nCellIndex].retained[nti] = false;
					}
				}

        if(nCellIndex != -1)
				if (tetrahedronCells[nCellIndex].retained[neighbor[3]] && tetrahedronCells[nCellIndex].retained[neighbor[4]] && tetrahedronCells[nCellIndex].retained[neighbor[5]])
				{
					float circumRad[4];
					circumRad[0] = distance(CGAL::circumcenter(Point(tetrahedronCells[nCellIndex].triangle[0].vertex(0)), Point(tetrahedronCells[nCellIndex].triangle[0].vertex(1)), Point(tetrahedronCells[nCellIndex].triangle[0].vertex(2))), Point(tetrahedronCells[nCellIndex].triangle[0].vertex(0)));
					circumRad[1] = distance(CGAL::circumcenter(Point(tetrahedronCells[nCellIndex].triangle[1].vertex(0)), Point(tetrahedronCells[nCellIndex].triangle[1].vertex(1)), Point(tetrahedronCells[nCellIndex].triangle[1].vertex(2))), Point(tetrahedronCells[nCellIndex].triangle[1].vertex(0)));
					circumRad[2] = distance(CGAL::circumcenter(Point(tetrahedronCells[nCellIndex].triangle[2].vertex(0)), Point(tetrahedronCells[nCellIndex].triangle[2].vertex(1)), Point(tetrahedronCells[nCellIndex].triangle[2].vertex(2))), Point(tetrahedronCells[nCellIndex].triangle[2].vertex(0)));
					circumRad[3] = distance(CGAL::circumcenter(Point(tetrahedronCells[nCellIndex].triangle[3].vertex(0)), Point(tetrahedronCells[nCellIndex].triangle[3].vertex(1)), Point(tetrahedronCells[nCellIndex].triangle[3].vertex(2))), Point(tetrahedronCells[nCellIndex].triangle[3].vertex(0)));
          float maximum = circumRad[0];
					int location = 0;
					for (int c = 1; c < 4; c++)
					{
						if (circumRad[c] > maximum)
						{
							maximum = circumRad[c];
							location = c;
						}
					}
					tetrahedronCells[nCellIndex].retained[location] = false;
          faceId[faceN].deleted = true;
					int nti = tetrahedronCells[nCellIndex].neighborTriangleIndex[location];
					int nci = tetrahedronCells[nCellIndex].neighborCellIndex[location];
					if (nci != -1)
					{
						tetrahedronCells[nci].retained[nti] = false;
					}
				}
			}
		}
		else
		{
			tetrahedronCells[cIndex].retained[triIndex] = false;
      faceId[faceN].deleted = true;
			if(nCellIndex != -1)
			{
				tetrahedronCells[nCellIndex].retained[nTriIndex] = false;
			}
		}
	}


  int count = 0;
  bool triangleDeleted = true;
  //removing hanging triangles

  while(triangleDeleted  && count < 0)
  {
    count ++;
    std::cout << count << std::endl;
    triangleDeleted = false;


    for(int nx = 0; nx < numberOfCells ; nx++)
    {

      for(int i = 0; i < 4; i++)
      {
        if((tetrahedronCells[nx].retained[i]))
        {
          int index = nx;
          int numberOfTriangles = 0;
          bool nonHangingTriangle = false;

          Triangle tempTriangle = tetrahedronCells[index].triangle[i];

          Point_3 midPoint0 = calculateMidPoint(tempTriangle.vertex(1), tempTriangle.vertex(2));
          float queryPoint1[3] = {midPoint0.x(), midPoint0.y(), midPoint0.z()};
          edgeIndex.knnSearch(&queryPoint1[0], num_results, &ret_index[0], &out_dist_sqr[0]);
          int edgeId = ret_index[0];

          for(FaceNoIterator fni = edgeIdentifier[edgeId].faceNoList.begin(); fni != edgeIdentifier[edgeId].faceNoList.end(); fni++)
          {
            int ti = *fni;

            if(tetrahedronCells[faceId[ti].cellId].retained[faceId[ti].triangleId])
            numberOfTriangles++;
          }

          if(numberOfTriangles > 1)
          {
            numberOfTriangles = 0;
            Point_3 midPoint1 = calculateMidPoint(tempTriangle.vertex(0), tempTriangle.vertex(2));
            float queryPoint2[3] = {midPoint1.x(), midPoint1.y(), midPoint1.z()};
            edgeIndex.knnSearch(&queryPoint2[0], num_results, &ret_index[0], &out_dist_sqr[0]);
            edgeId = ret_index[0];
            for(FaceNoIterator fni = edgeIdentifier[edgeId].faceNoList.begin(); fni != edgeIdentifier[edgeId].faceNoList.end(); fni++)
            {
              int ti = *fni;
              if(tetrahedronCells[faceId[ti].cellId].retained[faceId[ti].triangleId])
              numberOfTriangles++;
            }

            if(numberOfTriangles > 1)
            {
              numberOfTriangles = 0;
              Point_3 midPoint2 = calculateMidPoint(tempTriangle.vertex(0), tempTriangle.vertex(1));
              float queryPoint3[3] = {midPoint2.x(), midPoint2.y(), midPoint2.z()};
              edgeIndex.knnSearch(&queryPoint3[0], num_results, &ret_index[0], &out_dist_sqr[0]);
              edgeId = ret_index[0];
              for(FaceNoIterator fni = edgeIdentifier[edgeId].faceNoList.begin(); fni != edgeIdentifier[edgeId].faceNoList.end(); fni++)
              {
                int ti = *fni;
                if(tetrahedronCells[faceId[ti].cellId].retained[faceId[ti].triangleId])
                numberOfTriangles++;
              }
              if(numberOfTriangles > 1)
              nonHangingTriangle = true;
            }
          }
          if(!nonHangingTriangle)
          {
            tetrahedronCells[index].retained[i] = false;
            faceId[tetrahedronCells[index].faceNo[i]].deleted = true;
            if(tetrahedronCells[index].neighborCellIndex[i] != -1)
              tetrahedronCells[tetrahedronCells[index].neighborCellIndex[i]].retained[tetrahedronCells[index].neighborTriangleIndex[i]] = false;
            triangleDeleted = true;
          }
        }
      }
    }
  }


timestamp_t t1 = get_timestamp();

outFile << "solid testsphere" << std::endl;

  for (int n = 0; n < numberOfCells; n++)
  {
    for (int j = 0; j < 4; j++)
    {
      if (tetrahedronCells[n].retained[j])
      {
        writeTriangle(tetrahedronCells[n].triangle[j], outFile);
      }
    }
  }
  std::cout << "Finished Writing STL" << std::endl;
  outFile << " endsolid" << std::endl;
  outFile.close();












  double secs = (t1 - t0) / 1000000.0L;
  std::cout << "The execution time : " << secs << std::endl;
}
