#ifndef Model_H_
#define Model_H_

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "Octree.h"
#include <map>
#include <cstdlib>
#include <igl/readOBJ.h>
using namespace std;
/*************************************************************************** 
  OBJ Loading 
 ***************************************************************************/
 
class Model
{
public: 

	Model();			

  int Load(char *filename);	// Loads the model
  void Save(const char* filename, bool color);
  void SaveOBJ(const char* filename);
  void SavePLY(const char* filename, const std::vector<glm::dvec3>& v, const std::vector<glm::ivec3>& f, const std::vector<glm::dvec3>& c);

 	void Calc_Bounding_Box();
  void Process_Manifold(int resolution, float min_hole_size);
  void Build_Tree(int resolution, float min_hole_size);
  void Construct_Manifold();
  void Project_Manifold();
  glm::dvec3 Closest_Point( const glm::dvec3 *triangle, const glm::dvec3 &sourcePosition );
  glm::dvec3 Find_Closest(int i);
  void is_manifold(vector<pair<int,int> >& nonmanifold_edges, vector<pair<int,int> >& inconsistent_edges, vector<int>& outofbound_faces);
  void Split_Grid(map<Grid_Index,int>& grid_vind, vector<glm::dvec3>& overtices, vector<glm::ivec4>& ofaces, vector<set<int> >& vertices_orig_faces, vector<glm::ivec3>& triangles);

  double clamp(double d1, double l, double r)
  {
    if (d1 < l)
      return l;
    if (d1 > r)
      return l;
    return d1;
  }
  
  vector<glm::dvec3> vertices, orig_vertices;
  vector<glm::dvec3> colors;
	vector<set<int> > vertices_orig_faces;
  vector<Grid_Index > vertices_gridinds;
  
  vector<glm::ivec3> faces, orig_faces;
  vector<glm::dvec3> face_normals;
  
  Octree* tree;
  glm::dvec3 leaf_length;
  int leaf_level;
	glm::dvec3 min_corner, max_corner;
  
  char* filename;
};

#endif
