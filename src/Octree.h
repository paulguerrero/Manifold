#ifndef OCTREE_H_
#define OCTREE_H_

#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include <vector>
#include <list>
#include <set>
#include <map>
#include <cmath>
#include "Intersection.h"

using namespace std;

// int int_pow(int x, int p)
// {
//   if (p == 0) return 1;
//   if (p == 1) return x;

//   int tmp = myPow(x, p/2);
//   if (p%2 == 0) return tmp * tmp;
//   else return x * tmp * tmp;
// }

class Grid_Index
{
public:
    Grid_Index(){}
    Grid_Index(int x, int y, int z)
    : id(x,y,z)
    {}
    Grid_Index(const glm::ivec3& id)
    : id(id)
    {}
    bool operator<(const Grid_Index& ind) const
    {
        int i = 0;
        while (i < 3 && id[i] == ind.id[i])
            i++;
        return (i < 3 && id[i] < ind.id[i]);
    }
    Grid_Index operator+(const Grid_Index& ind) const
    {
        Grid_Index grid(*this);
        grid.id += ind.id;
        return grid;
    }
    Grid_Index operator/(int x) const
    {
        return Grid_Index(id[0]/x,id[1]/x,id[2]/x);
    }
    glm::ivec3 id;
};

// indices of child octants:
//   3    7    z  y
// 1 2  5 6    |/
// 0    4      --x
// ----------------------
// ind: (x,y,z)
// 0: (0,0,0)
// 1: (0,0,1)
// 2: (0,1,0)
// 3: (0,1,1)
// 4: (0,1,1)
// 5: (1,0,0)
// 6: (1,0,1)
// 7: (1,1,0)
// 8: (1,1,1)
//
// indices of octant sides:
// 0: x, 1: y, 2: z, 3: -x, 4: -y, 5: -z
class Octree
{
public:
    Octree()
    : children(8, NULL), empty_neighbors_side(6), depth(0), level(0), number(1), occupied(true), exterior(false), interior(false)
    {}
    Octree(glm::dvec3& min_c, glm::dvec3& max_c, vector<glm::ivec3>& faces, float thickness)
    : children(8, NULL), empty_neighbors_side(6), depth(0), level(0), number(1), occupied(true), exterior(false), interior(false)
    {
        min_corner = min_c;
        length = max_c - min_corner;
        int ind = 0;
        for (int i = 1; i < 3; ++i)
            if (length[i] > length[ind])
                ind = i;
        for (int i = 0; i < 3; ++i)
        {
            min_corner[i] -= (length[ind] - length[i]) * 0.5 + thickness * 0.5;
        }
        length = glm::dvec3(length[ind]+thickness, length[ind]+thickness, length[ind]+thickness);
        this->faces = faces;
        face_ind.resize(this->faces.size());
        for (int i = 0; i < (int)this->faces.size(); ++i)
            face_ind[i] = i;
    }

    Octree(glm::dvec3& min_c, glm::dvec3& length_)
    : children(8, NULL), empty_neighbors_side(6), depth(0), level(0), number(1), occupied(true), exterior(false), interior(false)
    {
        min_corner = min_c;
        length = length_;
    }

    ~Octree()
    {
        for (int i = 0; i < 8; ++i)
        {
            if (children[i])
                delete children[i];
        }
    }

    bool Is_Exterior(const glm::dvec3 &p)
    {
        // check if point is outside the octant
        for (int i = 0; i < 3; ++i)
            if (p[i] < min_corner[i] || p[i] > min_corner[i] + length[i])
                return true;

        // empty leaf nodes
        if (!occupied)
            return exterior;

        // occupied leaf nodes
        if (depth == 0)
            return false;

        // recurse to children
        int index = 0;
        for (int i = 0; i < 3; ++i)
        {
            index *= 2;
            if (p[i] > min_corner[i] + length[i] / 2)
                index += 1;
        }
        return children[index]->Is_Exterior(p);
    }

    bool Intersection(int face_index, glm::dvec3& min_corner, glm::dvec3& length, vector<glm::dvec3>& vertices)
    {
        float boxcenter[3];
        float boxhalfsize[3];
        float triverts[3][3];
        for (int i = 0; i < 3; ++i)
        {
            boxhalfsize[i] = length[i] * 0.5;
            boxcenter[i] = min_corner[i] + boxhalfsize[i];
        }
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                triverts[i][j] = vertices[faces[face_index][i]][j];
            }
        }
        return triBoxOverlap(boxcenter, boxhalfsize, triverts);
    }

    float min_hole_fit_size()
    {
        // the second smallest side of the bounding box gives the size of the smallest
        // hole this bounding box can be fitted through
        float min_length_xy = std::min(length[0], length[1]);
        float max_length_xy = std::max(length[0], length[1]);
        if (length[2] <= min_length_xy) {
            return min_length_xy;
        } else if (length[2] >= max_length_xy) {
            return max_length_xy;
        } else {
            return length[2];
        }
    }
    float max_length()
    {
        return std::max(length[0], std::max(length[1], length[2]));
    }

    void Split(vector<glm::dvec3>& vertices)
    {
        depth += 1;
        number = 0;
        if (depth > 1) {
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    for (int k = 0; k < 2; ++k) {
                        int ind = i * 4 + j * 2 + k;
                        if (children[ind] && children[ind]->occupied) {
                            children[ind]->Split(vertices);
                            number += children[ind]->number;
                        }
                    }
                }
            }
            faces.clear();
            face_ind.clear();
            return;
        }
        glm::dvec3 halfsize = length * glm::dvec3(0.5, 0.5, 0.5);
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 2; ++k) {
                    int ind = i * 4 + j * 2 + k;

                    glm::dvec3 startpoint = min_corner;
                    startpoint[0] += i * halfsize[0];
                    startpoint[1] += j * halfsize[1];
                    startpoint[2] += k * halfsize[2];

                    children[ind] = new Octree(startpoint, halfsize);
                    children[ind]->level = level+1;
                    children[ind]->occupied = false;
                    children[ind]->number = 0;

                    for (int face = 0; face < (int)faces.size(); ++face) {
                        if (Intersection(face, startpoint, halfsize, vertices)) {
                            children[ind]->faces.push_back(faces[face]);
                            children[ind]->face_ind.push_back(face_ind[face]);
                            if (!children[ind]->occupied) {
                                children[ind]->occupied = true;
                                number += 1;
                                children[ind]->number = 1;
                            }
                        }
                    }
                }
            }
        }
        faces.clear();
        face_ind.clear();
    }

    void FindEmptyNeighbors()
    {
        if (depth == 0)
            return;

        for (int i = 0; i < 8; ++i)
        {
            if (children[i]->occupied)
            {
                children[i]->FindEmptyNeighbors();
            }
        }
        // connect children of same parent
        int y_index[] = {0, 1, 4, 5};
        for (int i = 0; i < 4; ++i)
        {
            // dim=2 (z-axis)
            if (children[i * 2] && children[i * 2 + 1])
                ConnectEmptyNeighbors(children[i * 2], children[i * 2 + 1], 2);
            // dim=1 (y-axis)
            if (children[y_index[i]] && children[y_index[i] + 2])
                ConnectEmptyNeighbors(children[y_index[i]], children[y_index[i] + 2], 1);
            // dim=0 (x-axis)
            if (children[i] && children[i + 4])
                ConnectEmptyNeighbors(children[i], children[i + 4], 0);
        }
    }

    void ConnectEmptyNeighbors(Octree* l, Octree* r, int dim)
    {
        int y_index[] = {0, 1, 4, 5};
        if (l->occupied && r->occupied) // both are occupied (possibly non-leafs), connect children of different parents (l and r)
        {
            if (l->depth == 0)
                return;
            if (dim == 2) // z-axis
            {
                for (int i = 0; i < 4; ++i) {
                    ConnectEmptyNeighbors(l->children[i * 2 + 1], r->children[i * 2], dim);
                }
            } else
            if (dim == 1) // y-axis
            {
                for (int i = 0; i < 4; ++i) {
                    ConnectEmptyNeighbors(l->children[y_index[i] + 2], r->children[y_index[i]], dim);
                }
            } else
            if (dim == 0) // x-axis
            {
                for (int i = 0; i < 4; ++i) {
                    ConnectEmptyNeighbors(l->children[i + 4], r->children[i], dim);
                }
            }
        } else {
            if (!l->occupied && !r->occupied) // both are empty leafs
            {
                l->empty_neighbors.push_back(r);
                r->empty_neighbors.push_back(l);
            }
            if (!l->occupied) // left is an empty leaf
            {
                if (dim == 2)
                {
                    r->empty_neighbors_side[5].push_back(l);
                    if (r->depth > 0)
                    {
                        for (int i = 0; i < 4; ++i)
                        {
                            ConnectEmptyNeighbors(l, r->children[i * 2], dim);
                        }
                    }
                } else if (dim == 1)
                {
                    r->empty_neighbors_side[4].push_back(l);
                    if (r->depth > 0)
                    {
                        for (int i = 0; i < 4; ++i)
                        {
                            ConnectEmptyNeighbors(l, r->children[y_index[i]], dim);
                        }
                    }
                } else if (dim == 0)
                {
                    r->empty_neighbors_side[3].push_back(l);
                    if (r->depth > 0)
                    {
                        for (int i = 0; i < 4; ++i)
                        {
                            ConnectEmptyNeighbors(l, r->children[i], dim);
                        }
                    }
                }
            }
            if (!r->occupied) // right is an empty leaf
            {
                if (dim == 2)
                {
                    l->empty_neighbors_side[2].push_back(r);
                    if (l->depth > 0)
                    {
                        for (int i = 0; i < 4; ++i)
                        {
                            ConnectEmptyNeighbors(l->children[i * 2 + 1], r, dim);
                        }
                    }
                } else if (dim == 1)
                {
                    l->empty_neighbors_side[1].push_back(r);
                    if (l->depth > 0)
                    {
                        for (int i = 0; i < 4; ++i)
                        {
                            ConnectEmptyNeighbors(l->children[y_index[i] + 2], r, dim);
                        }
                    }
                } else if (dim == 0)
                {
                    l->empty_neighbors_side[0].push_back(r);
                    if (l->depth > 0)
                    {
                        for (int i = 0; i < 4; ++i)
                        {
                            ConnectEmptyNeighbors(l->children[i + 4], r, dim);
                        }
                    }
                }
            }
        }
    }

    void ExteriorSeeds(list<Octree*>& exterior_list, int dim)
    {
        if (!occupied)
        {
            if (!exterior)
            {
                exterior = true;
                exterior_list.push_back(this);
            }
            return;
        }
        if (depth == 0)
            return;
        int y_index[] = {0, 1, 4, 5};
        if (dim == 2 || dim == 5)
        {
            for (int i = 0; i < 4; ++i)
            {
                children[i * 2 + (dim == 5)]->ExteriorSeeds(exterior_list, dim);
            }
            return;
        }
        if (dim == 1 || dim == 4)
        {
            for (int i = 0; i < 4; ++i)
            {
                children[y_index[i] + 2 * (dim == 4)]->ExteriorSeeds(exterior_list, dim);
            }
            return;
        }
        for (int i = 0; i < 4; ++i)
        {
            children[i + 4 * (dim == 3)]->ExteriorSeeds(exterior_list, dim);
        }
    }

    void InteriorSeeds(list<Octree*>& interior_list, float min_hole_size)
    {
        if (exterior || this->min_hole_fit_size() < min_hole_size)
            return;

        if (!occupied)
        {
            if (!interior)
            {
                interior = true;
                interior_list.push_back(this);
            }
            return;
        }
        for (int i = 0; i < 8; ++i) {
            if (children[i]) {
                children[i]->InteriorSeeds(interior_list, min_hole_size);
            }
        }
    }

    void OccupiedAndEmptyLeafs(list<Octree*>& occupied_leafs, list<Octree*>& empty_leafs)
    {
        if (!occupied)
        {
            empty_leafs.push_back(this);
            return;
        }
        if (depth == 0)
        {
            occupied_leafs.push_back(this);
            return;
        }
        for (int i = 0; i < 8; ++i) {
            children[i]->OccupiedAndEmptyLeafs(occupied_leafs, empty_leafs);
        }
    }

    void ConstructOctreeMesh(
        map<Grid_Index,int>& grid_vind, const glm::ivec3& start,
        vector<glm::dvec3>& overtices, vector<glm::ivec4>& ofaces, vector<set<int> >& vertices_orig_faces,
        const glm::dvec3& face_length, int face_level)
    {

        if ((occupied && depth == 0) || interior || exterior)
        {
            std::vector<std::vector<glm::ivec3>> offsets;
            // faces are cloclwise when viewed from front (in right-handed coordinates)
            if (exterior) {
                // vertex offsets for flipped faces (reversed ordering)
                offsets = {
                    {glm::ivec3(1,0,0),glm::ivec3(1,1,0),glm::ivec3(1,1,1),glm::ivec3(1,0,1)},
                    {glm::ivec3(0,1,0),glm::ivec3(0,1,1),glm::ivec3(1,1,1),glm::ivec3(1,1,0)},
                    {glm::ivec3(0,0,1),glm::ivec3(1,0,1),glm::ivec3(1,1,1),glm::ivec3(0,1,1)},
                    {glm::ivec3(0,0,0),glm::ivec3(0,0,1),glm::ivec3(0,1,1),glm::ivec3(0,1,0)},
                    {glm::ivec3(0,0,0),glm::ivec3(1,0,0),glm::ivec3(1,0,1),glm::ivec3(0,0,1)},
                    {glm::ivec3(0,0,0),glm::ivec3(0,1,0),glm::ivec3(1,1,0),glm::ivec3(1,0,0)}};
            } else {
                // regular vertex offsets
                offsets = {
                    {glm::ivec3(1,0,0),glm::ivec3(1,0,1),glm::ivec3(1,1,1),glm::ivec3(1,1,0)},
                    {glm::ivec3(0,1,0),glm::ivec3(1,1,0),glm::ivec3(1,1,1),glm::ivec3(0,1,1)},
                    {glm::ivec3(0,0,1),glm::ivec3(0,1,1),glm::ivec3(1,1,1),glm::ivec3(1,0,1)},
                    {glm::ivec3(0,0,0),glm::ivec3(0,1,0),glm::ivec3(0,1,1),glm::ivec3(0,0,1)},
                    {glm::ivec3(0,0,0),glm::ivec3(0,0,1),glm::ivec3(1,0,1),glm::ivec3(1,0,0)},
                    {glm::ivec3(0,0,0),glm::ivec3(1,0,0),glm::ivec3(1,1,0),glm::ivec3(0,1,0)}};
            }

            int step_count = std::exp2(face_level - level);
            step_count = std::max(1, step_count);
            glm::ivec3 scaled_start = start * step_count;

            for (int side_i=0; side_i<6; ++side_i)
            {
                for (auto neighb_it = empty_neighbors_side[side_i].begin();
                        neighb_it != empty_neighbors_side[side_i].end(); ++neighb_it) {

                    // for a direct interior/exterior boundary, create the face
                    // only in the smaller of the two adjacent nodes
                    if ((occupied && (*neighb_it)->exterior) ||
                        (interior && (*neighb_it)->exterior && level >= (*neighb_it)->level) ||
                        (exterior && (*neighb_it)->interior && level > (*neighb_it)->level))
                    {
                        glm::ivec3 side_start_offset = offsets[side_i][0] * step_count;
                        glm::ivec3 side_end_offset = offsets[side_i][2] * step_count;

                        glm::ivec3 face_offset = side_start_offset;
                        for (face_offset[0]=side_start_offset[0]; face_offset[0]<max(side_end_offset[0], side_start_offset[0]+1); ++face_offset[0]) {
                            for (face_offset[1]=side_start_offset[1]; face_offset[1]<max(side_end_offset[1], side_start_offset[1]+1); ++face_offset[1]) {
                                for (face_offset[2]=side_start_offset[2]; face_offset[2]<max(side_end_offset[2], side_start_offset[2]+1); ++face_offset[2]) {

                                    glm::ivec4 oface;
                                    for (int vi = 0; vi < 4; ++vi)
                                    {
                                        glm::ivec3 vert_offset = offsets[side_i][vi] - offsets[side_i][0];

                                        Grid_Index grid_ind((scaled_start + face_offset + vert_offset) * 2);
                                        map<Grid_Index,int>::iterator it = grid_vind.find(grid_ind);

                                        if (it == grid_vind.end()) {
                                            // grid point has not been used as vertex yet

                                            glm::dvec3 vertex = min_corner;
                                            for (int dim = 0; dim < 3; ++dim)
                                                vertex[dim] += static_cast<float>(face_offset[dim]+vert_offset[dim]) * face_length[dim];

                                            overtices.push_back(vertex);
                                            vertices_orig_faces.push_back(set<int>());

                                            oface[vi] = overtices.size()-1;

                                            grid_vind.insert(make_pair(grid_ind, oface[vi]));

                                        } else {
                                            // grid point has already been used as vertex in another face
                                            oface[vi] = it->second;
                                        }

                                        // list of original faces corresponding to the octree vertex
                                        // (for projection later on), = all original faces in the octant
                                        for (vector<int>::iterator orig_face_ind_it = face_ind.begin();
                                            orig_face_ind_it != face_ind.end(); ++orig_face_ind_it)
                                            vertices_orig_faces[oface[vi]].insert(*orig_face_ind_it);
                                    }
                                    ofaces.push_back(oface);
                                }
                            }
                        }
                    }
                }
            }
        } else {
            for (int i = 0; i < 8; ++i)
            {
                if (children[i])
                {
                    int x = i / 4;
                    int y = (i - x * 4) / 2;
                    int z = i - x * 4 - y * 2;
                    glm::ivec3 child_start = start * 2 + glm::ivec3(x,y,z);
                    children[i]->ConstructOctreeMesh(
                        grid_vind, child_start, overtices, ofaces, vertices_orig_faces, face_length, face_level);
                }
            }
        }
    }

    std::vector<Octree*> children;
    list<Octree*> empty_neighbors;
    vector<list<Octree*>> empty_neighbors_side;
    int depth; // octree depth (depth of subtree starting at this node)
    int level; // octree level (graph distance to root)
    int number;
    bool occupied;
    bool exterior;
    bool interior;

    glm::dvec3 min_corner, length;

    vector<glm::ivec3> faces;
    vector<int> face_ind;
};

#endif
