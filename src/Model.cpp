#include "Model.h"
#include <queue>
#include <math.h>
#include <algorithm>
#define ITER_NUM 20
int g_sharp = 0;
Model::Model()
{
}

int Model::Load(char* filename)
{
    using namespace Eigen;
    using namespace std;

    // Load a mesh in OBJ format
    MatrixXd V;
    MatrixXi F;
    igl::readOBJ(filename, V, F);
    this->filename = filename;

    // Make the example deterministic
    srand(0);
    vertices.resize(V.rows());
    faces.resize(F.rows());
    for (int i = 0; i < V.rows(); ++i)
        vertices[i] = glm::dvec3(V(i,0),V(i,1),V(i,2));
    for (int i = 0; i < F.rows(); ++i) {
        faces[i] = glm::ivec3(F(i,0),F(i,1),F(i,2));
    }

    return 0;
}

void Model::Calc_Bounding_Box()
{
    min_corner = glm::dvec3(1e30,1e30,1e30);
    max_corner = -min_corner;
    for (int i = 0; i < (int)vertices.size(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            if (vertices[i][j] < min_corner[j])
            {
                min_corner[j] = vertices[i][j];
            }
            if (vertices[i][j] > max_corner[j])
            {
                max_corner[j] = vertices[i][j];
            }
        }
    }
    glm::dvec3 length = max_corner - min_corner;
    min_corner -= length * 0.2;
    max_corner += length * 0.2;
}

void Model::Build_Tree(int resolution, float min_hole_size)
{
    Calc_Bounding_Box();
    tree = new Octree(min_corner, max_corner, faces, 0.01);

    float max_bbox_length = tree->max_length();
    min_hole_size = min_hole_size * max_bbox_length; // from fraction of largest bounding box side to absolute size
    
    float leaf_min_hole_fit_size = tree->min_hole_fit_size();
    leaf_length = tree->length;
    leaf_level = 0;
    int max_depth = 99;

    while (tree->number < resolution && tree->depth <= max_depth)
    {
        tree->Split(vertices);

        leaf_min_hole_fit_size *= 0.5;
        leaf_length *= 0.5;
        ++leaf_level;
    }

    tree->FindEmptyNeighbors();

    // add boundary empty (= leaf) nodes as seeds to grow the exterior
    list<Octree*> exterior_list;
    for (int i = 0; i < 6; ++i)
    {
        tree->ExteriorSeeds(exterior_list, i);
    }

    // grow exterior using empty_neighbors, but grow only in the octree depths >= hole_min_depth,
    // so that the exterior does not grow through small holes
    list<Octree*> exterior_front(exterior_list);
    while (!exterior_front.empty())
    {
        Octree* empty = exterior_front.front();
        for (list<Octree*>::iterator it = empty->empty_neighbors.begin();
            it != empty->empty_neighbors.end(); ++it)
        {
            if (!(*it)->exterior && (*it)->min_hole_fit_size() >= min_hole_size)
            // if (!(*it)->exterior)
            {
                (*it)->exterior = true;
                exterior_list.push_back(*it);
                exterior_front.push_back(*it);
            }
        }
        exterior_front.pop_front();
    }

    // add empty (= leaf) nodes that are not exterior and larger than min_hole_size
    list<Octree*> interior_list;
    tree->InteriorSeeds(interior_list, min_hole_size);

    // propagate interior and exterior front in interleaved steps
    exterior_front = exterior_list;
    list<Octree*> interior_front(interior_list);
    list<Octree*> next_exterior_front;
    list<Octree*> next_interior_front;
    while(!exterior_front.empty() || !interior_front.empty()) {

        for (list<Octree*>::iterator exterior_front_it = exterior_front.begin();
            exterior_front_it != exterior_front.end(); ++exterior_front_it)
        {
            for (list<Octree*>::iterator it = (*exterior_front_it)->empty_neighbors.begin();
                it != (*exterior_front_it)->empty_neighbors.end(); ++it)
            {
                // if ((*exterior_front_it)->exterior && (*it)->interior) {
                //     std::cout << "prop exterior interior connection" << std::endl << std::flush;
                // }
                if (!(*it)->exterior && !(*it)->interior)
                {
                    (*it)->exterior = true;
                    exterior_list.push_back(*it);
                    next_exterior_front.push_back(*it);
                }
            }
        }

        for (list<Octree*>::iterator interior_front_it = interior_front.begin();
            interior_front_it != interior_front.end(); ++interior_front_it)
        {
            for (list<Octree*>::iterator it = (*interior_front_it)->empty_neighbors.begin();
                it != (*interior_front_it)->empty_neighbors.end(); ++it)
            {
                if (!(*it)->exterior && !(*it)->interior)
                {
                    (*it)->interior = true;
                    interior_list.push_back(*it);
                    next_interior_front.push_back(*it);
                }
            }
        }

        exterior_front = next_exterior_front;
        next_exterior_front.clear();
        interior_front = next_interior_front;
        next_interior_front.clear();
    }

    list<Octree*> occupied_leafs;
    list<Octree*> empty_leafs;
    tree->OccupiedAndEmptyLeafs(occupied_leafs, empty_leafs);

    
    // set all empty leafs that have not been reached by the interoir/exterior growing to interior
    // (these should be leafs in small interoir voids)
    for (list<Octree*>::iterator empty_leaf_it = empty_leafs.begin();
        empty_leaf_it != empty_leafs.end(); ++empty_leaf_it)
    {
        if ((*empty_leaf_it)->exterior == 0 && (*empty_leaf_it)->interior == 0)
        {
            (*empty_leaf_it)->interior = 1;
            interior_list.push_back(*empty_leaf_it);
        }
    }
}

glm::dvec3 Model::Closest_Point( const glm::dvec3 *triangle, const glm::dvec3 &sourcePosition )
{
    glm::dvec3 edge0 = triangle[1] - triangle[0];
    glm::dvec3 edge1 = triangle[2] - triangle[0];
    glm::dvec3 v0 = triangle[0] - sourcePosition;

    double a = glm::dot(edge0, edge0 );
    double b = glm::dot(edge0, edge1 );
    double c = glm::dot(edge1, edge1 );
    double d = glm::dot(edge0, v0 );
    double e = glm::dot(edge1, v0 );

    double det = a*c - b*b;
    double s = b*e - c*d;
    double t = b*d - a*e;

    if ( s + t < det )
    {
        if ( s < 0.f )
        {
            if ( t < 0.f )
            {
                if ( d < 0.f )
                {
                    s = clamp( -d/a, 0.f, 1.f );
                    t = 0.f;
                }
                else
                {
                    s = 0.f;
                    t = clamp( -e/c, 0.f, 1.f );
                }
            }
            else
            {
                s = 0.f;
                t = clamp( -e/c, 0.f, 1.f );
            }
        }
        else if ( t < 0.f )
        {
            s = clamp( -d/a, 0.f, 1.f );
            t = 0.f;
        }
        else
        {
            float invDet = 1.f / det;
            s *= invDet;
            t *= invDet;
        }
    }
    else
    {
        if ( s < 0.f )
        {
            float tmp0 = b+d;
            float tmp1 = c+e;
            if ( tmp1 > tmp0 )
            {
                float numer = tmp1 - tmp0;
                float denom = a-2*b+c;
                s = clamp( numer/denom, 0.f, 1.f );
                t = 1-s;
            }
            else
            {
                t = clamp( -e/c, 0.f, 1.f );
                s = 0.f;
            }
        }
        else if ( t < 0.f )
        {
            if ( a+d > b+e )
            {
                float numer = c+e-b-d;
                float denom = a-2*b+c;
                s = clamp( numer/denom, 0.f, 1.f );
                t = 1-s;
            }
            else
            {
                s = clamp( -e/c, 0.f, 1.f );
                t = 0.f;
            }
        }
        else
        {
            float numer = c+e-b-d;
            float denom = a-2*b+c;
            s = clamp( numer/denom, 0.f, 1.f );
            t = 1.f - s;
        }
    }

    return triangle[0] + s * edge0 + t * edge1;
}

void Model::Construct_Manifold()
{
    map<Grid_Index,int> grid_vind;
    vector<glm::dvec3> overtices;
    vector<glm::ivec4> ofaces;
    vector<glm::ivec3> triangles;
    vector<glm::dvec3> temp_pts;
    vector<glm::ivec3> temp_faces;
    vector<glm::dvec3> temp_pt_colors;
    tree->ConstructOctreeMesh(
        grid_vind, glm::ivec3(0,0,0), overtices, ofaces, vertices_orig_faces, leaf_length, leaf_level);
    Split_Grid(grid_vind, overtices, ofaces, vertices_orig_faces, triangles);
    
    // remove unused vertices
    vector<int> hash_v(overtices.size(), 0);
    for (int i = 0; i < (int)triangles.size(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            hash_v[triangles[i][j]] = 1;
        }
    }
    vertices.clear();
    for (int i = 0; i < (int)hash_v.size(); ++i)
    {
        if (hash_v[i])
        {
            hash_v[i] = (int)vertices.size();
            vertices_orig_faces[vertices.size()] = vertices_orig_faces[i];
            vertices_gridinds[vertices.size()] = vertices_gridinds[i];
            vertices.push_back(overtices[i]);
            colors.push_back(glm::dvec3(1,1,1));
        }
    }
    for (int i = 0; i < (int)triangles.size(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            triangles[i][j] = hash_v[triangles[i][j]];
        }
    }
    faces = triangles;
}

glm::dvec3 Model::Find_Closest(int i)
{
    glm::dvec3 cpoint = glm::dvec3(1e20,1e20,1e20);
    glm::dvec3 tris[3];
    glm::dvec3 normal;
    for (set<int>::iterator it = vertices_orig_faces[i].begin();
        it != vertices_orig_faces[i].end(); ++it)
    {
        int face_ind = *it;
        for (int j = 0; j < 3; ++j)
            tris[j] = orig_vertices[orig_faces[face_ind][j]];
        glm::dvec3 p = Closest_Point(tris, vertices[i]);
        if (glm::length(p-vertices[i]) < glm::length(cpoint-vertices[i]))
        {
            normal = glm::normalize(glm::cross(tris[1]-tris[0],tris[2]-tris[0]));
            if (glm::dot(normal,vertices[i]-cpoint)<0)
                normal = -normal;
            cpoint = p;
        }
    }
    return cpoint + normal * 5e-4;
}

void Model::Project_Manifold()
{
    double len = glm::length(vertices[faces[0][1]] - vertices[faces[0][0]]);
    double min_len = glm::length(vertices[faces[0][2]] - vertices[faces[0][0]]);
    if (min_len < len)
        len = min_len;
    colors.clear();
    colors.resize(vertices.size(),glm::dvec3(1,1,1));

    vector<vector<int> > vertex_faces(vertices.size());
    face_normals.resize(faces.size());
    for (int i = 0; i < (int)faces.size(); ++i)
    {
        int id[3];
        id[0] = faces[i][0];
        id[1] = faces[i][1];
        id[2] = faces[i][2];
        for (int j = 0; j < 3; ++j)
        {
            vertex_faces[id[j]].push_back(i);
        }
    }
    vector<int> vertices_hash(vertices.size(), 0);
    double min_step = 2.0 / ITER_NUM;

    for (int iter = 0; iter < ITER_NUM; ++iter) {

        for (int i = 0; i < (int)faces.size(); ++i)
        {
            int id[3];
            id[0] = faces[i][0];
            id[1] = faces[i][1];
            id[2] = faces[i][2];
            face_normals[i] = glm::normalize(glm::cross(vertices[id[1]] - vertices[id[0]],
                vertices[id[2]] - vertices[id[0]]));
        }

        vector<int> invalid_vertices;
        vector<int> invalid_indices(vertices.size(), -1);

        for (int i = 0; i < (int)vertices.size(); ++i)
        {
            if (vertices_hash[i])
                continue;
            glm::dvec3 cpoint = Find_Closest(i);
            glm::dvec3 move_dir = cpoint - vertices[i];
            double orig_step = glm::length(move_dir);
            move_dir /= orig_step;
            double step = orig_step;
            bool flag = step < 1e15;
            if (g_sharp) {
                vertices[i] = cpoint;
                continue;
            }
            glm::dvec3 normal(0,0,0);
            for (int j = 0; j < (int)vertex_faces[i].size(); ++j)
            {
                normal += face_normals[vertex_faces[i][j]];
            }
            normal = glm::normalize(normal);
            if (flag) {
                bool convex = true;
                for (int j = 0; j < (int)vertex_faces[i].size(); ++j)
                {
                    for (int k = 0; k < 3; ++k) {
                        if (glm::dot(vertices[faces[vertex_faces[i][j]][k]] - vertices[i], normal) > 0)
                            convex = false;
                    }
                    if (!convex)
                        break;
                }
                if (convex)
                {
                    for (int j = 0; j < (int)vertex_faces[i].size(); ++j)
                    {
                        if (glm::dot(face_normals[vertex_faces[i][j]],move_dir)>0)
                        {
                            flag = false;
                            break;
                        }
                    }
                } else
                {
                    flag = false;
                    for (int j = 0; j < (int)vertex_faces[i].size(); ++j)
                    {
                        if (glm::dot(face_normals[vertex_faces[i][j]],move_dir)<0)
                        {
                            flag = true;
                            break;
                        }
                    }
                }
            }
            if (flag)
            {
                if (step > min_step * len)
                    step = min_step * len;
                for (int j = 0; j < (int)vertex_faces[i].size(); ++j)
                {
                    glm::ivec3& face_index = faces[vertex_faces[i][j]];
                    int t = 0;
                    while (face_index[t] != i)
                        t += 1;
                    glm::dvec3 dir = glm::normalize(vertices[face_index[(t+2)%3]] - vertices[face_index[(t+1)%3]]);
                    glm::dvec3 h = vertices[face_index[(t+1)%3]]+glm::dot(vertices[i] - vertices[face_index[(t+1)%3]],dir) * dir - vertices[i];
                    double h_len = glm::length(h);
                    h /= h_len;
                    double h_step = glm::dot(h,move_dir) * step;
                    if (h_step > h_len * 0.7)
                    {
                        step *= (h_len * 0.7) / h_step;
                        invalid_indices[i] = (int)invalid_vertices.size();
                        invalid_vertices.push_back(i);
                        colors[i] = glm::dvec3(0,0,1);
                    }
                }
                if (fabs(step - orig_step)<1e-6) {
                    vertices[i] = cpoint + len * normal;
                    if (step > 1e-4)
                        step -= 1e-4;
                    else
                        step = 0;
                    vertices_hash[i] = 1;
                } else
                vertices[i] += step * move_dir;
                for (int j = 0; j < (int)vertex_faces[i].size(); ++j)
                {
                    int face = vertex_faces[i][j];
                    face_normals[face] = glm::normalize(glm::cross(vertices[faces[face][1]] - vertices[faces[face][0]],
                        vertices[faces[face][2]] - vertices[faces[face][0]]));
                }
            } else
            {
                invalid_indices[i] = (int)invalid_vertices.size();
                invalid_vertices.push_back(i);
                colors[i] = glm::dvec3(0,0,1);
            }
        }
        vector<int> invalid_colors(invalid_vertices.size(), -1);
        int c = 0;
        for (int i = 0; i < (int)invalid_vertices.size(); ++i)
        {
            if (invalid_colors[i] == -1)
            {
                invalid_colors[i] = c;
                vector<int> queue;
                int f = 0;
                queue.push_back(i);
                while (f != (int)queue.size())
                {
                    int id = invalid_vertices[queue[f]];
                    for (int j = 0; j < (int)vertex_faces[id].size(); ++j)
                    {
                        for (int k = 0; k < 3; ++k)
                        {
                            int index = invalid_indices[faces[vertex_faces[id][j]][k]];
                            if (index != -1 && invalid_colors[index] == -1)
                            {
                                invalid_colors[index] = c;
                                queue.push_back(index);
                            }
                        }
                    }
                    f++;
                }
                for (vector<int>::reverse_iterator it = queue.rbegin();
                    it != queue.rend(); ++it)
                {
                    glm::dvec3 midpoint(0,0,0);
                    int count = 0;
                    int id = invalid_vertices[*it];
                    for (int j = 0; j < (int)vertex_faces[id].size(); ++j)
                    {
                        for (int k = 0; k < 3; ++k)
                        {
                            int vind = faces[vertex_faces[id][j]][k];
                            if (invalid_indices[vind] == -1 || invalid_colors[invalid_indices[vind]] == -1)
                            {
                                midpoint += vertices[vind];
                                count += 1;
                            }
                        }
                    }
                    glm::dvec3 move_dir = midpoint / (double)count - vertices[id];
                    invalid_colors[*it] = -1;
                    double l = glm::length(move_dir);
                    if (l == 0 || count == 0)
                        continue;
                    move_dir /= l;
                    for (int j = 0; j < (int)vertex_faces[id].size(); ++j)
                    {
                        glm::ivec3& face_index = faces[vertex_faces[id][j]];
                        int t = 0;
                        while (face_index[t] != id)
                            t += 1;
                        glm::dvec3 dir = glm::normalize(vertices[face_index[(t+2)%3]] - vertices[face_index[(t+1)%3]]);
                        glm::dvec3 h = vertices[face_index[(t+1)%3]]+glm::dot(vertices[id] - vertices[face_index[(t+1)%3]],dir) * dir - vertices[id];
                        double h_len = glm::length(h);
                        h /= h_len;
                        double h_step = glm::dot(h,move_dir) * l;
                        if (h_step > h_len * 0.7)
                        {
                            l *= (h_len * 0.7) / h_step;
                        }
                    }
                    move_dir *= l;
                    vertices[id] += move_dir;
                }
                for (int i = 0; i < (int)queue.size(); ++i)
                {
                    invalid_colors[queue[i]] = c;
                }
            }
            c++;
        }
    }
}

void Model::Process_Manifold(int resolution, float min_hole_size)
{
    orig_vertices = vertices;
    orig_faces = faces;
    Build_Tree(resolution, min_hole_size);
    Construct_Manifold();
    Project_Manifold();

    // smoothing (each smoothed vertex is average of its 1-ring)
    for (int iter = 0; iter < 1; ++iter)
    {
        vector<glm::dvec3> dis(vertices.size());
        vector<int> dis_weight(vertices.size());
        for (int i = 0; i < (int)faces.size(); ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                int x = faces[i][j];
                int y = faces[i][(j + 1) % 3];
                dis[x] += vertices[y];
                dis[y] += vertices[x];
                dis_weight[x] += 1;
                dis_weight[y] += 1;
            }
        }
        for (int i = 0; i < (int)vertices.size(); ++i)
        {
            if (dis_weight[i] > 0)
                vertices[i] = dis[i] * (1.0 / dis_weight[i]);
        }
    }

    vector<pair<int,int> > nonmanifold_edges;
    vector<pair<int,int> > inconsistent_edges;
    vector<int> outofbound_faces;
    is_manifold(nonmanifold_edges, inconsistent_edges, outofbound_faces);
    if (!nonmanifold_edges.empty() || !inconsistent_edges.empty() || !outofbound_faces.empty()) {
        ofstream os("error.txt");
        os << filename << "\n";
        os.close();
        std::cout << "Not a Manifold!" << "\n";
        std::cout << "nonmanifold edge count: " << nonmanifold_edges.size() << "\n";
        std::cout << "inconsistent edge count: " << inconsistent_edges.size() << "\n";
        std::cout << "out of bound face count: " << outofbound_faces.size() << "\n";
        // exit(0);
    }
}

void Model::Split_Grid(map<Grid_Index,int>& grid_vind, vector<glm::dvec3>& overtices, vector<glm::ivec4>& ofaces, vector<set<int> >& vertices_orig_faces, vector<glm::ivec3>& triangles)
{
    double unit_len = 0;
    vertices_gridinds.resize(grid_vind.size());
    for (map<Grid_Index,int>::iterator it = grid_vind.begin();
        it != grid_vind.end(); ++it)
    {
        vertices_gridinds[it->second] = it->first;
    }
    
    // create list of unique edges with adjecent faces
    map<pair<int,int>,list<pair<int,int> > > edge_info;
    for (int i = 0; i < (int)ofaces.size(); ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            int x = ofaces[i][j];
            int y = ofaces[i][(j + 1) % 4];
            if (x > y)
            {
                int temp = x;
                x = y;
                y = temp;
            }
            pair<int,int> edge = make_pair(x,y);
            map<pair<int,int>, list<pair<int,int> > >::iterator it = edge_info.find(edge);
            if (it != edge_info.end())
            {
                it->second.push_back(make_pair(i,j));
            } else
            {
                list<pair<int,int> > buf;
                buf.push_back(make_pair(i,j));
                edge_info.insert(make_pair(edge, buf));
            }
        }
    }
    
    // mark vertices adjacent to non-manifold edges (edges with more than two adjacent faces)
    set<int> marked_v;
    for (map<pair<int,int>,list<pair<int,int> > >::iterator it = edge_info.begin();
        it != edge_info.end(); ++it)
    {
        if (it->second.size() > 2) {
            marked_v.insert(it->first.first);
            marked_v.insert(it->first.second);
        }
    }

    // create triangles from octant faces and 
    triangles.clear();
    double half_len = glm::length(overtices[ofaces[0][1]] - overtices[ofaces[0][0]]) * 0.5;
    for (int i = 0; i < (int)ofaces.size(); ++i)
    {
        int t = 0;
        while (t < 4 && marked_v.find(ofaces[i][t]) == marked_v.end())
            ++t;
        if (t == 4)
        {
            // none of the vertices of this face is non-manifold
            triangles.push_back(glm::ivec3(ofaces[i][0],ofaces[i][2],ofaces[i][1]));
            triangles.push_back(glm::ivec3(ofaces[i][0],ofaces[i][3],ofaces[i][2]));
            continue;
        }

        int ind[4];
        for (int j = 0; j < 4; ++j)
            ind[j] = ofaces[i][(t+j)%4];
        bool flag1 = marked_v.find(ind[1]) != marked_v.end();
        bool flag2 = marked_v.find(ind[2]) != marked_v.end();
        bool flag3 = marked_v.find(ind[3]) != marked_v.end();
        Grid_Index pt1 = (vertices_gridinds[ind[0]] + vertices_gridinds[ind[1]]) / 2;
        Grid_Index pt2 = (vertices_gridinds[ind[0]] + vertices_gridinds[ind[3]]) / 2;
        Grid_Index pt3 = (vertices_gridinds[ind[2]] + vertices_gridinds[ind[3]]) / 2;
        Grid_Index pt4 = (vertices_gridinds[ind[1]] + vertices_gridinds[ind[2]]) / 2;
        int ind1, ind2, ind3, ind4;
        map<Grid_Index,int>::iterator it = grid_vind.find(pt1);

        if (it == grid_vind.end())
        {
            grid_vind.insert(make_pair(pt1,overtices.size()));
            vertices_gridinds.push_back(pt1);
            ind1 = (int)overtices.size();
            overtices.push_back((overtices[ind[0]]+overtices[ind[1]])*0.5);
            vertices_orig_faces.push_back(vertices_orig_faces[ind[0]]);
        } else {
            ind1 = it->second;
        }

        it = grid_vind.find(pt2);
        if (it == grid_vind.end())
        {
            grid_vind.insert(make_pair(pt2,overtices.size()));
            vertices_gridinds.push_back(pt2);
            ind2 = (int)overtices.size();
            vertices_orig_faces.push_back(vertices_orig_faces[ind[0]]);
            overtices.push_back((overtices[ind[0]]+overtices[ind[3]])*0.5);
        } else {
            ind2 = it->second;
        }

        if (flag1 || flag2)
        {
            it = grid_vind.find(pt4);
            if (it == grid_vind.end())
            {
                grid_vind.insert(make_pair(pt4,overtices.size()));
                vertices_gridinds.push_back(pt4);
                ind4 = (int)overtices.size();
                overtices.push_back((overtices[ind[1]]+overtices[ind[2]])*0.5);
                if (flag1)
                    vertices_orig_faces.push_back(vertices_orig_faces[ind[1]]);
                else
                    vertices_orig_faces.push_back(vertices_orig_faces[ind[2]]);
            } else {
                ind4 = it->second;
            }
        }

        if (flag2 || flag3)
        {
            it = grid_vind.find(pt3);
            if (it == grid_vind.end())
            {
                grid_vind.insert(make_pair(pt3,overtices.size()));
                vertices_gridinds.push_back(pt3);
                ind3 = (int)overtices.size();
                overtices.push_back((overtices[ind[2]]+overtices[ind[3]])*0.5);
                if (flag2)
                    vertices_orig_faces.push_back(vertices_orig_faces[ind[2]]);
                else
                    vertices_orig_faces.push_back(vertices_orig_faces[ind[3]]);
            } else {
                ind3 = it->second;
            }
        }

        if (!flag1 && !flag2 && !flag3) {
            triangles.push_back(glm::ivec3(ind1,ind[2],ind[1]));
            triangles.push_back(glm::ivec3(ind2,ind[2],ind1));
            triangles.push_back(glm::ivec3(ind[3],ind[2],ind2));
        } else if (!flag1 && !flag2 && flag3) {
            triangles.push_back(glm::ivec3(ind1,ind2,ind3));
            triangles.push_back(glm::ivec3(ind1,ind3,ind[2]));
            triangles.push_back(glm::ivec3(ind1,ind[2],ind[1]));
        } else if (!flag1 && flag2 && !flag3) {
            triangles.push_back(glm::ivec3(ind1,ind4,ind[1]));
            triangles.push_back(glm::ivec3(ind1,ind2,ind4));
            triangles.push_back(glm::ivec3(ind2,ind[3],ind3));
            triangles.push_back(glm::ivec3(ind2,ind3,ind4));
        } else if (!flag1 && flag2 && flag3) {
            triangles.push_back(glm::ivec3(ind1,ind4,ind[1]));
            triangles.push_back(glm::ivec3(ind1,ind2,ind4));
            triangles.push_back(glm::ivec3(ind2,ind3,ind4));
        } else if (flag1 && !flag2 && !flag3) {
            triangles.push_back(glm::ivec3(ind1,ind2,ind4));
            triangles.push_back(glm::ivec3(ind4,ind2,ind[3]));
            triangles.push_back(glm::ivec3(ind4,ind[3],ind[2]));
        } else if (flag1 && !flag2 && flag3) {
            triangles.push_back(glm::ivec3(ind1,ind2,ind4));
            triangles.push_back(glm::ivec3(ind4,ind2,ind3));
            triangles.push_back(glm::ivec3(ind4,ind3,ind[2]));
        } else if (flag1 && flag2 && !flag3) {
            triangles.push_back(glm::ivec3(ind1,ind2,ind4));
            triangles.push_back(glm::ivec3(ind2,ind3,ind4));
            triangles.push_back(glm::ivec3(ind2,ind[3],ind3));
        } else if (flag1 && flag2 && flag3) {
            triangles.push_back(glm::ivec3(ind1,ind2,ind3));
            triangles.push_back(glm::ivec3(ind1,ind3,ind4));
        }
    }

    for (set<int>::iterator it = marked_v.begin();
        it != marked_v.end(); ++it)
    {
        glm::dvec3 p = overtices[*it];
        for (int dimx = -1; dimx < 2; dimx += 2) {
            for (int dimy = -1; dimy < 2; dimy += 2) {
                for (int dimz = -1; dimz < 2; dimz += 2) {
                    glm::dvec3 p1 = p + glm::dvec3(dimx * half_len, dimy * half_len, dimz * half_len);
                    if (tree->Is_Exterior(p1))
                    {
                        Grid_Index ind = vertices_gridinds[*it];
                        Grid_Index ind1 = ind;
                        Grid_Index ind2 = ind;
                        Grid_Index ind3 = ind;
                        ind1.id[0] += dimx;
                        ind2.id[1] += dimy;
                        ind3.id[2] += dimz;
                        if (grid_vind.find(ind1) == grid_vind.end())
                        {
                            grid_vind.insert(make_pair(ind1, overtices.size()));
                            vertices_gridinds.push_back(ind1);

                            overtices.push_back(glm::dvec3(p[0]+half_len*dimx,p[1],p[2]));
                            vertices_orig_faces.push_back(vertices_orig_faces[*it]);
                        }
                        if (grid_vind.find(ind2) == grid_vind.end())
                        {
                            grid_vind.insert(make_pair(ind2, overtices.size()));
                            vertices_gridinds.push_back(ind2);

                            overtices.push_back(glm::dvec3(p[0],p[1]+half_len*dimy,p[2]));
                            vertices_orig_faces.push_back(vertices_orig_faces[*it]);
                        }
                        if (grid_vind.find(ind3) == grid_vind.end())
                        {
                            grid_vind.insert(make_pair(ind3, overtices.size()));
                            vertices_gridinds.push_back(ind3);

                            overtices.push_back(glm::dvec3(p[0],p[1],p[2]+half_len*dimz));
                            vertices_orig_faces.push_back(vertices_orig_faces[*it]);
                        }
                        int id1 = grid_vind[ind1];
                        int id2 = grid_vind[ind2];
                        int id3 = grid_vind[ind3];
                        glm::dvec3 norm = glm::cross(overtices[id2]-overtices[id1], overtices[id3]-overtices[id1]);
                        if (glm::dot(norm, glm::dvec3(dimx,dimy,dimz)) < 0)
                            triangles.push_back(glm::ivec3(id1,id3,id2));
                        else
                            triangles.push_back(glm::ivec3(id1,id2,id3));
                    }
                }
            }
        }
    }

    // get vertices at odd grid positions (any grid coordinate is odd)
    map<int,set<pair<int,int> > > ocs, ecs;
    set<int> odds;
    set<int> evens;
    for (int i = 0; i < (int)overtices.size(); ++i)
    {
        bool flag = false;
        for (int k = 0; k < 3; ++k)
            if (vertices_gridinds[i].id[k] % 2 == 1)
                flag = true;
        if (flag) {
            odds.insert(i);
            ocs.insert(make_pair(i,set<pair<int,int> >()));
        }
    }

    // get vertices at even grid positions (all grid coordinates are even)
    for (int i = 0; i < (int)overtices.size(); ++i)
    {
        Grid_Index ind = vertices_gridinds[i];
        int flag = 0;
        while (flag < 3 && ind.id[flag] % 2 == 0)
        {
            flag++;
        }
        if (flag < 3)
            continue;
        for (int j = -2; j < 5; j += 4)
        {
            if (flag < 3)
                break;
            for (int k = 0; k < 3; ++k)
            {
                Grid_Index ind1 = ind;
                ind1.id[k] += j;
                map<Grid_Index,int>::iterator it = grid_vind.find(ind1);
                if (it == grid_vind.end())
                {
                    flag = 0;
                    break;
                }
                int y = it->second;
                unit_len = glm::length(overtices[y] - overtices[i]);
                pair<int,int> edge_id;
                if (i < y)
                    edge_id = make_pair(i,y);
                else
                    edge_id = make_pair(y,i);
                if (edge_info.find(edge_id) == edge_info.end())
                {
                    flag = 0;
                    break;
                }
            }
        }

        if (flag < 3)
            continue;

        evens.insert(i);
        ecs.insert(make_pair(i,set<pair<int,int> >()));
    }

    for (int i = 0; i < (int)triangles.size(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            int x = triangles[i][j];
            if (odds.find(x) != odds.end())
            {
                ocs[x].insert(make_pair(i,j));
            }
            if (evens.find(x) != evens.end())
            {
                ecs[x].insert(make_pair(i,j));
            }
        }
    }

    for (set<int>::iterator it = evens.begin();
        it != evens.end(); ++it)
    {
        int i = *it;
        glm::dvec3 dir;
        int count = 0;
        for (int j = 0; j < 8; ++j)
        {
            glm::dvec3 d((j&0x04)>0,(j&0x02)>0,(j&0x01)>0);
            d = d * 2.0 - glm::dvec3(1,1,1);
            d = glm::normalize(d) * (unit_len * 0.5);
            if (!tree->Is_Exterior(overtices[i] + d))
            {
                dir = glm::normalize(d);
                count += 1;
            }
        }
        if (count > 2)
            continue;
        set<pair<int,int> >& p = ecs[i];
        for (set<pair<int,int> >::iterator it1 = p.begin();
            it1 != p.end(); ++it1)
        {
            assert(triangles[it1->first][it1->second] == i);
            if (glm::dot(overtices[triangles[it1->first][(it1->second+1)%3]]-overtices[i],dir)<0)
            {
                triangles[it1->first][it1->second] = (int)overtices.size();
            }
        }
        overtices[i] += dir * (0.5 * unit_len);
        vertices_orig_faces.push_back(vertices_orig_faces[i]);
        overtices.push_back(overtices[i]);
        overtices.back() -= unit_len * dir;

    }

    for (set<int>::iterator it = odds.begin();
        it != odds.end(); ++it)
    {
        int i = *it;
        int k = 0;
        while (vertices_gridinds[i].id[k] % 2 == 0)
            k += 1;
        Grid_Index id1, id2;
        id1 = vertices_gridinds[i];
        id2 = vertices_gridinds[i];
        id1.id[k] -= 1;
        id2.id[k] += 1;
        int x = grid_vind[id1];
        int y = grid_vind[id2];
        if (x > y)
        {
            int temp = x;
            x = y;
            y = temp;
        }
        if (edge_info[make_pair(x,y)].size() > 2)
        {
            glm::dvec3 vert = overtices[x]-overtices[y];
            double len = glm::length(vert);
            vert /= len;
            glm::dvec3 dir(len*0.5,len*0.5,len*0.5);
            dir = dir - glm::dot(dir,vert)*vert;
            if (!tree->Is_Exterior(overtices[i]+dir))
            {
                dir = glm::cross(vert,dir);
            }
            dir = glm::normalize(dir);
            set<pair<int,int> >& p = ocs[i];
            for (set<pair<int,int> >::iterator it1 = p.begin();
                it1 != p.end(); ++it1)
            {
                assert(triangles[it1->first][it1->second] == i);
                if (glm::dot(overtices[triangles[it1->first][(it1->second+1)%3]]-overtices[i],dir)<0)
                {
                    triangles[it1->first][it1->second] = (int)overtices.size();
                }
            }
            overtices[i] += dir * (0.5 * len);
            vertices_orig_faces.push_back(vertices_orig_faces[i]);
            overtices.push_back(overtices[i]);
            overtices.back() -= len * dir;
        }
    }
}

void Model::is_manifold(vector<pair<int,int> >& nonmanifold_edges, vector<pair<int,int> >& inconsistent_edges, vector<int>& outofbound_faces)
{
    // edge id is ordered pair of vertex indices,
    // list of adjacent faces contains for each face:
    // a face index and edge flipped flag (flipped w.r.t. the ordered vertex pair) 
    map<pair<int,int>, list<pair<int, bool> > > edges;

    for (int i = 0; i < (int)faces.size(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            int vind1 = faces[i][j];
            int vind2 = faces[i][(j + 1) % 3];
            bool flipped = false;
            if (vind1 > vind2)
            {
                int temp = vind1;
                vind1 = vind2;
                vind2 = temp;
                flipped = true;
            }
            if (vind1 >= (int)vertices.size() || vind2 >= (int)vertices.size())
            {
                // a vertex index of the face is out of bounds
                colors[vind1] = glm::dvec3(0,1,0);
                colors[vind2] = glm::dvec3(0,1,0);
                outofbound_faces.push_back(i);
            }
            pair<int,int> edge_id = make_pair(vind1, vind2);
            auto it = edges.find(edge_id);
            if (it == edges.end())
            {
                list<pair<int, bool> > edge_faces;
                edge_faces.push_back(make_pair(i, flipped));
                edges.insert(make_pair(edge_id, edge_faces));
            } else
            {
                if (it->second.size() == 2)
                {
                    // this is the third triangle adjacent to the same edge => non-manifold
                    colors[vind1] = glm::dvec3(1,0,0);
                    colors[vind2] = glm::dvec3(1,0,0);
                    nonmanifold_edges.push_back(edge_id);
                } else
                {
                    it->second.push_back(make_pair(i, flipped));
                    bool flipped1 = it->second.front().second;
                    bool flipped2 = it->second.back().second;

                    if (flipped1 == flipped2) {
                        // the edge has the same orientation in both faces
                        // => the faces are not oriented consistently
                        colors[vind1] = glm::dvec3(0,0,1);
                        colors[vind2] = glm::dvec3(0,0,1);
                        inconsistent_edges.push_back(edge_id);
                    }
                }
            }
        }
    }

    for (auto edge_it = edges.begin(); edge_it != edges.end(); ++edge_it)
    {
        if (edge_it->second.size() == 1)
        {
            // single face adjecent to the edge => open boundary
            pair<int,int> edge_id = edge_it->first;
            int vind1 = edge_id.first;
            int vind2 = edge_id.second;
            colors[vind1] = glm::dvec3(1,0,0);
            colors[vind2] = glm::dvec3(1,0,0);
            nonmanifold_edges.push_back(edge_id);
            // std::cout << "open edge" << std::endl << std::flush;
        }
    }
}

void Model::SaveOBJ(const char* filename)
{
    std::ofstream os(filename);
    for (int i = 0; i < (int)vertices.size(); ++i)
    {
        os << "v " << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << "\n";
    }
    for (int i = 0; i < (int)faces.size(); ++i)
    {
        os << "f " << faces[i][0] + 1 << " " << faces[i][1] + 1 << " " << faces[i][2] + 1 << "\n";
    }
    os.close();
}

void Model::SavePLY(const char* filename, const std::vector<glm::dvec3>& v, const std::vector<glm::ivec3>& f, const std::vector<glm::dvec3>& c)
{
    if (v.size() != c.size()) {
        std::cout << "number of vertices and colors don't match.\n";
        exit(0);
    }

    std::ofstream os(filename);
    os << "ply\n";
    os << "format ascii 1.0\n";
    os << "element vertex " << v.size() << "\n";
    os << "property float x\n";
    os << "property float y\n";
    os << "property float z\n";
    os << "property uchar red\n";
    os << "property uchar green\n";
    os << "property uchar blue\n";
    os << "element face " << f.size() << "\n";
    os << "property list uchar int vertex_index\n";
    os << "end_header\n";
    for (int i = 0; i < (int)v.size(); ++i)
    {
        os << v[i][0] << " " << v[i][1] << " " << v[i][2]
            << " " << static_cast<int>(c[i][0]*255)
            << " " << static_cast<int>(c[i][1]*255)
            << " " << static_cast<int>(c[i][2]*255)
            << "\n";
    }
    for (int i = 0; i < (int)f.size(); ++i) {
        os << 3 << " " << f[i][0] << " " << f[i][1] << " " << f[i][2] << "\n";
    }
    os.close();
}

void Model::Save(const char* filename, bool color)
{
    std::ofstream os(filename);
    if (color)
        os << "COFF\n";
    else
        os << "OFF\n";
    os << vertices.size() << " " << faces.size() << " " << 0 << "\n";
    for (int i = 0; i < (int)vertices.size(); ++i)
    {
        if (color)
            os << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << " " << (int)(colors[i][2]*255) << " " << (int)(colors[i][1]*255) << " " << (int)(colors[i][0]*255) << " 255\n";
        else
            os << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << "\n";
    }
    double min_len = 1e30, max_len = -1e30;
    for (int i = 0; i < (int)faces.size(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            int x = faces[i][j];
            int y = faces[i][(j + 1) % 3];
            double len = glm::length(vertices[x] - vertices[y]);
            if (len < min_len)
                min_len = len;
            if (len > max_len)
                max_len = len;
        }
        os << "3 " << faces[i][0] << " " << faces[i][1] << " " << faces[i][2] << "\n";
    }
    os.close();
    std::cout << min_len << " " << max_len << "\n";
}


