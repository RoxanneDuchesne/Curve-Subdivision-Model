#include "obj_mesh.hpp"

#include <algorithm>
#include <cassert>

using namespace math;

namespace geometry {

std::vector<math::Vec3f>
calculateTriangleNormals(std::vector<IndicesTriangle> const &indexTriangles,
                         std::vector<math::Vec3f> const &vertices) {
   std::vector<Vec3f> tmp;

   // reverse vector using operator[]:
   for (auto i : indexTriangles)
   {
     unsigned int vertex_id_a;
     unsigned int vertex_id_b;
     unsigned int vertex_id_c;
     vertex_id_a = i[0][0];
     vertex_id_b = i[1][0];
     vertex_id_c = i[2][0];

     Vec3f A = vertices[vertex_id_a];
     Vec3f B = vertices[vertex_id_b];
     Vec3f C = vertices[vertex_id_c];

     Vec3f AB = B - A;
     Vec3f AC = C - A;


     tmp.push_back(normalized(AB^AC));
   }

   return tmp;
}

std::vector<math::Vec3f>
calculateVertexNormals(std::vector<IndicesTriangle> const &indexTriangles,
                       std::vector<math::Vec3f> const &vertices) {

    std::vector<Vec3f> tmp (vertices.size());
    std::vector<math::Vec3f> TriangleNormals = calculateTriangleNormals(indexTriangles, vertices);

    unsigned int cur_index = 0;

    // reverse vector using operator[]:
    for (auto i : indexTriangles)
    {
      unsigned int vertex_id_a;
      unsigned int vertex_id_b;
      unsigned int vertex_id_c;
      vertex_id_a = i[0][0];
      vertex_id_b = i[1][0];
      vertex_id_c = i[2][0];

      tmp[vertex_id_a] = normalized(tmp[vertex_id_a] + TriangleNormals[cur_index]);
      tmp[vertex_id_b] = normalized(tmp[vertex_id_b] + TriangleNormals[cur_index]);
      tmp[vertex_id_c] = normalized(tmp[vertex_id_c] + TriangleNormals[cur_index]);

      cur_index += 1;
    }

    return tmp;
}

std::vector<math::Vec3f>
calculateVertexNormals(std::vector<IndicesTriangle> const &indexTriangles,
                       std::vector<math::Vec3f> const &vertices,
                       std::vector<math::Vec3f> const &triangleNormals) {

    std::vector<Vec3f> tmp (vertices.size());

    unsigned int cur_index = 0;

    // reverse vector using operator[]:
    for (auto i : indexTriangles)
    {
      unsigned int vertex_id_a;
      unsigned int vertex_id_b;
      unsigned int vertex_id_c;
      vertex_id_a = i[0][0];
      vertex_id_b = i[1][0];
      vertex_id_c = i[2][0];

      tmp[vertex_id_a] = normalized(tmp[vertex_id_a] + triangleNormals[cur_index]);
      tmp[vertex_id_b] = normalized(tmp[vertex_id_b] + triangleNormals[cur_index]);
      tmp[vertex_id_c] = normalized(tmp[vertex_id_c] + triangleNormals[cur_index]);

      cur_index += 1;
    }

    return tmp;
}

} // namespace geometry
