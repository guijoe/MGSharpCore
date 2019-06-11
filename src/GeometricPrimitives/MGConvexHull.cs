using System;
using System.Linq;
using System.Collections.Generic;
using MGSharp.MIConvexHull;
//using MGSharp.MIConvexHull;

namespace MGSharp.Core.GeometricPrimitives
{
    class MGConvexHull
    {
        public static Mesh GenerateConvexHull(Mesh[] meshes) {

            Mesh CV = new Mesh();

            List<MIVertex> vertices = new List<MIVertex>();
            for(int i=0; i<meshes.Length; i++)
            {
                //Console.WriteLine(meshes[i].vertexCount());
                for (int j=0; j < meshes[i].vertexCount(); j++)
                {
                    vertices.Add(new MIVertex(
                        meshes[i].vertices[j].v.x,
                        meshes[i].vertices[j].v.y,
                        meshes[i].vertices[j].v.z
                    ));

                    
                }
            }

            ConvexHull<MIVertex, MIFace> convexHull = ConvexHull.Create<MIVertex, MIFace>(vertices);

            List<MIVertex> convexHullVertices = convexHull.Points.ToList();
            List<MIFace> faces = convexHull.Faces.ToList();

            for(int i=0; i<faces.Count; i++)
            {
                CV.AddVertex(new Vector(faces[i].Vertices[0].Position[0],
                                        faces[i].Vertices[0].Position[1],
                                        faces[i].Vertices[0].Position[2]
                ));
                CV.AddVertex(new Vector(faces[i].Vertices[1].Position[0],
                                        faces[i].Vertices[1].Position[1],
                                        faces[i].Vertices[1].Position[2]
                ));
                CV.AddVertex(new Vector(faces[i].Vertices[2].Position[0],
                                        faces[i].Vertices[2].Position[1],
                                        faces[i].Vertices[2].Position[2]
                ));
                CV.AddTriangle(3*i,3*i+1,3*i+2);
            }

            return CV;
        }

        public static Mesh GenerateConvexHull(Vector[] points)
        {

            Mesh CV = new Mesh();

            List<MIVertex> vertices = new List<MIVertex>();
            for (int j = 0; j < points.Length; j++)
            {
                vertices.Add(new MIVertex(
                    points[j].x,
                    points[j].y,
                    points[j].z
                ));
            }

            ConvexHull<MIVertex, MIFace> convexHull = ConvexHull.Create<MIVertex, MIFace>(vertices);
            Console.WriteLine(convexHull.Points.ToList().Count);

            List<MIVertex> convexHullVertices = convexHull.Points.ToList();
            List<MIFace> faces = convexHull.Faces.ToList();

            for (int i = 0; i < faces.Count; i++)
            {
                //Console.WriteLine(faces[i].Vertices.ToList().Count);
                CV.AddVertex(new Vector(faces[i].Vertices[0].Position[0],
                                        faces[i].Vertices[0].Position[1],
                                        faces[i].Vertices[0].Position[2]
                ));
                CV.AddVertex(new Vector(faces[i].Vertices[1].Position[0],
                                        faces[i].Vertices[1].Position[1],
                                        faces[i].Vertices[1].Position[2]
                ));
                CV.AddVertex(new Vector(faces[i].Vertices[2].Position[0],
                                        faces[i].Vertices[2].Position[1],
                                        faces[i].Vertices[2].Position[2]
                ));
                CV.AddTriangle(3 * i, 3 * i + 1, 3 * i + 2);
            }

            return CV;
        }

        public static Mesh GenerateConvexHull(List<Vertex> vertices)
        {

            Mesh CV = new Mesh();

            //var convexHull = ConvexHull.Create<MIVertex, Face>(vertices);
            //convexHullVertices = convexHull.Points.ToList();
            //faces = convexHull.Faces.ToList();


            return CV;
        }
    }
}
