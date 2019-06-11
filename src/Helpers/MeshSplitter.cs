using System;
using System.Collections.Generic;

namespace MGSharp.Core.GeometricPrimitives
{
    class MeshSplitter
    {
        // Mesh info
        public Mesh mesh;
        public Mesh meshUp;
        public Mesh meshDown;
        
        // mesh copy
        public int vertexCount;
        public Vector[] vertices;
        public Vector[] normals;
        
        public VertexSet wsVerts;
        public FaceSet triangles;

        public VertexSet vertsUp;
        public VertexSet vertsDown;
        
        protected Plane splitPlane;

        public MeshSplitter(Mesh m, Plane plane)
        {
            mesh = m;
            splitPlane = plane;
            meshUp = new Mesh();
            meshDown = new Mesh();
        }

        public void MeshInitialize()
        {
            // store mesh data
            vertexCount = mesh.vertexCount();
            
            vertices = new Vector[vertexCount];
            normals = new Vector[vertexCount];
            Vector mean = Vector.zero;
            
            for (int i = 0; i < vertexCount; i++)
            {
                vertices[i] = mesh.vertices[i].v;
                normals[i] = mesh.vertices[i].normal;
                mean += vertices[i];
            }
            
            wsVerts = mesh.vertices;
            triangles = mesh.faces;
            //normals = mesh.no;
            
            
            if (vertices.Length != 0) vertsUp = new VertexSet();
            if (vertices.Length != 0) vertsDown = new VertexSet();
        }

        public void MeshSplit()
        {
            Vector[] normalsUp = new Vector[vertices.Length];
            Vector[] normalsDown = new Vector[vertices.Length];

            VertexSet verticesUp = new VertexSet();
            VertexSet verticesDown = new VertexSet();

            for (int i = 0; i < mesh.vertexCount(); i++)
            {
                verticesUp.add(vertices[i]);
                verticesDown.add(vertices[i]);
                normalsUp[i] = normalsDown[i] = normals[i];

                meshUp.AddVertex(mesh.getVertex(i).Clone().v);
                meshDown.AddVertex(mesh.getVertex(i).Clone().v);
            }

            for (int i = 0; i < mesh.vertexCount(); i++)
            {
                if (splitPlane.PointSide(mesh.getVertex(i).v) > 0)
                {
                    //meshDown.SetVertex(i, splitPlane.PointOrthogonalProjection(meshDown.getVertex(i).v) - splitPlane.Normal * MGModels.MGModel.DCol/2);
                    meshDown.SetVertex(i, splitPlane.PointOrthogonalProjection(meshDown.getVertex(i).v));
                }
                else if (splitPlane.PointSide(mesh.getVertex(i).v) < 0)
                {
                    //meshUp.SetVertex(i, splitPlane.PointOrthogonalProjection(meshUp.getVertex(i).v) + splitPlane.Normal * MGModels.MGModel.DCol/2);
                    meshUp.SetVertex(i, splitPlane.PointOrthogonalProjection(meshUp.getVertex(i).v));
                }
            }
        }

        public bool IsMeshSplit()
        {
            return HasMeshUpper() && HasMeshLower();
        }

        public bool HasMeshUpper()
        {
            return vertsUp.getCount() > 0;
        }

        public bool HasMeshLower()
        {
            return vertsDown.getCount() > 0;
        }

        public Mesh CreateMeshUpper()
        {
            MyCreateMesh(true);
            return meshUp;
        }

        public Mesh CreateMeshLower()
        {
            MyCreateMesh(false);
            return meshDown;
        }

        private void MyCreateMesh(bool upperMesh)
        {
            if (upperMesh)
            {
                for (int i = 0; i < mesh.faceCount(); i++)
                {
                    meshUp.AddTriangle(mesh.faces[i].vertices[0].pos,
                                    mesh.faces[i].vertices[1].pos,
                                    mesh.faces[i].vertices[2].pos);
                }
                //Console.WriteLine("Mesh up Triangles: " + meshUp.faceCount());
            }
            else
            {
                for (int i = 0; i < mesh.faceCount(); i++)
                {
                    meshDown.AddTriangle(mesh.faces[i].vertices[0].pos,
                                mesh.faces[i].vertices[1].pos,
                                mesh.faces[i].vertices[2].pos);
                }
                //Console.WriteLine("Mesh down Triangles: " + meshDown.faceCount());
            }
        }
    }
}