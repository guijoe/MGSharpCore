using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MGSharp.Core.GeometricPrimitives
{
    class Ellipsoid : Mesh
    {
        public int gridSize = 2;
        public int ySize = 3;

        public float radius = 1;
        public float roundness = .5f;
        int quads = 0;

        public Ellipsoid()
        {
            CreateVertices();
            CreateTriangles();
            //RoundCube();
            ScaleCube(.5f);
        }

        private void CreateVertices()
        {
            int cornerVertices = 8;
            int edgeVertices = (gridSize + ySize + gridSize - 3) * 4;
            int faceVertices = (
                (gridSize - 1) * (gridSize - 1) +
                (ySize - 1) * (ySize - 1) +
                (gridSize - 1) * (gridSize - 1)) * 2;

            //vertices = new List<Vector>(cornerVertices + edgeVertices + faceVertices); //8 + 12(s-1)+6(s-1)^2
            //normals = new List<Vector>(cornerVertices + edgeVertices + faceVertices);
            //cubeUV = new Color32[cornerVertices + edgeVertices + faceVertices];

            int v = 0;
            for (int y = 0; y <= ySize; y++)
            {
                for (int x = 0; x <= gridSize; x++)
                {
                    SetVertex(v++, x, y, 0);
                }
                for (int z = 1; z <= gridSize; z++)
                {
                    SetVertex(v++, gridSize, y, z);
                }
                for (int x = gridSize - 1; x >= 0; x--)
                {
                    SetVertex(v++, x, y, gridSize);
                }
                for (int z = gridSize - 1; z > 0; z--)
                {
                    SetVertex(v++, 0, y, z);
                }
            }

            for (int z = 1; z < gridSize; z++)
            {
                for (int x = 1; x < gridSize; x++)
                {
                    SetVertex(v++, x, ySize, z);
                }
            }

            for (int z = 1; z < gridSize; z++)
            {
                for (int x = 1; x < gridSize; x++)
                {
                    SetVertex(v++, x, 0, z);
                }
            }

            //mesh.vertices = vertices.ToArray();
            //mesh.normals = normals.ToArray();
        }

        private void SetVertex(int i, int x, int y, int z)
        {
            Vector v = new Vector(x, y, z); // * 2f / gridSize - Vector.one;
            AddVertex(v);
        }

        private void CreateTriangles()
        {
            int quads = (gridSize * gridSize + ySize * ySize + gridSize * gridSize) * 2;
            int[] triangles = new int[quads * 6];

            int ring = (gridSize + gridSize) * 2;
            int t = 0, v = 0;

            for (int y = 0; y < ySize; y++, v++)
            {
                for (int q = 0; q < ring - 1; q++, v++)
                {
                    t = SetQuad(triangles, t, v, v + 1, v + ring, v + ring + 1);
                }
                t = SetQuad(triangles, t, v, v - ring + 1, v + ring, v + 1);
            }

            t = CreateTopFace(triangles, t, ring);
            t = CreateBottomFace(triangles, t, ring);

            /*
	        if(gridSize == 2)
	        {
	            t = SetQuad(triangles, t, 17, 18, 24, 19);
	        }
	        else
	        {
	            t = SetQuad(triangles, t, 38, 39, 49, 40);
	        }
            //*/
            //*
            if (gridSize == 2 && ySize == 3)
            {
                t = SetQuad(triangles, t, 25, 26, 32, 27);
            }
            //*/
            //mesh.triangles = triangles;
        }

        private int CreateTopFace(int[] triangles, int t, int ring)
        {
            int v = ring * ySize;
            for (int x = 0; x < gridSize - 1; x++, v++)
            {
                t = SetQuad(triangles, t, v, v + 1, v + ring - 1, v + ring);
            }
            int vMin = ring * (ySize + 1) - 1;
            int vMid = vMin + 1;
            int vMax = v + 2;

            for (int z = 1; z < gridSize - 1; z++, vMin--, vMid++, vMax++)
            {
                t = SetQuad(triangles, t, vMin, vMid, vMin - 1, vMid + gridSize - 1);
                for (int x = 1; x < gridSize - 1; x++, vMid++)
                {
                    t = SetQuad(
                        triangles, t,
                        vMid, vMid + 1, vMid + gridSize - 1, vMid + gridSize);
                }
                t = SetQuad(triangles, t, vMid, vMax, vMid + gridSize - 1, vMax + 1);
            }

            int vTop = vMin - 2;
            t = SetQuad(triangles, t, vMin, vMid, vTop + 1, vTop);
            for (int x = 1; x < gridSize - 1; x++, vTop--, vMid++)
            {
                t = SetQuad(triangles, t, vMid, vMid + 1, vTop, vTop - 1);
            }
            t = SetQuad(triangles, t, vMid, vTop - 2, vTop, vTop - 1);
            return t;
        }

        private int CreateBottomFace(int[] triangles, int t, int ring)
        {
            int v = 1;
            int vMid = vertices.getCount() - (gridSize - 1) * (gridSize - 1);
            t = SetQuad(triangles, t, ring - 1, vMid, 0, 1);
            for (int x = 1; x < gridSize - 1; x++, v++, vMid++)
            {
                t = SetQuad(triangles, t, vMid, vMid + 1, v, v + 1);
            }
            t = SetQuad(triangles, t, vMid, v + 2, v, v + 1);

            int vMin = ring - 2;
            vMid -= gridSize - 2;
            int vMax = v + 2;

            for (int z = 1; z < gridSize - 1; z++, vMin--, vMid++, vMax++)
            {
                t = SetQuad(triangles, t, vMin, vMid + gridSize - 1, vMin + 1, vMid);
                for (int x = 1; x < gridSize - 1; x++, vMid++)
                {
                    t = SetQuad(
                        triangles, t,
                        vMid + gridSize - 1, vMid + gridSize, vMid, vMid + 1);
                }
                t = SetQuad(triangles, t, vMid + gridSize - 1, vMax + 1, vMid, vMax);
            }

            int vTop = vMin - 1;
            t = SetQuad(triangles, t, vTop + 1, vTop, vTop + 2, vMid);
            for (int x = 1; x < gridSize - 1; x++, vTop--, vMid++)
            {
                t = SetQuad(triangles, t, vTop, vTop - 1, vMid, vMid + 1);
            }
            t = SetQuad(triangles, t, vTop, vTop - 1, vMid, vTop - 2);

            return t;
        }

        private int SetQuad(int[] triangles, int i, int v00, int v10, int v01, int v11)
        {
            //if (v00 == 39) {
            //Debug.Log("("+v00+", "+v10+", "+v01+", "+v11+")");
            //}

            AddTriangle(v00, v01, v10);
            AddTriangle(v10, v01, v11);

            //triangles[i] = v00;
            //triangles[i + 1] = triangles[i + 4] = v01;
            //triangles[i + 2] = triangles[i + 3] = v10;
            //triangles[i + 5] = v11;

            return i + 6;
        }

        //*
        private void ScaleCube(float r)
        {
            innerRadius = new Vector(r, 1f, r);
            Console.WriteLine(innerRadius);
            for (int i = 0; i < vertices.getCount(); i++)
            {
                Vector vec = vertices[i].v;
                vec.x = 2 * innerRadius.x * vec.x / gridSize - innerRadius.x;
                vec.y = 2 * innerRadius.y * vec.y / ySize - innerRadius.y;
                vec.z = 2 * innerRadius.z * vec.z / gridSize - innerRadius.z;
                vertices[i].v = vec;
            }
        }

        private void RoundCube()
        {
            Vector CI = Vector.zero, C = Vector.zero, M = Vector.zero, I = Vector.zero, V = Vector.zero, IV = Vector.zero;

            for (int i = 0; i < vertices.getCount(); i++)
            {
                V = vertices[i].v;
                M = vertices[i].v;

                float x = (float)vertices[i].v.x, y = (float)vertices[i].v.y, z = (float)vertices[i].v.z;

                if (x < roundness)
                {
                    CI.x = roundness;
                    C.x = 0;
                }
                else if (x > gridSize - roundness)
                {
                    CI.x = -roundness;
                    C.x = gridSize;
                }
                if (y < roundness)
                {
                    CI.y = roundness;
                    C.y = 0;
                }
                else if (y > ySize - roundness)
                {
                    CI.y = -roundness;
                    C.y = ySize;
                }
                if (z < roundness)
                {
                    CI.z = roundness;
                    C.z = 0;
                }
                else if (z > gridSize - roundness)
                {
                    CI.z = -roundness;
                    C.z = gridSize;
                }

                if (CI.x * CI.y * CI.z != 0)
                {
                    IV = V - C - CI; //CV - CI
                    IV.normalize();
                    M = C + CI + roundness * IV;
                    vertices[i].v = M;
                }
                CI = Vector.zero;
            }
            //mesh.vertices = vertices.ToArray();
        }
        //*/
    }
}