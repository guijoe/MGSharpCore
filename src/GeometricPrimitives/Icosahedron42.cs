using System;
using System.Collections.Generic;


namespace MGSharp.Core.GeometricPrimitives
{
	public class Icosahedron42 : Mesh
    {
        private int numSides = 6;
	    private float frontRadius = 1f;
	    private float length = 1f;


        public Icosahedron42()
        {
            CreateIcosahedron(2);
        }

        public Icosahedron42(int n) {
            CreateIcosahedron(n);
        }
	    public Icosahedron42(int n, float r)
        {
            CreateIcosahedron(n);
            ScaleIcosahedron(r);
        }

        public Icosahedron42(int n, Vector r)
        {
            CreateIcosahedron(n);
            ScaleIcosahedron(r);
        }

        public void CreateIcosahedron(int n)
        {
            SetVertices();
            SetTriangles();
            innerRadius = new Vector(1f, 1f, 1f);

            Mesh defaultMesh = new Mesh(this);

            for (int i = 0; i < defaultMesh.faceCount(); i++)
            {
                Vector p1 = defaultMesh.faces[i].vertices[0].v;
                Vector p2 = defaultMesh.faces[i].vertices[1].v;
                Vector p3 = defaultMesh.faces[i].vertices[2].v;

                Vector[][] faces = Subdivide(n, p1, p2, p3);

                for (int j = 0; j < faces.Length; j++)
                {
                    int i0 = vertices.ContainsComponents(faces[j][0]);
                    if (i0 == -1)
                    {
                        i0 = AddVertex(faces[j][0]);
                    }

                    int i1 = vertices.ContainsComponents(faces[j][1]);
                    if (i1 == -1)
                    {
                        i1 = AddVertex(faces[j][1]);
                    }

                    int i2 = vertices.ContainsComponents(faces[j][2]);
                    if (i2 == -1)
                    {
                        i2 = AddVertex(faces[j][2]);
                    }

                    AddTriangle(i0, i1, i2);
                }
            }

            for (int i = 0; i < vertexCount(); i++)
            {
                vertices[i].v.normalize();
            }
        }

        public void SetVertices()
        {
            double goldenRatio = (1 + Math.Sqrt(5)) / 2;

            double theta = 26.56505117707799f * Math.PI / 180f;

            double stheta = Math.Sin(theta);
            double ctheta = Math.Cos(theta);

            AddVertex(new Vector(0.0f, 0.0f, -1.0f));

            double phi = Math.PI / 5f;
            for (int i = 1; i < 6; ++i)
            {
                AddVertex(new Vector(
                  ctheta * Math.Cos(phi), ctheta * Math.Sin(phi), -stheta));

                phi += 2f * Math.PI / 5f;
            }

            phi = 0f;
            for (int i = 6; i < 11; ++i)
            {
                AddVertex(new Vector(
                  ctheta * Math.Cos(phi), ctheta * Math.Sin(phi), stheta));

                phi += 2f * Math.PI / 5f;
            }

            AddVertex(new Vector(0f, 0f, 1f));
        }

        public void SetTriangles()
        {
            AddTriangle(0, 2, 1);
            AddTriangle(0, 3, 2);
            AddTriangle(0, 4, 3);
            AddTriangle(0, 5, 4);
            AddTriangle(0, 1, 5);

            AddTriangle(1, 2, 7);
            AddTriangle(2, 3, 8);
            AddTriangle(3, 4, 9);
            AddTriangle(4, 5, 10);
            AddTriangle(5, 1, 6);

            AddTriangle(1, 7, 6);
            AddTriangle(2, 8, 7);
            AddTriangle(3, 9, 8);
            AddTriangle(4, 10, 9);
            AddTriangle(5, 6, 10);

            AddTriangle(6, 7, 11);
            AddTriangle(7, 8, 11);
            AddTriangle(8, 9, 11);
            AddTriangle(9, 10, 11);
            AddTriangle(10, 6, 11);
        }

        public Vector[][] Subdivide(int n, Vector p1, Vector p2, Vector p3)
        {
            //n must be greater than or equal to two
            Vector[] side12 = new Vector[n + 1]; side12[0] = p1; side12[n] = p2;
            Vector[] side13 = new Vector[n + 1]; side13[0] = p1; side13[n] = p3;
            Vector[,] span = new Vector[n + 1, n];

            //subdivision points
            for (int i = 1; i < n; i++)
            {
                side12[i] = new Vector(p1.x + ((float)i / n) * (p2.x - p1.x),
                                      p1.y + ((float)i / n) * (p2.y - p1.y), p1.z + ((float)i / n) * (p2.z - p1.z));

                side13[i] = new Vector(p1.x + ((float)i / n) * (p3.x - p1.x),
                                      p1.y + ((float)i / n) * (p3.y - p1.y), p1.z + ((float)i / n) * (p3.z - p1.z));
            }

            for (int i = 0; i < n + 1; i++)
            {
                for (int j = 1; j < n - i; j++)
                {
                    span[i, j] = new Vector(side12[n - i].x + ((float)j / (n - i)) * (side13[n - i].x - side12[n - i].x),
                                           side12[n - i].y + ((float)j / (n - i)) * (side13[n - i].y - side12[n - i].y),
                                           side12[n - i].z + ((float)j / (n - i)) * (side13[n - i].z - side12[n - i].z));
                }
            }

            Vector[][] faces = new Vector[n * n][];
            for (int i = 0; i < n * n; i++)
            {
                faces[i] = new Vector[3];
            }
            Vector[] ps = new Vector[3];
            int nfaces = 0;

            //top four subfaces
            faces[nfaces][0] = side12[0]; faces[nfaces][1] = side13[1]; faces[nfaces][2] = side12[1];
            nfaces++;
            //Debug.Log(faces[0][0] + ", " + faces[0][1] + ", " + faces[0][2]);
            faces[nfaces][0] = side12[1]; faces[nfaces][2] = side12[2]; faces[nfaces][1] = span[n - 2, 1];
            nfaces++;
            //Debug.Log(faces[0][0] + ", " + faces[0][1] + ", " + faces[0][2]);
            faces[nfaces][0] = side12[1]; faces[nfaces][2] = span[n - 2, 1]; faces[nfaces][1] = side13[1];
            nfaces++;
            //Debug.Log(faces[0][0] + ", " + faces[0][1] + ", " + faces[0][2]);
            faces[nfaces][0] = side13[1]; faces[nfaces][1] = side13[2]; faces[nfaces][2] = span[n - 2, 1];
            nfaces++;
            //Debug.Log(faces[0][0] + ", " + faces[0][1] + ", " + faces[0][2]);


            //the rest of the subfaces
            for (int i = 2; i < n; i++)
            {
                //side 12
                faces[nfaces][0] = side12[i]; faces[nfaces][2] = side12[i + 1]; faces[nfaces][1] = span[n - i - 1, 1];
                nfaces++;
                faces[nfaces][0] = side12[i]; faces[nfaces][1] = span[n - i, 1]; faces[nfaces][2] = span[n - i - 1, 1];
                nfaces++;

                //center
                for (int j = 1; j < n - i + 1; j++)
                {
                    faces[nfaces][0] = span[i - 2, j + 1]; faces[nfaces][1] = span[i - 2, j]; faces[nfaces][2] = span[i - 1, j];
                    nfaces++;
                    if (i > 2)
                    {
                        faces[nfaces][0] = span[i - 2, j + 1]; faces[nfaces][2] = span[i - 2, j]; faces[nfaces][1] = span[i - 3, j + 1];
                        nfaces++;
                    }
                }

                //side 13
                faces[nfaces][0] = side13[i]; faces[nfaces][1] = side13[i + 1]; faces[nfaces][2] = span[n - i - 1, i];
                nfaces++;
                faces[nfaces][0] = side13[i]; faces[nfaces][2] = span[n - i, i - 1]; faces[nfaces][1] = span[n - i - 1, i];
                nfaces++;
            }

            for (int i = 0; i < faces.Length; i++)
            {
                for (int j = 0; j < faces[i].Length; j++)
                {
                    faces[i][j] = new Vector((float)faces[i][j].x, (float)faces[i][j].y, (float)faces[i][j].z);
                }
            }

            return faces;
        }

        protected void ScaleIcosahedron(float r)
        {
            innerRadius0 = new Vector(r, 1f, r);
            for (int i=0; i<vertexCount(); i++)
            {
                vertices[i].v.x *= r;
                vertices[i].v.z *= r;
            }
            innerRadius.x *= r;
            innerRadius.z *= r;
        }

        protected void ScaleIcosahedron(Vector r)
        {
            innerRadius0 = r;
            for (int i = 0; i < vertexCount(); i++)
            {
                vertices[i].v.x *= r.x;
                vertices[i].v.y *= r.y;
                vertices[i].v.z *= r.z;
            }
            innerRadius.x *= r.x;
            innerRadius.y *= r.y;
            innerRadius.z *= r.z;
        }        
    }
}