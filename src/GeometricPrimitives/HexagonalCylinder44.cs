using System;
using System.Collections.Generic;


namespace MGSharp.Core.GeometricPrimitives
{
	public class HexagonalCylinder44 : Mesh
    {
        private int numSides = 6;
	    private float frontRadius = 1f;
	    private float length = 1f;

        public HexagonalCylinder44() { }
	    public HexagonalCylinder44(float r)
        {
            SetVertices();
            SetTriangles();
            innerRadius = new Vector(Math.Sqrt(3f) / 2f, 1f, 0.75f);
            ScaleHexagon(r);
        }

        public HexagonalCylinder44(Vector r)
        {
            SetVertices();
            SetTriangles();
            innerRadius = new Vector(Math.Sqrt(3f) / 2f, 1f, 0.75f);
            ScaleHexagon(r);
        }

        protected void SetVertices()
	    {
	        //building block vertices
	        Vector[] vs = new Vector[7 * numSides + 2];

	        //Set the vs
            for(int j=0; j<7; j+=2)
            {
                float y = -j * length / 3 + length / 2f;
                for (int i = 0; i < numSides; i++)
                {
                    float angle = (float)(2 * Math.PI * (i + 0.5f) / numSides);
                    vs[j * numSides + i] = new Vector(frontRadius * Math.Cos(angle), y, frontRadius * Math.Sin(angle));
                    //vs[j * numSides + i + 1] = new Vector(frontRadius * Math.Cos(angle), y, frontRadius * Math.Sin(angle));
                }
            }

            for (int j = 1; j < 7; j += 2)
            {
                for (int i = 0; i < numSides; i++)
                {
                    //*
                    vs[j * numSides + i] = 0.25f * 
                                        (vs[(j-1) * numSides + i] + 
                                        vs[(j - 1) * numSides + (i == (numSides - 1) ? 0 : i + 1)] + 
                                        vs[(j + 1) * numSides + i] + 
                                        vs[(j + 1) * numSides + (i == (numSides - 1) ? 0 : i + 1)]);
                    //*/

                    /*
                    vs[j * numSides + i + 1] = 0.25f *
                                       (vs[(j - 1) * numSides + i + 1] +
                                        vs[(j - 1) * numSides + (i == (numSides - 1) ? 0 : i + 1) + 1] +
                                        vs[(j + 1) * numSides + i + 1] +
                                        vs[(j + 1) * numSides + (i == (numSides - 1) ? 0 : i + 1) + 1]);
                    */
                }
            }

            vs[42] = (vs[1] + vs[2] + vs[3] + vs[4] + vs[5] + vs[6]) / 6;
            vs[43] = (vs[41] + vs[40] + vs[39] + vs[38] + vs[37] + vs[36]) / 6;

            //for (int i = 0; i < 42; i++)
            for (int i = 0; i < 44; i++)
            {
                AddVertex(vs[i]);
            }
  	    }

        protected void SetTriangles()
        {
            //first end
            for (int i = 0; i < numSides - 1; i++)
            //for (int i = 0; i < numSides; i++)
            {
                AddTriangle(42, i + 1, i);
                //AddTriangle(0, i + 1, i);
            }
           AddTriangle(42, 0, 5);

            //middle
            for (int j = 1; j < 7; j += 2)
            {
                for (int i = 0; i < numSides; i++)
                {
                    //*
                    AddTriangle(j * numSides + i, 
                                (j - 1) * numSides + i, 
                                (j - 1) * numSides + (i == (numSides - 1) ? 0 : i + 1));

                    AddTriangle(j * numSides + i, 
                                (j - 1) * numSides + (i == (numSides - 1) ? 0 : i + 1), 
                                (j + 1) * numSides + (i == (numSides - 1) ? 0 : i + 1));

                    AddTriangle(j * numSides + i,
                                (j + 1) * numSides + (i == (numSides - 1) ? 0 : i + 1),
                                (j + 1) * numSides + i);

                    AddTriangle(j * numSides + i,
                                (j + 1) * numSides + i,
                                (j - 1) * numSides + i);
                    //*/

                    /*
                    AddTriangle(j * numSides + i + 1,
                                (j - 1) * numSides + i + 1,
                                (j - 1) * numSides + (i == (numSides - 1) ? 0 : i + 1) + 1);

                    AddTriangle(j * numSides + i + 1,
                                (j - 1) * numSides + (i == (numSides - 1) ? 0 : i + 1 + 1),
                                (j + 1) * numSides + (i == (numSides - 1) ? 0 : i + 1) + 1);

                    AddTriangle(j * numSides + i + 1,
                                (j + 1) * numSides + (i == (numSides - 1) ? 0 : i + 1) + 1,
                                (j + 1) * numSides + i + 1);

                    AddTriangle(j * numSides + i + 1,
                                (j + 1) * numSides + i + 1,
                                (j - 1) * numSides + i + 1);
                    */
                }
            }

            //other end - opposite way round so face points outwards
            for (int i = 0; i < numSides - 1; i++)
            //for (int i = 1; i < numSides; i++)
            {
                //There are numSides triangles in the first end, 4*numSides triangles in the middle, so this starts on 5*numSides
                
                /*
                AddTriangle(6 * numSides,
                            6 * numSides + i,
                            6 * numSides + i + 1);
                */
                
                //*
                AddTriangle(43,
                            6 * numSides + i,
                            6 * numSides + i + 1);
                //*/
            }
            AddTriangle(43, 41, 36);
        }

        protected void ScaleHexagon(float r)
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

        protected void ScaleHexagon(Vector r)
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