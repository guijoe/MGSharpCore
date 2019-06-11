using System;

namespace MGSharp.Core.GeometricPrimitives
{
    public class Face
    {
        public Face()
        {
            label = 0;
            edges = new EdgeSet();
            vertices = new VertexSet();
            normal = new Vector();
            vertices.reserve(3);
            edges.reserve(3);
        }

        public Face Clone()
        {
            Face f = new Face();
            f.edges = edges.Clone();
            f.label = label;
            f.normal = (Vector)normal.Clone();
            f.vertices = vertices.Clone();
            return f;
        }

        public Edge adjacent(Face f)
        {
            int i, j;
            for (i = 0; i < edges.getCount(); i++)
            {
                for (j = 0; j < f.edges.getCount(); j++)
                {
                    if (edges[i].Equals(f.edges[j]))
                        return edges[i];
                }
            }
            return null;
        }

        public Vertex OffEdge(Edge e)
        {
            int i;
            for (i = 0; i < vertices.getCount(); i++)
            {
                if ((vertices[i] != e.ends[0]) &&
                    (vertices[i] != e.ends[1]))
                {
                    return vertices[i];
                }
            }
            return null;
        }

        public void replace(Edge eold, Edge enew)
        {
            int i;
            for (i = 0; i < edges.getCount(); i++)
            {
                if (edges[i].Equals(eold))
                {
                    eold.remove(this);
                    edges[i] = enew;
                    enew.add(this);
                }
            }
        }

        public void replace(Vertex vold, Vertex vnew)
        {
            int i;
            for (i = 0; i < vertices.getCount(); i++)
            {
                if (vertices[i].Equals(vold))
                {
                    vold.faces.remove(this);
                    vertices[i] = vnew;
                    vnew.faces.add(this);
                }
            }
        }
        public void disconnect()
        {
            int i;
            for (i = 0; i < vertices.getCount(); i++)
            {
                vertices[i].faces.remove(this);
            }
            for (i = 0; i < edges.getCount(); i++)
            {
                edges[i].remove(this);
            }
            edges.clear();
            vertices.clear();
        }
        public void ComputeNormal()
        {
            normal = (vertices[1].v - vertices[0].v) ^
                (vertices[2].v - vertices[0].v);
            normal.normalize();
        }

        public void ComputeNormal2()
        {
            normal = (vertices[1].v - vertices[0].v) ^
                (vertices[2].v - vertices[0].v);
            normal.normalize();
        }
        public void ComputeNormalTex()
        {
            normal = (vertices[1].tex - vertices[0].tex) ^
                (vertices[2].tex - vertices[0].tex);
            normal.normalize();
        }
        public double area()
        {
            Vector e1, e2;
            e1 = vertices[1].v - vertices[0].v;
            e2 = vertices[2].v - vertices[0].v;
            return Math.Abs((e1 ^ e2).norm() / 2);
        }
        public double texarea()
        {
            Vector e1, e2;
            e1 = vertices[1].tex - vertices[0].tex;
            e2 = vertices[2].tex - vertices[0].tex;
            return Math.Abs((e1 ^ e2).norm() / 2);
        }
        //See: Texture mapping progressive meshes,
        // by Sander et al.
        public double stretch()
        {
            Vector Ss, St;
            double a, b, c;
            double A;
            Vertex q1, q2, q3;
            q1 = vertices[0];
            q2 = vertices[1];
            q3 = vertices[2];

            // x= s coordinate
            // y= t coordinate
            //Area of parameterised triangle
            A = ((q2.tex.x - q1.tex.x) * (q3.tex.y - q1.tex.y) -
               (q3.tex.x - q1.tex.x) * (q2.tex.y - q1.tex.y)) / 2.0;
            A = Math.Abs(A);
            if (A < 1e-10) A = 1e-10;
            Ss = q1.v * (q2.tex.y - q3.tex.y) +
                q2.v * (q3.tex.y - q1.tex.y) +
                q3.v * (q1.tex.y - q2.tex.y);
            St = q1.v * (q3.tex.x - q2.tex.x) +
                q2.v * (q1.tex.x - q3.tex.x) +
                q3.v * (q2.tex.x - q1.tex.x);
            Ss /= 2 * A;
            St /= 2 * A;
            a = Ss * Ss;
            b = Ss * St;
            c = St * St;
            //L^2 stretch metric from Sander et al.
            return Math.Sqrt((a + c) / 2.0);
        }
        public bool boundary()
        {
            if (edges[0].boundary()) return true;
            if (edges[1].boundary()) return true;
            if (edges[2].boundary()) return true;
            return false;
        }
        //neglects to check if vertices are in the
        // same order. But who wants a mesh where
        //	the same vertices occur in several faces
        //	but cross/overlap in some way?
        public override bool Equals(object obj)
        {
            int i, j, sum;
            int[] tag;
            Face f2;
            f2 = (Face)obj;
            tag = new int[vertices.getCount()];
            for (i = 0; i < vertices.getCount(); i++)
                tag[i] = 0;
            if (vertices.getCount() != f2.vertices.getCount())
                return false;
            for (i = 0; i < vertices.getCount(); i++)
            {
                for (j = 0; j < f2.vertices.getCount(); j++)
                {
                    if (vertices[i].Equals(f2.vertices[j]))
                        tag[i] = 1;
                }
            }
            sum = 0;
            for (i = 0; i < vertices.getCount(); i++)
                sum += tag[i];
            if (sum == vertices.getCount())
                return true;
            else return false;
        }
        //for fast access
        public VertexSet vertices;
        public EdgeSet edges;
        public Vector normal;
        public int label;

        public bool ArePointsOnSameSide(Vector p1, Vector p2)
        {
            Vector normal = (vertices[1].v - vertices[0].v)^(vertices[2].v - vertices[0].v);
            double p1Side = (normal * (p1 - vertices[0].v));
            double p2Side = (normal * (p2 - vertices[0].v));

            return (p1Side * p2Side >= -MGModels.MGModel.epsilon);
        }

        public Vertex ThirdVertex(Vertex v1, Vertex v2)
        {
            if((v1.Equals(vertices[0]) && v2.Equals(vertices[1])) ||
                (v1.Equals(vertices[1]) && v2.Equals(vertices[0])))
            {
                return vertices[2];
            }
            if ((v1.Equals(vertices[1]) && v2.Equals(vertices[2])) ||
                (v1.Equals(vertices[2]) && v2.Equals(vertices[1])))
            {
                return vertices[0];
            }
            if ((v1.Equals(vertices[2]) && v2.Equals(vertices[0])) ||
                (v1.Equals(vertices[0]) && v2.Equals(vertices[2])))
            {
                return vertices[1];
            }
            else
            {
                return new Vertex();
            }
        }

        public override string ToString()
        {
            string face = "";
            for(int i=0; i<vertices.getCount(); i++)
            {
                face += vertices[i].pos + ",";
            }
            face = face.Remove(face.Length - 1);
            return face;
        }

        public Vector centre { get
            {
                return (vertices[0].v + vertices[1].v + vertices[2].v) / 3;
            }
        }

        
    }
}
