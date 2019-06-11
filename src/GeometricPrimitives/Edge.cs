using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MGSharp.Core.GeometricPrimitives
{
    public class Edge
    {
        public Edge()
        {
            ends = new Vertex[2];
            ends[0] = ends[1] = null;
            label = 0;
            faces = new FaceSet();
            faces.reserve(2);    
        }

        public Edge(Vertex v1, Vertex v2)
        {
            ends = new Vertex[2];
            ends[0] = v1;  ends[1] = v2;
            label = 0;
            faces = new FaceSet();
            faces.reserve(2);
        }

        public Edge Clone()
        {
            Edge e = new Edge();
            e.ends = (Vertex[])ends.Clone();
            e.faces = faces.Clone();
            e.label = label;
            return e;
        }
        public void add(Face f)
        {
            faces.add(f);
        }

        public void remove(Face f)
        {
            faces.remove(f);
        }

        public double length()
        {
            return (ends[0].v - ends[1].v).norm();
        }

        public Vector UnitVector()
        {
            Vector n = (ends[0].v - ends[1].v);
            n.normalize();
            return n;
        }

        public bool boundary()
        {
            if (faces.getCount() < 2) return true;
            return false;
        }

        public bool orphan()
        {
            if (faces.getCount() == 0) return true;
            return false;
        }

        public void disconnect()
        {
            //only called if orphaned
            if (ends[0] != null) ends[0].edges.remove(this);
            if (ends[1] != null) ends[1].edges.remove(this);
            ends[0] = ends[1] = null;
            faces.clear();
        }

        public void replace(Vertex vold, Vertex vnew)
        {
            int i;
            for (i = 0; i < 2; i++)
            {
                if (ends[i] == vold)
                {
                    vold.edges.remove(this);
                    ends[i] = vnew;
                    vnew.edges.add(this);
                }
            }
        }

        public override bool Equals(object obj)
        {
            Edge e1 = this;
            Edge e2 = (Edge)obj;
            if (e1.ends[0].Equals(e2.ends[0]))
            {
                if (e1.ends[1].Equals(e2.ends[1]))
                    return true;
                else
                    return false;
            }
            else if (e1.ends[0].Equals(e2.ends[1]))
            {
                if (e1.ends[1].Equals(e2.ends[0]))
                    return true;
                else
                    return false;
            }
            else return false;
        }

        public Vertex[] ends;
        public FaceSet faces;
        public int label;

        public int id;
        public double l0;
        public double force;
    }
}
