using System.Collections.Generic;
using MGSharp.Core.MGModels;

namespace MGSharp.Core.GeometricPrimitives
{
    public class Vertex
    {
        public Vertex()
        {
            label = 0;
            pos = 0;
            dist = 1e15;
            next = null;
            edges = new EdgeSet();
            faces = new FaceSet();
            normal = new Vector();
            tex = new Vector();
            n = new Vector();
            v = new Vector();
        }
        public Vertex(Vector v)
        {
            this.v = v;
        }
        public Vertex Clone()
        {
            Vertex v = new Vertex();
            v.dist = dist;
            if(faces.getCount() > 0) v.faces = faces.Clone();
            v.label = label;
            v.n = (Vector)n.Clone();
            v.normal = (Vector)normal.Clone();
            v.next = next;
            v.pos = pos;
            v.tex = (Vector)tex.Clone();
            v.v = (Vector)this.v.Clone();
            v.wij = null;
            return v;
        }
        public int valence()
        {
            return edges.getCount();
        }
        public void ComputeNormal()
        {
            int i;
            normal = new Vector(0, 0, 0, 1);
            for (i = 0; i < faces.getCount(); i++)
            {
                normal += faces[i].normal;
            }
            normal.normalize();
        }

        public void ComputeNormal2()
        {
            int i;
            normal = new Vector(0, 0, 0, 1);
            for (i = 0; i < faces.getCount(); i++)
            {
                faces[i].ComputeNormal2();
                normal += faces[i].normal;
            }
            normal.normalize();
        }
        public bool adjacent(Vertex v)
        {
            int i;
            for (i = 0; i < edges.getCount(); i++)
            {
                if (edges[i].ends[0] == v) return true;
                if (edges[i].ends[1] == v) return true;
            }
            return false;
        }
        public void add(Face f)
        {
            faces.add(f);
        }
        public void remove(Face f)
        {
            faces.remove(f);
        }
        public void add(Edge e)
        {
            edges.add(e);
        }
        public void remove(Edge e)
        {
            edges.remove(e);
        }
        public bool boundary()
        {
            int i;
            for (i = 0; i < edges.getCount(); i++)
            {
                if (edges[i].boundary())
                    return true;
            }
            return false;
        }
        public override bool Equals(object o)
        {
            if (this == o) return true;
            else return false;
        }

        public void translate(Vector u)
        {
            v += u;
        }

        public void Move()
        {
            velocity = force * MGModel.damping;
            v += velocity * MGModel.dT;
        }
        
        public Vector GetPosition()
        {
            return v;
        }
        
        public int id { get; set; }
        public int cellId { get; set; }
        public Vector v { get; set; }
        //vertex normal (may also be used for other purposes)
        public Vector n;
        public Vector tex;
        //actual normal to be used
        public Vector normal;// { get; set; }
        public EdgeSet edges;
        public FaceSet faces;
        public int label;
        public int pos;
        //geodesic dist;
        public double dist;
        //weights for parameterisation
        // per vertex, so that weights on opposite
        // ends of an edge can differ
        public double[] wij;
        //for the creation of a path
        public Vertex next;

        //Vertex Dynamics
        public bool nullForces = false;
        public Vector velocity;// { get; set; }
        public Vector force;// { get; set; }
        public Vector externalForces;
        public Vector internalForces;
        public Vector nucleusForce0;
        public Vector nucleusForce1;
        public Vector globalForces;
        public float signedNucleusForce0;
        public List<int[]> externalNeighbours;
        public List<Vertex> internalNeighbours;
        public List<Vertex> collidingNeighbours;
    }
}