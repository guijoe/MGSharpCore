using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MGSharp.Core.GeometricPrimitives
{
    public class VertexSet : Set
    {
        public new Vertex match(object x)
        {
            int i = index(x);
            if (i != -1) return (Vertex)set[i];
            else return null;
        }

        public new Vertex this[int i]
        {
            get { return (Vertex)set[i]; }
            set { set[i] = value; }
        }

        public VertexSet Clone()
        {
            VertexSet v = new VertexSet();
            v.set = (object[])set.Clone();
            v.n = n;
            v.size = size;
            return v;
        }

        public int ContainsComponents(Vector u)
        {
            int i = 0;
            int index = -1;
            while (i<n)
            {
                Vertex vertex = (Vertex) set[i];

                if(u.x == vertex.v.x && u.y == vertex.v.y && u.z == vertex.v.z)
                {
                    index = vertex.pos;
                    i = n;
                }

                i++;
            }
            return index;
        }
    }
}