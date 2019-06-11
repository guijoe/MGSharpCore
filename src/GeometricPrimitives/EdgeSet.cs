using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MGSharp.Core.GeometricPrimitives
{
    public class EdgeSet : Set
    {
        public new Edge match(object x)
        {
            int i = index(x);
            if (i != -1) return (Edge)set[i];
            else return null;
        }

        public EdgeSet() : base()
        {
            
        }

        public EdgeSet(int s) : base()
        {
            size = s;
        }

        public new Edge this[int i]
        {
            get { return (Edge)set[i]; }
            set { set[i] = value; }
        }
        public EdgeSet Clone()
        {
            EdgeSet e = new EdgeSet();
            e.set = (object[])set.Clone();
            e.n = n;
            e.size = size;
            return e;
        }
    }
}
