using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MGSharp.Core.GeometricPrimitives
{
    public class FaceSet : Set
    {
        public new Face match(object x)
        {
            int i = index(x);
            if (i != -1) return (Face)set[i];
            else return null;
        }
        public new Face this[int i]
        {
            get { return (Face)set[i]; }
            set { set[i] = value; }
        }
        public FaceSet Clone()
        {
            FaceSet f = new FaceSet();
            f.set = (object[])set.Clone();
            f.n = n;
            f.size = size;
            return f;
        }
    }
}
