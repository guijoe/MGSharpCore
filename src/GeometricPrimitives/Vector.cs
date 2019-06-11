using System;

namespace MGSharp.Core.GeometricPrimitives
{
    [Serializable()]
    public class Vector
    {
        public double x { get; set; }
        public double y { get; set; }
        public double z { get; set; }
        public double w;

        public static Vector zero = new Vector();
        public static Vector up = new Vector(0,1,0,1);
        public static Vector down = new Vector(0,-1,0,1);
        public static Vector left = new Vector(-1,0,0,1);
        public static Vector right = new Vector(1,0,0,1);

        public Vector()
        {
            x = 0.0; y = 0.0; z = 0.0; w = 1.0; 
        }
        
        public Vector(Vector v)
        {
            x = v.x; y = v.y; z = v.z; w = v.w;
        }

        public Vector(double x, double y)
        {
            this.x = x; this.y = y; this.z = 0; this.w = 1;
        }

        public Vector(double x, double y, double z)
        {
            this.x = x; this.y = y; this.z = z; this.w = 1;
        }

        public Vector(double x, double y, double z, double w)
        {
            this.x = x; this.y = y; this.z = z; this.w = w;
        }

        public object Clone()
        {
            Object o = new Vector(x, y, z, w);
            return o;
        }

        public override bool Equals(object b)
        {
            if ((((Vector)b).x == x) &&
            (((Vector)b).y == y) &&
            (((Vector)b).z == z) &&
            (((Vector)b).w == w)) return true;
            else return false;
        }

        public static Vector operator +(Vector a, Vector v)
        {
            return new Vector(a.x + v.x, a.y + v.y, a.z + v.z, a.w);
        }

        public static Vector operator -(Vector a, Vector v)
        {
            return new Vector(a.x - v.x, a.y - v.y, a.z - v.z, a.w);
        }

        //dot product
        public static double operator *(Vector a, Vector v)
        {
            return a.x * v.x + a.y * v.y + a.z * v.z;
        }

        public static Vector operator *(double c, Vector v)
        {
            return new Vector(c * v.x, c * v.y, c * v.z, v.w);
        }

        //cross product
        public static Vector operator ^(Vector a, Vector v)
        {
            Vector r = new Vector();

            r.x = a.y * v.z - a.z * v.y;
            r.y = a.z * v.x - a.x * v.z;
            r.z = a.x * v.y - a.y * v.x;
            r.w = a.w;
            return r;
        }

        //componentwise multiply
        public static Vector operator %(Vector a, Vector v)
        {
            Vector r = new Vector();
            r.x = a.x * v.x;
            r.y = a.y * v.y;
            r.z = a.z * v.z;
            r.w = a.w * v.w;
            return r;
        }

        public static Vector operator *(Vector a, double c)
        {
            return new Vector(c * a.x, c * a.y, c * a.z, a.w);
        }

        public static Vector operator /(Vector a, double c)
        {
            return new Vector(a.x / c, a.y / c, a.z / c, a.w);
        }

        //unary minus
        public static Vector operator -(Vector v)
        {
            return new Vector(-v.x, -v.y, -v.z, v.w);
        }

        //unary plus
        public static Vector operator +(Vector v)
        {
            return v;
        }

        public double norm()
        {
            return Math.Sqrt(x * x + y * y + z * z);
        }

        public double sqrNorm()
        {
            return x * x + y * y + z * z;
        }

        public void normalize()
        {
            double s = norm();
            if (s != 0.0) { x /= s; y /= s; z /= s; }
        }
        
        public void Copy(Vector v) {
            x = v.x;
            y = v.y;
            z = v.z;
            w = v.w;
        }

        public static Vector reflect(Vector v, Vector n)
        {
            Vector r;
            r = 2.0 * (v * n) * n - v;
            return r;
        }

        public static Vector transmit(Vector i, Vector n, double n1, double n2)
        {
            Vector t;
            double w, r, k;

            i = (Vector)i.Clone();
            n = (Vector)n.Clone();
            i.normalize();
            n.normalize();
            if (n * i > 0)
            {
                n = -n; r = n1;
                n1 = n2; n2 = r;
            }
            r = n1 / n2;
            w = -(i * n) * r;
            k = 1.0 + (w - r) * (w + r);
            //check for total internal reflection
            if (k < 0.0)
            {
                return reflect(-i, n);
            }
            k = Math.Sqrt(k);
            t = r * i + (w - k) * n;
            t.normalize();
            return t;
        }

        public static double Distance(Vector a, Vector v)
        {
            return (a-v).norm();
        }

        public static double SqrDistance(Vector a, Vector v)
        {
            return (a.x - v.x) * (a.x - v.x) + (a.y - v.y) * (a.y - v.y) + (a.z - v.z)* (a.z - v.z);
        }

        public static Vector Lerp(Vector a, Vector b, float t)
        {
            Vector barycentre = new Vector();
            barycentre = a + (b - a) * t;
            return barycentre;
        }

        //*
        public override string ToString()
        {
            //return x + ", " + y + ", " + z;
            return String.Format("{0:N4},{1:N4},{2:N4}", x, y, z);
        }
        //*/
    }
}
