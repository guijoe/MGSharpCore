using System;
using System.Collections.Generic;


namespace MGSharp.Core.GeometricPrimitives
{
	public class CustomHexagonal : Icosahedron42
    {
        public CustomHexagonal() : base(){ }
        
        public CustomHexagonal(int n) : base(n)
        {

        }

        public CustomHexagonal(int n, float r, float h) : base(n)
        {
            Default2Hexagonal(r, h);

            innerRadius = new Vector(r * Math.Sqrt(3) / 2, h / 2, r * 0.75);
        }
        
        public void Default2Hexagonal(float r, float h)
        {
            for (int i = 0; i < vertexCount(); i++)
            {
                double coef = 0f;
                double x = vertices[i].v.x;
                double z = vertices[i].v.z;

                if (Math.Abs(x) > 0 || Math.Abs(z) > 0)
                {
                    coef = r / (Math.Sqrt(x * x + z * z));
                }

                vertices[i].v.x = coef * vertices[i].v.x;
                vertices[i].v.z = coef * vertices[i].v.z;
                vertices[i].v.y = h * vertices[i].v.y;
            }
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