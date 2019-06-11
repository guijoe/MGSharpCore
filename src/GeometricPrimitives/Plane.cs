

namespace MGSharp.Core.GeometricPrimitives
{
    public class Plane
    {
        public Vector Point;
        public Vector Normal;
        static float epsilon = 0.5f;

        public Plane()
        {
            Point = Vector.zero;
            Normal = Vector.up;
        }

        public Plane(Plane plane)
        {
            Point = plane.Point;
            Normal = plane.Normal;
        }

        public Plane(Vector point, Vector normal)
        {
            Point = point;
            Normal = normal;
        }

        public float LineIntersect(Vector lineStart, Vector lineEnd)
        {
            return (float) ((Normal * (Point - lineStart)) / (Normal * (lineEnd - lineStart)));
        }

        public float PointSide(Vector point)
        {
            return (float)(Normal * (point - Point));
        }

        public float PointSideNormalized(Vector pt)
        {
            Vector dif = (pt - Point);
            dif.normalize(); 
            return (float) (Normal * dif);
        }

        public Vector PointOrthogonalProjection(Vector PointA)
        {
            float t = (float) (Normal * (Point - PointA) / (Normal * Normal));
            return new Vector(
                                Normal.x * t + PointA.x,
                                Normal.y * t + PointA.y,
                                Normal.z * t + PointA.z
                              );
        }
        
        /*
        public static Plane RandomSphereSplittingPlane()
        {
            Vector n = Random.onUnitSphere;
            Plane P0 = new PlaneMath(Vector.zero, n);

            Vector Pt = Random.insideUnitSphere;
            while (Vector.Distance(Pt, P0.PointOrthogonalProjection(Pt)) > epsilon)
            {
                Pt = Random.insideUnitSphere;
            }
            return new Plane(Pt, n);
        }
        */
    }
}