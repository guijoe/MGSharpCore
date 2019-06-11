using System;
using System.Collections.Generic;
using MGSharp.Core.GeometricPrimitives;

namespace MGSharp.Core.MGModels
{
    public class MGModel
    {
        public static float epsilon = 0.000001f;

        //public static int nbOfParticles;
        public static float Rcell { get { return 1f; } set { } }
        public static Vector innerRcell;
        public static int neighbourPeriod = 100;
        public static float maximumNeighbourDistance;// { get { return 2 * Rcell; } }
        public static bool searchNeighbours = false;
        public static int nextNeighboursSearchFrame = 0;
        public static List<Mesh> modelMeshes;
        public static float dT {get; set;}
        public static float delta { get; set; }
        public static float T { get; set; }
        public static float rho { get; set; }
        public static float alpha { get; set; }
        public static float u0 { get; set; }
        public static float u1 { get; set; }
        public static float springForce { get; set; }
        public static float damping { get; set; }
        public static float DCol { get; set; }
        public static float DInt { get; set; }
        public static float divisionRate { get; set; }
        public static int cellCyclePeriod { get; set; }
        public static bool elasticExternalSpring { get; set; }
        public static bool staticNeighbourhood { get; set; }
        public static bool mooreNeighbourhoodForCells { get; set; }
        public static bool hexagonalNeighbourhoodForCells { get; set; }
        public static int maxExternalNeighbours { get; set; }
        
        public static float[,] DInteraction;
        public static float[,] J;

        public static double ForceAdj(float d, float u)
        {
            if (d == 0)
            {
                return 0f;
            }

            float exponent = rho * (DCol - d);
            return u * rho * (2 * Math.Exp(2 * exponent) - alpha * Math.Exp(exponent));
        }

        public static double Force(double d, float u, double cellReq)
        {
            if (d == cellReq)
            {
                return 0f;
            }
            
            double exponent = rho * (cellReq - d);
            return u * rho * (2f * Math.Exp(2 * exponent) - alpha * Math.Exp(exponent));
        }

        public static double NucleusToMembraneForce(double d, float u, double cellReqNucleus)
        {
            if (Math.Abs(d - cellReqNucleus) < MGModel.epsilon)
            {
                return 0f;
            }
            
            double exponent = rho * (cellReqNucleus - d);
            return u * rho * (2f * Math.Exp(2 * exponent) - alpha * Math.Exp(exponent));
        }

        public static float[,] ScaleInteractionDistances(float coef)
        {
            int dim = (int)Math.Sqrt(DInteraction.Length);
            for (int i = 0; i < dim; i++)
            {
                for (int j = 0; j < dim; j++)
                {
                    DInteraction[i, j] *= coef;
                }
            }
            return DInteraction;
        }
    }
}