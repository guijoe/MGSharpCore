using MGSharp.Core.GeometricPrimitives;
using System;
using System.Text;

namespace MGSharp.Core.Helpers
{
    public class PersistantMetrics
    {
        public int frame { get; set; }
        public float elasticEnergy { get; set; }

        public PersistantMetrics() { }

        public PersistantMetrics(int frame, float E)
        {
            this.frame = frame;
            this.elasticEnergy = E;
        }
    }
}
