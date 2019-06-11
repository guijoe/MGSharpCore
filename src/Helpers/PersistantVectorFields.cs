using MGSharp.Core.GeometricPrimitives;
using System;
using System.Text;

namespace MGSharp.Core.Helpers
{
    [Serializable()]
    public class PersistantVectorFields
    {
        //*
        public int I { get; set; }
        public int J { get; set; }
        public int K { get; set; }

        public Vector X { get; set; }
        public float vPHI { get; set; }
        public Vector gradPHI { get; set; }
        public float vPeronaMalik { get; set; }
        public Vector gradPeronaMalik { get; set; }
        
        public PersistantVectorFields() { }

        public PersistantVectorFields(int i, int j, int k, Vector X, double vPHI, Vector gradPHI, double vPeronaMalik, Vector gradPeronaMalik)
        {
            this.I = i;
            this.J = j;
            this.K = k;
            this.X = X;
            this.vPHI = (float)vPHI;
            this.gradPHI = gradPHI;
            this.vPeronaMalik = (float)vPeronaMalik;
            this.gradPeronaMalik = gradPeronaMalik;
        }

        

        public String ToString(int frame, int id, int cellId, Vector v)
        {
            return frame + ";" + id + ";" + cellId + ";" + v.ToString();
        }

        public String ToString()
        {
            //return frame + ";" + id + ";" + cellId + ";" + v.ToString();
            return "";
        }
    }
}
