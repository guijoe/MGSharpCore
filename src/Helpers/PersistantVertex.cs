using MGSharp.Core.GeometricPrimitives;
using System;
using System.Text;

namespace MGSharp.Core.Helpers
{
    [Serializable()]
    public class PersistantVertex
    {
        //*
        public ushort tissueId;
        public ushort frame { get; set; }
        public ushort id { get; set; }
        public ushort cellId { get; set; }
        public Vector v { get; set; }
        //public float x { get; set; }
        //public float y { get; set; }
        //public float z { get; set; }
        //*/

        //public string str { get; set; }
        //public Byte[] strBytes { get; set; }

        public PersistantVertex() { }

        public PersistantVertex(int frame, int id, int cellId, Vector v)
        {
            ///*
            this.frame = (ushort) frame;
            this.id = (ushort) id;
            this.cellId = (ushort) cellId;
            //this.x = (float)v.x;
            //this.y = (float)v.y;
            //this.z = (float)v.z;
            this.v = v;
            //*/
            //strBytes = Encoding.ASCII.GetBytes(ToString(frame, id, cellId, v));
            //str = ToString(frame, id, cellId, v);
        }

        public PersistantVertex Clone()
        {
            
            PersistantVertex v = new PersistantVertex();
            //*
            v.v = (Vector)this.v.Clone();
            v.id = id;
            v.cellId = cellId;
            v.frame = frame;
            //*/
            return v;
            
        }

        public String ToString(int frame, int id, int cellId, Vector v)
        {
            return frame + ";" + id + ";" + cellId + ";" + v.ToString();
        }

        public String ToString()
        {
            return frame + ";" + id + ";" + cellId + ";" + v.ToString();
        }
    }
}
