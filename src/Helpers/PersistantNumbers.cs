using System;

namespace MGSharp.Core.Helpers
{
    public class PersistantNumbers : Persistant
    {
        public ushort frame { get; set; }
        public ushort n1 { get; set; }
        public ushort n2 { get; set; }
        public ushort tissue { get; set; }

        public PersistantNumbers() { }

        public PersistantNumbers(int frame, int n1, int n2)
        {
            this.frame = (ushort)frame;
            this.n1 = (ushort)n1;
            this.n2 = (ushort)n2;
        }

        public PersistantNumbers(int frame, int n1, int n2, int tissue)
        {
            this.frame = (ushort) frame;
            this.n1 = (ushort)n1;
            this.n2 = (ushort)n2;
            this.tissue = (ushort)tissue;
        }
    }
}
