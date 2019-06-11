using MGSharp.Core.Helpers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MGSharp
{
    class Mathf
    {
        public const float E = 2.71828183f;
        public const float PI = 3.14159265f;
        public const float EPS = .0000001f;
        public static float epsilon = 0.000001f;

        public static List<PersistantVertex> floats;
        public static List<PersistantVertex> doubles;
        public static string logDir = "../../Logs/Mathf";

        static void Main(string[] args)
        {
            double timeStamp = (DateTime.Now.ToUniversalTime() - new DateTime(2000, 1, 1)).TotalMilliseconds;
            logDir += timeStamp.ToString();
            Directory.CreateDirectory(logDir);

            floats = new List<PersistantVertex>();
            doubles = new List<PersistantVertex>();
            /*
            float u1 = 0.02f;
            double r = Math.Exp(.5f * (Math.Log(2) + 2 * Math.Log(u1)));
            double a = Math.Log(r) - Math.Log(u1);

            double one = r * Math.Exp(a);
            double mOne = r * Math.Exp(-a);

            //Console.WriteLine($"{r}, {a}, {one}, {mOne}");
            Console.ReadLine();
            */

            float sinfX = 0f, xf=0;
            double sinX = 0, x=0;

            Stopwatch stopwatch = new Stopwatch();
            Console.WriteLine("Single precision ...");
            stopwatch.Start();

            Stopwatch s1 = new Stopwatch();
            for (int j=0; j<10000; j++)
            {
                
                for (int i = 0; i < 100; i++)
                {
                    xf = i * 36f * Mathf.PI / 180;
                    sinfX = Mathf.Exp(xf);
                }
                //s1.Stop();
                floats.Add(new PersistantVertex(j, 0, 0, new Core.GeometricPrimitives.Vector(xf, sinfX, s1.ElapsedMilliseconds)));
                //s1.Reset();
            }

            FileStream fs;
            string floatNumbersFile = "Floats.mg";
            floatNumbersFile = logDir + "/" + floatNumbersFile;
            fs = new FileStream(floatNumbersFile, FileMode.OpenOrCreate);
            CsvSerializer<PersistantVertex> serializer = new CsvSerializer<PersistantVertex>();
            serializer.Separator = ';';
            serializer.Serialize(fs, floats);
            fs.Close();

            stopwatch.Stop();
            Console.WriteLine("Single precision ended succesfully in " + stopwatch.ElapsedMilliseconds + " milliseconds.");
            Console.WriteLine();


            stopwatch.Reset();
            Console.WriteLine("Double precision ...");
            stopwatch.Start();

            Stopwatch s2 = new Stopwatch();
            for (int j=0; j<10000; j++)
            {
                s2.Start();
                for (int i = 0; i < 100; i++)
                {
                    x = i * 36 * Math.PI / 180;
                    sinX = Math.Exp(x);
                }
                //s1.Stop();
                doubles.Add(new PersistantVertex(j, 0, 0, new Core.GeometricPrimitives.Vector(x, sinX, s2.ElapsedMilliseconds)));
                //s1.Reset();
            }
            stopwatch.Stop();
            Console.WriteLine("Double precision ended succesfully in " + stopwatch.ElapsedMilliseconds + " milliseconds.");

            string doubleNumbersFile = "Doubles.mg";
            doubleNumbersFile = logDir + "/" + doubleNumbersFile;
            fs = new FileStream(doubleNumbersFile, FileMode.OpenOrCreate);
            serializer = new CsvSerializer<PersistantVertex>();
            serializer.Separator = ';';
            serializer.Serialize(fs, doubles);
            fs.Close();


            Console.ReadLine();
        }

        public static float Abs(float x)
        {
            if (x < 0f) return -x;
            else return x;
        }

        public static float Sin(float x)
        {
            //return x - x * x * x / 6f + x * x * x * x * x / 120f - x*x*x*x*x*x*x/5040f;
            //*
            bool first=true;
            int n = 1;
            float Fa = 1f;
            float y0 = 0f, y1 = 0f;
            while(Abs(y1 - y0) > EPS || first)
            {
                y0 = y1;
                y1 += Fa * Taylor(x, n);
                n += 2;
                Fa = -Fa;

                if (first) first = false;
            }
            return y1;
            //*/
        }

        public static float Exp(float x)
        {
            //return x - x * x * x / 6f + x * x * x * x * x / 120f - x*x*x*x*x*x*x/5040f;
            //*
            bool first = true;
            int n = 1;
            float y0 = 0f, y1 = 0f;
            while (Abs(y1 - y0) > EPS || first)
            {
                y0 = y1;
                y1 += Taylor(x, n);
                n++;
                //Fa = -Fa;

                if (first) first = false;
            }
            return y1;
            //*/
        }

        private static float Taylor(float x, int n)
        {
            //float y = 1;
            //if (x == 0) return 1f;
            if (n == 0) return 1f;
            else
            {
                return x / n * Taylor(x, n - 1);
            }
        }
    }
}
