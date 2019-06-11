using MGSharp.Core.Helpers;
using System;
using System.Collections.Generic;
using System.IO;

namespace MGSharp.Core.GeometricPrimitives
{
    public class Bounds2D
    {
        public Vector centre;
        public Vector size;
        public Vector extents;
        public Vector max;
        public Vector min;
        public Mesh convexHull;

        float deltaK = .0f;
        float damping = .1f;
        float lambda = 1f;
        float KPM = 100f;
        //public int resolution = 25;
        public Vector resolution = new Vector(25,25);
        //public int NConvolutions = 10;

        public Vector thau;
        public Vector[,] X;
        public double[,] vPHI;
        Vector[,] gradPHI;
        double[,] vPeronaMalik;
        public Vector[,] gradPeronaMalik;
        Vector[,] velocities;
        double[,] vGauss;

        public Bounds2D()
        {
            centre = new Vector();
            extents = new Vector();
            size = new Vector();
        }

        public Bounds2D(Vector min, Vector max)
        {
            this.min = min;
            this.max = max;

            this.size = max - min;
            this.extents = size / 2;
            //Console.WriteLine(size + ", " + extents);
            this.centre = (min + max) / 2;
        }

        public Bounds2D(Vector size)
        {

            this.size = size;
            this.extents = size / 2;

            this.min = -extents ;
            this.max = extents;

            //Console.WriteLine(size + ", " + extents);
            this.centre = (min + max) / 2;
        }

        public Bounds2D(Bounds bounds)
        {
            centre = bounds.centre;
            size = bounds.size;
            extents = bounds.extents;
            max = bounds.max;
            min = bounds.min;
            convexHull = new Mesh();
            convexHull.copy(bounds.convexHull);
        }

        public Vector[,] CreateSpatialMesh(int resolution)
        {
            Vector[,] X = new Vector[resolution,resolution];

            for (int i = 0; i < resolution; i++)
            {
                for (int j = 0; j < resolution; j++)
                {
                    X[i, j] = new Vector();
                    if (i == 0 || j == 0)
                    {
                        if (i == 0)
                        {
                            X[0, j].x = centre.x - extents.x + extents.x / resolution;
                            X[0, j].y = X[i, 0].y + j * size.y / resolution;
                        }
                        if (j == 0)
                        {
                            X[i, 0].y = centre.y - extents.y + extents.y / resolution;
                            X[i, 0].x = X[0, j].x + i * size.x / resolution;
                        }
                    }
                    else
                    {
                        X[i, j] = new Vector(
                            X[0, j].x + i * size.x / resolution,
                            X[i, 0].y + j * size.y / resolution
                        );
                    }
                    
                }
                //Console.WriteLine(X[i, 0]);
            }

            return X;
        }

        //*
        public Vector[,] CreateSpatialMesh(Vector resolution)
        {
            Vector[,] X = new Vector[(int) resolution.x, (int) resolution.y];

            for (int i = 0; i < resolution.x; i++)
            {
                for (int j = 0; j < resolution.y; j++)
                {
                    X[i, j] = new Vector();
                    if (i == 0 || j == 0)
                    {
                        if (i == 0)
                        {
                            X[0, j].x = centre.x - extents.x + extents.x / resolution.x;
                            X[0, j].y = X[i, 0].y + j * size.y / resolution.y;
                        }
                        if (j == 0)
                        {
                            X[i, 0].y = centre.y - extents.y + extents.y / resolution.x;
                            X[i, 0].x = X[0, j].x + i * size.x / resolution.y;
                        }
                    }
                    else
                    {
                        X[i, j] = new Vector(
                            X[0, j].x + i * size.x / resolution.x,
                            X[i, 0].y + j * size.y / resolution.y
                        );
                    }
                }
                //Console.WriteLine(X[i, 0]);
            }

            return X;
        }
        //*/

        /*
        public Vector[,] CreateSpatialMesh(Vector resolution)
        {
            Vector[,] X = new Vector[(int)resolution.x, (int)resolution.y];

            for (int i = 0; i < resolution.x; i++)
            {
                for (int j = 0; j < resolution.y; j++)
                {
                    X[i, j] = new Vector();
                    if (i == 0 || j == 0)
                    {
                        if (i == 0)
                        {
                            X[0, j].x = min.x + j;
                            X[0, j].y = min.y;
                        }
                        if (j == 0)
                        {
                            X[i, 0].y = min.y + i;
                            X[i, 0].x = min.x;
                        }
                    }
                    else
                    {
                        X[i, j] = new Vector(
                            min.x + j,
                            min.y + i
                        );
                    }
                }
            }

            return X;
        }
        */

        /*
        public void DiscretiseSpace(int resolution, Func<Vector, Bounds2D, double> PHI)
        {
            this.resolution = resolution;
            this.thau = size / resolution;
            this.thau = new Vector(1, 1, 1);

            X = new Vector[resolution, resolution];
            X = CreateSpatialMesh(resolution);


            double sum = 0;
            vPHI = new double[resolution + 1, resolution + 1];
            for (int i = 0; i < resolution; i++)
            {
                for (int j = 0; j < resolution; j++)
                {
                    vPHI[i, j] = PHI(X[i, j], this);
                    //sum += vPHI[i, j];
                }
            }
            //Console.WriteLine(sum);
        }
        */

        public void DiscretiseSpace(Vector resolution, Func<Vector, Bounds2D, double> PHI)
        {
            this.resolution = resolution;
            //this.thau = size / resolution;
            //this.thau = new Vector();
            //this.thau.x = size.x / resolution.x;
            //this.thau.y = size.y / resolution.y;

            //Console.WriteLine("Resolution: " + resolution);

            this.thau = new Vector(1, 1);

            X = new Vector[(int) resolution.x, (int)resolution.y];
            X = CreateSpatialMesh(resolution);


            double sum = 0;
            vPHI = new double[(int) resolution.x + 1, (int) resolution.y + 1];
            for (int i = 0; i < resolution.x; i++)
            {
                for (int j = 0; j < resolution.y; j++)
                {
                    vPHI[i, j] = PHI(X[i, j], this);
                    //sum += vPHI[i, j];
                }
            }
            //Console.WriteLine(sum);
        }

        public void ComputeGaussCoefficients()
        {
            vGauss = new double[3, 3];
            double sumGauss = 0;
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    vGauss[i, j] = Helper.Gauss2D(this,
                        new Vector(thau.x * (i - 1),
                            thau.y * (j - 1)
                        )
                    );
                    sumGauss += vGauss[i, j];
                    
                }
            }

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    vGauss[i, j] /= sumGauss;
                    //Console.WriteLine(i + "," + j + ":" + vGauss[i, j]);
                }
            }
        }

        /*
        public void GaussianConvolution(int NConvolutions)
        {
            ComputeGaussCoefficients();
            for (int t = 0; t < NConvolutions; t++)
            {
                for (int i = 0; i < resolution; i++)
                {
                    for (int j = 0; j < resolution; j++)
                    {
                        vPHI[i, j] *= vGauss[1, 1];
                        

                        List<int> neighbours = GetMooreNeighbourhood2D(i, j);
                        for (int l = 0; l < neighbours.Count; l++)
                        {
                            //int z = neighbours[l] / (resolution * resolution);
                            int xy = neighbours[l];//% (resolution * resolution);
                            int x = xy / resolution;
                            int y = xy % resolution;
                            vPHI[i, j] += vGauss[x - i + 1, y - j + 1] * vPHI[x, y];
                        }
                    }
                }
            }
        }
        */

        public void GaussianConvolution(int NConvolutions)
        {
            ComputeGaussCoefficients();
            for (int t = 0; t < NConvolutions; t++)
            {
                for (int i = 0; i < resolution.x; i++)
                {
                    for (int j = 0; j < resolution.y; j++)
                    {
                        vPHI[i, j] *= vGauss[1, 1];

                        List<int> neighbours = GetMooreNeighbourhood2D(i, j);
                        for (int l = 0; l < neighbours.Count; l++)
                        {
                            //int z = neighbours[l] / (resolution * resolution);
                            int xy = neighbours[l];
                            if(xy>=0 && xy < resolution.x * resolution.y)
                            {
                                int y = xy / (int)resolution.x;
                                int x = xy % (int)resolution.x;
                                vPHI[i, j] += vGauss[x - i + 1, y - j + 1] * vPHI[x, y];
                            }
                        }
                    }
                }
            }
        }

        public void ComputeVectorFields1(int resolution, double KPM, double deltaK)
        {
            string statesFile = "VectorFields.mg";
            statesFile = Simulator.logDir + "/" + statesFile;
            FileStream fs = new FileStream(statesFile, FileMode.OpenOrCreate);
            List<PersistantVectorFields> vectorFields = new List<PersistantVectorFields>();

            CsvSerializer<PersistantVectorFields> serializer = new CsvSerializer<PersistantVectorFields>();
            serializer.Separator = ';';
            

            velocities = new Vector[resolution, resolution];
            gradPHI = new Vector[resolution, resolution];
            vPeronaMalik = new double[resolution + 1, resolution + 1];
            gradPeronaMalik = new Vector[resolution, resolution];

            int count = 0;
            for (int i = 0; i < resolution; i++)
            {
                for (int j = 0; j < resolution; j++)
                {        
                    PersistantVectorFields pvf = new PersistantVectorFields();
                    pvf.I = i; pvf.J = j;
                    pvf.X = X[i, j];
                    pvf.vPHI = (float)vPHI[i,j];
                        

                    gradPHI[i, j] = new Vector(
                        (vPHI[i + 1, j] - vPHI[i, j]) / thau.x,
                        (vPHI[i, j + 1] - vPHI[i, j]) / thau.y
                    );
                    pvf.gradPHI = gradPHI[i,j];

                    vPeronaMalik[resolution, j] = vPeronaMalik[1, j];
                    vPeronaMalik[i, resolution] = vPeronaMalik[i, 1];
                    
                    vPeronaMalik[i, j] = Helper.PeronaMalik(gradPHI[i, j].sqrNorm(), KPM);
                    pvf.vPeronaMalik = (float)vPeronaMalik[i, j];
                    
                    vectorFields.Add(pvf);
                    count++;
                }
            }

            count = 0;
            for (int i = 0; i < resolution; i++)
            {
                for (int j = 0; j < resolution; j++)
                {   
                    gradPeronaMalik[i, j] = new Vector(
                        (vPeronaMalik[i + 1, j] - vPeronaMalik[i, j]) / thau.x,
                        (vPeronaMalik[i, j + 1] - vPeronaMalik[i, j]) / thau.y
                    );

                    vectorFields[count].gradPeronaMalik = gradPeronaMalik[i, j];

                    gradPeronaMalik[resolution - 1, j] = -gradPeronaMalik[0, j];
                    gradPeronaMalik[i, resolution - 1] = -gradPeronaMalik[i, 0];
                        
                    count++;
                }
            }

            for (int i = 0; i < resolution; i++)
            {
                for (int j = 0; j < resolution; j++)
                {
                    velocities[i, j] = -gradPeronaMalik[i, j];
                }
            }

            serializer.Serialize(fs, vectorFields);
            fs.Close();
        }

        /*
        public void ComputeVectorFields(int resolution, double KPM, double deltaK)
        {

            string statesFile = "VectorFields.mg";
            statesFile = Simulator.logDir + "/" + statesFile;
            FileStream fs = new FileStream(statesFile, FileMode.OpenOrCreate);
            List<PersistantVectorFields> vectorFields = new List<PersistantVectorFields>();

            CsvSerializer<PersistantVectorFields> serializer = new CsvSerializer<PersistantVectorFields>();
            serializer.Separator = ';';


            velocities = new Vector[resolution, resolution];
            gradPHI = new Vector[resolution, resolution];
            vPeronaMalik = new double[resolution + 1, resolution + 1];
            gradPeronaMalik = new Vector[resolution, resolution];

            int halfResolution = resolution / 2;

            int count = 0;
            for (int i = 0; i < resolution; i++)
            {
                for (int j = 0; j < resolution; j++)
                {
                    PersistantVectorFields pvf = new PersistantVectorFields();
                    pvf.I = i; pvf.J = j;
                    pvf.X = X[i, j];
                    pvf.vPHI = (float)vPHI[i, j];

                    gradPHI[i, j] = new Vector();
                    if (i == halfResolution)
                    {
                        gradPHI[i, j].x = 0;
                    }
                    else if (i < halfResolution)
                    {
                        gradPHI[i, j].x = (vPHI[i + 1, j] - vPHI[i, j]) / thau.x;
                    }
                    else
                    {
                        gradPHI[i, j].x = -(vPHI[i - 1, j] - vPHI[i, j]) / thau.x;
                    }
                    if (j == halfResolution)
                    {
                        gradPHI[i, j].y = 0;
                    }
                    else if (j < halfResolution)
                    {
                        gradPHI[i, j].y = (vPHI[i, j + 1] - vPHI[i, j]) / thau.y;
                    }
                    else
                    {
                        gradPHI[i, j].y = -(vPHI[i, j - 1] - vPHI[i, j]) / thau.y;
                    }

                    pvf.gradPHI = gradPHI[i, j];

                    vPeronaMalik[i, j] = Helper.PeronaMalik(gradPHI[i, j].sqrNorm(), KPM);
                    pvf.vPeronaMalik = (float)vPeronaMalik[i, j];
                    //Console.WriteLine(i +","+ j + "," + k + ":" + vPeronaMalik[i, j]);

                    vectorFields.Add(pvf);
                    count++;

                }
            }

            count = 0;
            for (int i = 0; i < resolution; i++)
            {
                for (int j = 0; j < resolution; j++)
                {
                    gradPeronaMalik[i, j] = new Vector();

                    if (i == halfResolution)
                    {
                        gradPeronaMalik[i, j].x = 0;
                    }
                    else if (i < halfResolution)
                    {
                        gradPeronaMalik[i, j].x = (vPeronaMalik[i + 1, j] - vPeronaMalik[i, j]) / thau.x;
                    }
                    else
                    {
                        gradPeronaMalik[i, j].x = -(vPeronaMalik[i - 1, j] - vPeronaMalik[i, j]) / thau.x;
                    }

                    if (j == halfResolution)
                    {
                        gradPeronaMalik[i, j].y = 0;
                    }
                    else if (j < halfResolution)
                    {
                        gradPeronaMalik[i, j].y = (vPeronaMalik[i, j + 1] - vPeronaMalik[i, j]) / thau.y;
                    }
                    else
                    {
                        gradPeronaMalik[i, j].y = -(vPeronaMalik[i, j - 1] - vPeronaMalik[i, j]) / thau.y;
                    }

                    vectorFields[count].gradPeronaMalik = gradPeronaMalik[i, j];
                    count++;
                }
            }

            for (int i = 0; i < resolution; i++)
            {
                for (int j = 0; j < resolution; j++)
                {
                    velocities[i, j] = -gradPeronaMalik[i, j];
                }
            }

            serializer.Serialize(fs, vectorFields);
            fs.Close();
        }
        */

        public void ComputeVectorFields(Vector resolution, double KPM, double deltaK)
        {
            string statesFile = "VectorFields.mg";
            statesFile = Simulator.logDir + "/" + statesFile;
            FileStream fs = new FileStream(statesFile, FileMode.OpenOrCreate);
            List<PersistantVectorFields> vectorFields = new List<PersistantVectorFields>();

            CsvSerializer<PersistantVectorFields> serializer = new CsvSerializer<PersistantVectorFields>();
            serializer.Separator = ';';


            velocities = new Vector[(int)resolution.x, (int)resolution.y];
            gradPHI = new Vector[(int)resolution.x, (int)resolution.y];
            vPeronaMalik = new double[(int)resolution.x + 1, (int)resolution.y + 1];
            gradPeronaMalik = new Vector[(int)resolution.x, (int)resolution.y];

            Vector halfResolution = resolution / 2;

            int count = 0;
            for (int i = 0; i < resolution.x; i++)
            {
                for (int j = 0; j < resolution.y; j++)
                {
                    PersistantVectorFields pvf = new PersistantVectorFields();
                    pvf.I = i; pvf.J = j;
                    pvf.X = X[i, j];
                    pvf.vPHI = (float)vPHI[i, j];

                    gradPHI[i, j] = new Vector();
                    if (i == halfResolution.x)
                    {
                        gradPHI[i, j].x = 0;
                    }else if (i < halfResolution.x)
                    {
                        gradPHI[i, j].x = (vPHI[i + 1, j] - vPHI[i, j]) / thau.x;
                    }
                    else
                    {
                        gradPHI[i, j].x = -(vPHI[i - 1, j] - vPHI[i, j]) / thau.x;
                    }
                    if(j== halfResolution.y)
                    {
                        gradPHI[i, j].y = 0;
                    }
                    else if (j < halfResolution.y)
                    {
                        gradPHI[i, j].y = (vPHI[i, j+1] - vPHI[i, j]) / thau.y;
                    }
                    else
                    {
                        gradPHI[i, j].y = -(vPHI[i, j-1] - vPHI[i, j]) / thau.y;
                    }
                    
                    pvf.gradPHI = gradPHI[i, j];

                    vPeronaMalik[i, j] = Helper.PeronaMalik(gradPHI[i, j].sqrNorm(), KPM);
                    pvf.vPeronaMalik = (float)vPeronaMalik[i, j];
                    //Console.WriteLine(i +","+ j + "," + k + ":" + vPeronaMalik[i, j]);

                    vectorFields.Add(pvf);
                    count++;
                    
                }
            }

            count = 0;
            for (int i = 0; i < resolution.x; i++)
            {
                for (int j = 0; j < resolution.y; j++)
                {
                    gradPeronaMalik[i, j] = new Vector();

                    if (i == halfResolution.x)
                    {
                        gradPeronaMalik[i, j].x = 0;
                    }
                    else if (i < halfResolution.x)
                    {
                        gradPeronaMalik[i, j].x = (vPeronaMalik[i + 1, j] - vPeronaMalik[i, j]) / thau.x;
                    }
                    else
                    {
                        gradPeronaMalik[i, j].x = -(vPeronaMalik[i - 1, j] - vPeronaMalik[i, j]) / thau.x;
                    }

                    if (j == halfResolution.y)
                    {
                        gradPeronaMalik[i, j].y = 0;
                    }
                    else if (j < halfResolution.y)
                    {
                        gradPeronaMalik[i, j].y = (vPeronaMalik[i, j+1] - vPeronaMalik[i, j]) / thau.y;
                    }
                    else
                    {
                        gradPeronaMalik[i, j].y = -(vPeronaMalik[i, j-1] - vPeronaMalik[i, j]) / thau.y;
                    }

                    vectorFields[count].gradPeronaMalik = gradPeronaMalik[i, j];
                    count++;
                }
            }

            for (int i = 0; i < resolution.x; i++)
            {
                for (int j = 0; j < resolution.y; j++)
                {
                    velocities[i, j] = -gradPeronaMalik[i, j];
                }
            }

            serializer.Serialize(fs, vectorFields);
            fs.Close();
        }
        
        /*
        public List<int> GetMooreNeighbourhood2D(int x, int y)
        {
            int dim = resolution;
            int dim2 = resolution * resolution;
            int size = dim2;

            int cellIndex = x * dim + y;

            int[] neighbours = {
                
                (cellIndex - dim - 1) ,
                (cellIndex - dim),
                (cellIndex - dim + 1),
                (cellIndex - 1),
                (cellIndex + 1),
                (cellIndex + dim - 1),
                (cellIndex + dim),
                (cellIndex + dim + 1),
            };

            List<int> mooreNeighbours = new List<int>();
            string neighboursString = cellIndex + ": ";

            for (int j = 0; j < 8; j++)
            {
                int i = j;
                int nxy = neighbours[j] % dim2;
                int nx = nxy / dim;
                int ny = nxy % dim;
                if (neighbours[j] >= 0 && neighbours[j] < size)
                {
                    if (i < 3 && nx == x - 1) { 
                        mooreNeighbours.Add(neighbours[j]); neighboursString += neighbours[j] + ",";
                    }
                    else if (i >= 3 && i < 5 && nx == x) { 
                        mooreNeighbours.Add(neighbours[j]); neighboursString += neighbours[j] + ",";
                    }
                    else if (i >= 5 && nx == x + 1) {
                        mooreNeighbours.Add(neighbours[j]); neighboursString += neighbours[j] + ",";
                    }
                }
            }
            //Console.WriteLine(neighboursString);
            return mooreNeighbours;
        }
        */

        public List<int> GetMooreNeighbourhood2D(int x, int y)
        {
            int dimX = (int) resolution.x;
            int dimY = (int) resolution.y;
            int dim2 = dimX * dimY;
            int size = dim2;

            int cellIndex = y * dimX + x;

            int[] neighbours = {
                (cellIndex - dimX - 1),
                (cellIndex - dimX),
                (cellIndex - dimX + 1),
                (cellIndex - 1),
                (cellIndex + 1),
                (cellIndex + dimX - 1),
                (cellIndex + dimX),
                (cellIndex + dimX + 1),
            };

            List<int> mooreNeighbours = new List<int>();
            //mooreNeighbours.AddRange(neighbours);

            //*
            string neighboursString = cellIndex + ": ";

            for (int j = 0; j < 8; j++)
            {
                int i = j;
                int nxy = neighbours[j] % dim2;
                int nx = nxy / dimX;
                int ny = nxy % dimX;
                if (neighbours[j] >= 0 && neighbours[j] < size)
                {
                    if (i < 3 && nx == y - 1)
                    {
                        mooreNeighbours.Add(neighbours[j]); neighboursString += neighbours[j] + ",";
                    }
                    else if (i >= 3 && i < 5 && nx == y)
                    {
                        mooreNeighbours.Add(neighbours[j]); neighboursString += neighbours[j] + ",";
                    }
                    else if (i >= 5 && nx == y + 1)
                    {
                        mooreNeighbours.Add(neighbours[j]); neighboursString += neighbours[j] + ",";
                    }
                }
            }
            //Console.WriteLine(neighboursString);
            //*/
            return mooreNeighbours;
        }

        /*
        public Vector DichotomicSearch(Vector point)
        {
            Vector min = X[0, 0], max = X[resolution - 1, resolution - 1];

            //Console.WriteLine(min);

            Vector min_point = point - min;
            int I = Helper.RoundToInt(min_point.x / thau.x);
            int J = Helper.RoundToInt(min_point.y / thau.y);
            
            return new Vector(I, J);
        }
        */

        public Vector DichotomicSearch(Vector point)
        {
            Vector min = X[0, 0], max = X[(int)resolution.x - 1, (int)resolution.y - 1];

            Vector min_point = point - min;
            int I = Helper.RoundToInt(min_point.x / thau.x);
            int J = Helper.RoundToInt(min_point.y / thau.y);

            return new Vector(I, J);
        }
    }
}