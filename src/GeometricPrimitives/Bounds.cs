using MGSharp.Core.Helpers;
using System;
using System.Collections.Generic;
using System.IO;

namespace MGSharp.Core.GeometricPrimitives
{
    public class Bounds
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
        public int resolution = 25;
        //public int NConvolutions = 10;
        
        public Vector thau;
        Vector[,,] X;
        public double[,,] vPHI;
        Vector[,,] gradPHI;
        double[,,] vPeronaMalik;
        public Vector[,,] gradPeronaMalik;
        Vector[,,] velocities;
        double[,,] vGauss;

        public Bounds()
        {
            centre = new Vector();
            extents = new Vector();
            size = new Vector();
        }

        public Bounds(Bounds bounds)
        {
            centre = bounds.centre;
            size = bounds.size;
            extents = bounds.extents;
            max = bounds.max;
            min = bounds.min;
            convexHull = new Mesh();
            convexHull.copy(bounds.convexHull);
        }

        public Bounds(Mesh[] meshes)
        {

        }

        public Vector[,,] CreateSpatialMesh(int resolution)
        {
            Vector[,,] X = new Vector[resolution,resolution,resolution];
            for (int i = 0; i < resolution; i++)
            {
                for (int j = 0; j < resolution; j++)
                {
                    for (int k = 0; k < resolution; k++)
                    {
                        X[i, j, k] = new Vector();
                        if (i == 0 || j == 0 || k == 0)
                        {
                            if (i == 0)
                            {
                                X[0, j, k].x = centre.x - extents.x + extents.x / resolution;
                                X[0, j, k].y = X[i, 0, k].y + j * size.y / resolution;
                                X[0, j, k].z = X[i, j, 0].z + k * size.z / resolution;
                            }
                            if (j == 0)
                            {
                                X[i, 0, k].y = centre.y - extents.y + extents.y / resolution;
                                X[i, 0, k].x = X[0, j, k].x + i * size.x / resolution;
                                X[i, 0, k].z = X[i, j, 0].z + k * size.z / resolution;
                            }
                            if (k == 0)
                            {
                                X[i, j, 0].z = centre.z - extents.z + extents.z / resolution;
                                X[i, j, 0].x = X[0, j, k].x + i * size.x / resolution;
                                X[i, j, 0].y = X[i, 0, k].y + j * size.y / resolution;
                            }
                        }
                        else
                        {
                            X[i, j, k] = new Vector(
                                X[0, j, k].x + i * size.x / resolution,
                                X[i, 0, k].y + j * size.y / resolution,
                                X[i, j, 0].z + k * size.z / resolution
                            );
                        }
                    }
                }
            }
            return X;
        }

        public void DiscretiseSpace(int resolution, Func<Vector, Bounds, double> PHI)
        {
            this.resolution = resolution;
            this.thau = size / resolution;

            X = new Vector[resolution, resolution, resolution];
            X = CreateSpatialMesh(resolution);


            double sum = 0;
            vPHI = new double[resolution + 1, resolution + 1, resolution + 1];
            for (int i = 0; i < resolution; i++)
            {
                for (int j = 0; j < resolution; j++)
                {
                    for (int k = 0; k < resolution; k++)
                    {
                        vPHI[i, j, k] = PHI(X[i, j, k], this);
                        //sum += vPHI[i, j, k];
                        
                    }
                }
            }
            //Console.WriteLine(sum);
        }

        public void ComputeGaussCoefficients()
        {
            vGauss = new double[3, 3, 3];
            double sumGauss = 0;
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        vGauss[i, j, k] = Helper.Gauss(this,
                            new Vector(thau.x * (i - 1),
                                thau.y * (j - 1),
                                thau.z * (k - 1)
                            )
                        );
                        sumGauss += vGauss[i, j, k];
                    }
                }
            }

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        vGauss[i, j, k] /= sumGauss;
                        //Console.WriteLine(i + "," + j + "," + k + ":" + vGauss[i, j, k]);
                    }
                }
            }
        }

        public void GaussianConvolution(int NConvolutions)
        {
            ComputeGaussCoefficients();
            for (int t = 0; t < NConvolutions; t++)
            {
                for (int i = 0; i < resolution; i++)
                {
                    for (int j = 0; j < resolution; j++)
                    {
                        for (int k = 0; k < resolution; k++)
                        {
                            vPHI[i, j, k] *= vGauss[1, 1, 1];

                            List<int> neighbours = GetMooreNeighbourhood3D(i, j, k);
                            for (int l = 0; l < neighbours.Count; l++)
                            {
                                int z = neighbours[l] / (resolution * resolution);
                                int xy = neighbours[l] % (resolution * resolution);
                                int x = xy / resolution;
                                int y = xy % resolution;
                                vPHI[i, j, k] += vGauss[x - i + 1, y - j + 1, z - k + 1] * vPHI[x, y, z];
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
            

            velocities = new Vector[resolution, resolution, resolution];
            gradPHI = new Vector[resolution, resolution, resolution];
            vPeronaMalik = new double[resolution + 1, resolution + 1, resolution + 1];
            gradPeronaMalik = new Vector[resolution, resolution, resolution];

            int count = 0;
            for (int i = 0; i < resolution; i++)
            {
                for (int j = 0; j < resolution; j++)
                {
                    for (int k = 0; k < resolution; k++)
                    {
                        PersistantVectorFields pvf = new PersistantVectorFields();
                        pvf.I = i; pvf.J = j; pvf.K = k;
                        pvf.X = X[i, j, k];
                        pvf.vPHI = (float)vPHI[i,j,k];
                        

                        gradPHI[i, j, k] = new Vector(
                            (vPHI[i + 1, j, k] - vPHI[i, j, k]) / thau.x,
                            (vPHI[i, j + 1, k] - vPHI[i, j, k]) / thau.y,
                            (vPHI[i, j, k + 1] - vPHI[i, j, k]) / thau.z
                        );
                        pvf.gradPHI = gradPHI[i,j,k];

                        vPeronaMalik[resolution, j, k] = vPeronaMalik[1, j, k];
                        vPeronaMalik[i, resolution, k] = vPeronaMalik[i, 1, k];
                        vPeronaMalik[i, j, resolution] = vPeronaMalik[i, j, 1];

                        vPeronaMalik[i, j, k] = Helper.PeronaMalik(gradPHI[i, j, k].sqrNorm(), KPM);
                        pvf.vPeronaMalik = (float)vPeronaMalik[i, j, k];
                        //Console.WriteLine(i +","+ j + "," + k + ":" + vPeronaMalik[i, j, k]);

                        vectorFields.Add(pvf);
                        count++;
                    }
                }
            }

            count = 0;
            for (int i = 0; i < resolution; i++)
            {
                for (int j = 0; j < resolution; j++)
                {
                    for (int k = 0; k < resolution; k++)
                    {
                        gradPeronaMalik[i, j, k] = new Vector(
                            (vPeronaMalik[i + 1, j, k] - vPeronaMalik[i, j, k]) / thau.x,
                            (vPeronaMalik[i, j + 1, k] - vPeronaMalik[i, j, k]) / thau.y,
                            (vPeronaMalik[i, j, k + 1] - vPeronaMalik[i, j, k]) / thau.z
                        );

                        vectorFields[count].gradPeronaMalik = gradPeronaMalik[i, j, k];

                        gradPeronaMalik[resolution - 1, j, k] = -gradPeronaMalik[0, j, k];
                        gradPeronaMalik[i, resolution - 1, k] = -gradPeronaMalik[i, 0, k];
                        gradPeronaMalik[i, j, resolution - 1] = -gradPeronaMalik[i, j, 0];

                        count++;
                    }
                }
            }

            for (int i = 0; i < resolution; i++)
            {
                for (int j = 0; j < resolution; j++)
                {
                    for (int k = 0; k < resolution; k++)
                    {
                        velocities[i, j, k] = -gradPeronaMalik[i, j, k];
                    }
                }
            }

            serializer.Serialize(fs, vectorFields);
            fs.Close();
        }

        public void ComputeVectorFields(int resolution, double KPM, double deltaK)
        {

            string statesFile = "VectorFields.mg";
            statesFile = Simulator.logDir + "/" + statesFile;
            FileStream fs = new FileStream(statesFile, FileMode.OpenOrCreate);
            List<PersistantVectorFields> vectorFields = new List<PersistantVectorFields>();

            CsvSerializer<PersistantVectorFields> serializer = new CsvSerializer<PersistantVectorFields>();
            serializer.Separator = ';';


            velocities = new Vector[resolution, resolution, resolution];
            gradPHI = new Vector[resolution, resolution, resolution];
            vPeronaMalik = new double[resolution + 1, resolution + 1, resolution + 1];
            gradPeronaMalik = new Vector[resolution, resolution, resolution];

            int halfResolution = resolution / 2;

            int count = 0;
            for (int i = 0; i < resolution; i++)
            {
                for (int j = 0; j < resolution; j++)
                {
                    for (int k = 0; k < resolution; k++)
                    {
                        PersistantVectorFields pvf = new PersistantVectorFields();
                        pvf.I = i; pvf.J = j; pvf.K = k;
                        pvf.X = X[i, j, k];
                        pvf.vPHI = (float)vPHI[i, j, k];

                        gradPHI[i, j, k] = new Vector();
                        if (i == halfResolution)
                        {
                            gradPHI[i, j, k].x = 0;
                        }else if (i < halfResolution)
                        {
                            gradPHI[i, j, k].x = (vPHI[i + 1, j, k] - vPHI[i, j, k]) / thau.x;
                        }
                        else
                        {
                            gradPHI[i, j, k].x = -(vPHI[i - 1, j, k] - vPHI[i, j, k]) / thau.x;
                        }
                        if(j== halfResolution)
                        {
                            gradPHI[i, j, k].y = 0;
                        }
                        else if (j < halfResolution)
                        {
                            gradPHI[i, j, k].y = (vPHI[i, j+1, k] - vPHI[i, j, k]) / thau.y;
                        }
                        else
                        {
                            gradPHI[i, j, k].y = -(vPHI[i, j-1, k] - vPHI[i, j, k]) / thau.y;
                        }
                        if (k == halfResolution)
                        {
                            gradPHI[i, j, k].z = 0;
                        }
                        else if (k < halfResolution)
                        {
                            gradPHI[i, j, k].z = (vPHI[i, j, k+1] - vPHI[i, j, k]) / thau.z;
                        }
                        else
                        {
                            gradPHI[i, j, k].z = -(vPHI[i, j, k-1] - vPHI[i, j, k]) / thau.z;
                        }

                        //gradPHI[i, j, k] = new Vector(
                        //    (vPHI[i + 1, j, k] - vPHI[i, j, k]) / thau.x,
                        //    (vPHI[i, j + 1, k] - vPHI[i, j, k]) / thau.y,
                        //    (vPHI[i, j, k + 1] - vPHI[i, j, k]) / thau.z
                        //);
                        pvf.gradPHI = gradPHI[i, j, k];

                        vPeronaMalik[i, j, k] = Helper.PeronaMalik(gradPHI[i, j, k].sqrNorm(), KPM);
                        pvf.vPeronaMalik = (float)vPeronaMalik[i, j, k];
                        //Console.WriteLine(i +","+ j + "," + k + ":" + vPeronaMalik[i, j, k]);

                        vectorFields.Add(pvf);
                        count++;
                    }
                }
            }

            count = 0;
            for (int i = 0; i < resolution; i++)
            {
                for (int j = 0; j < resolution; j++)
                {
                    for (int k = 0; k < resolution; k++)
                    {
                        gradPeronaMalik[i, j, k] = new Vector();

                        if (i == halfResolution)
                        {
                            gradPeronaMalik[i, j, k].x = 0;
                        }
                        else if (i < halfResolution)
                        {
                            gradPeronaMalik[i, j, k].x = (vPeronaMalik[i + 1, j, k] - vPeronaMalik[i, j, k]) / thau.x;
                        }
                        else
                        {
                            gradPeronaMalik[i, j, k].x = -(vPeronaMalik[i - 1, j, k] - vPeronaMalik[i, j, k]) / thau.x;
                        }

                        if (j == halfResolution)
                        {
                            gradPeronaMalik[i, j, k].y = 0;
                        }
                        else if (j < halfResolution)
                        {
                            gradPeronaMalik[i, j, k].y = (vPeronaMalik[i, j+1, k] - vPeronaMalik[i, j, k]) / thau.y;
                        }
                        else
                        {
                            gradPeronaMalik[i, j, k].y = -(vPeronaMalik[i, j-1, k] - vPeronaMalik[i, j, k]) / thau.y;
                        }

                        if (k == halfResolution)
                        {
                            gradPeronaMalik[i, j, k].z = 0;
                        }
                        else if (k < halfResolution)
                        {
                            gradPeronaMalik[i, j, k].z = (vPeronaMalik[i, j, k+1] - vPeronaMalik[i, j, k]) / thau.z;
                        }
                        else
                        {
                            gradPeronaMalik[i, j, k].z = -(vPeronaMalik[i, j, k-1] - vPeronaMalik[i, j, k]) / thau.z;
                        }

                        //gradPeronaMalik[i, j, k] = new Vector(
                        //    (vPeronaMalik[i + 1, j, k] - vPeronaMalik[i, j, k]) / thau.x,
                        //    (vPeronaMalik[i, j + 1, k] - vPeronaMalik[i, j, k]) / thau.y,
                        //    (vPeronaMalik[i, j, k + 1] - vPeronaMalik[i, j, k]) / thau.z
                        //);

                        vectorFields[count].gradPeronaMalik = gradPeronaMalik[i, j, k];
                        count++;
                    }
                }
            }

            for (int i = 0; i < resolution; i++)
            {
                for (int j = 0; j < resolution; j++)
                {
                    for (int k = 0; k < resolution; k++)
                    {
                        velocities[i, j, k] = -gradPeronaMalik[i, j, k];
                    }
                }
            }

            serializer.Serialize(fs, vectorFields);
            fs.Close();
        }

        public List<int> GetMooreNeighbourhood3D(int x, int y, int z)
        {
            int dim = resolution;
            int dim2 = resolution * resolution;
            int maxSize = resolution * resolution * resolution;
            int size = maxSize;

            int cellIndex = z * dim2 + x * dim + y;

            int[] neighbours = {
                (cellIndex - dim2 - dim - 1),
                (cellIndex - dim2 - dim),
                (cellIndex - dim2 - dim + 1),
                (cellIndex - dim2 - 1),
                (cellIndex - dim2),
                (cellIndex - dim2 + 1),
                (cellIndex - dim2 + dim - 1),
                (cellIndex - dim2 + dim),
                (cellIndex - dim2 + dim + 1),
                (cellIndex - dim - 1) ,
                (cellIndex - dim),
                (cellIndex - dim + 1),
                (cellIndex - 1),
                (cellIndex + 1),
                (cellIndex + dim - 1),
                (cellIndex + dim),
                (cellIndex + dim + 1),
                (cellIndex + dim2 - dim - 1),
                (cellIndex + dim2 - dim),
                (cellIndex + dim2 - dim + 1),
                (cellIndex + dim2 - 1),
                (cellIndex + dim2),
                (cellIndex + dim2 + 1),
                (cellIndex + dim2 + dim - 1),
                (cellIndex + dim2 + dim),
                (cellIndex + dim2 + dim + 1),
            };

            List<int> mooreNeighbours = new List<int>();

            for (int j = 0; j < 9; j++)
            {
                int i = j;
                int nz = neighbours[j] / dim2;
                int nxy = neighbours[j] % dim2;
                int nx = nxy / dim;
                int ny = nxy % dim;
                if (neighbours[j] >= 0 && neighbours[j] < size)
                {
                    if (i < 3 && nz == z - 1 && nx == x - 1)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 3 && i < 6 && nz == z - 1 && nx == x)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 6 && nz == z - 1 && nx == x + 1)
                        mooreNeighbours.Add(neighbours[j]);
                }
            }
            for (int j = 9; j < 17; j++)
            {
                int i = j - 9;
                int nz = neighbours[j] / dim2;
                int nxy = neighbours[j] % dim2;
                int nx = nxy / dim;
                int ny = nxy % dim;
                if (neighbours[j] >= 0 && neighbours[j] < size)
                {
                    if (i < 3 && nz == z && nx == x - 1)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 3 && i < 5 && nz == z && nx == x)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 5 && nz == z && nx == x + 1)
                        mooreNeighbours.Add(neighbours[j]);
                }
            }
            for (int j = 17; j < 26; j++)
            {
                int i = j - 17;
                int nz = neighbours[j] / dim2;
                int nxy = neighbours[j] % dim2;
                int nx = nxy / dim;
                int ny = nxy % dim;
                if (neighbours[j] >= 0 && neighbours[j] < size)
                {
                    if (i < 3 && nz == z + 1 && nx == x - 1)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 3 && i < 6 && nz == z + 1 && nx == x)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 6 && nz == z + 1 && nx == x + 1)
                        mooreNeighbours.Add(neighbours[j]);
                }
            }

            return mooreNeighbours;
        }

        public Vector DichotomicSearch(Vector point)
        {
            Vector min = X[0, 0, 0], max = X[resolution - 1, resolution - 1, resolution - 1];

            Vector min_point = point - min;
            int I = Helper.RoundToInt(min_point.x / thau.x);
            int J = Helper.RoundToInt(min_point.y / thau.y);
            int K = Helper.RoundToInt(min_point.z / thau.z);
            
            return new Vector(I, J, K);
        }

        public double Kronecker(Vector u, Bounds bds)
        {
            double val = 0;
            if (convexHull.isPointInConvexMesh(u))
            {
                val = 1f;
            }
            return val;
        }
    }
}
