using System;
using MGSharp.Core.MGCellModels;
using MGSharp.Core.GeometricPrimitives;
using MGSharp.Core.MGModels;
using System.Threading.Tasks;
using System.Collections.Generic;
using MGSharp.Core.Helpers;
using System.Text;

namespace MGSharp.Core.MGCellPopulation
{
    class CellPopulation
    {
        public int populationSize;
        public int maxPopulationSize;
        public List<Tissue> tissues;
        public MGCell[] cells;
        public int[] sigma;
        protected Vector dim;
        public Vector halfDim;
        int appliedForces;
        int frame = 0;
        StringBuilder interCellEdges;
        public Vector reference = new Vector();
        public Bounds bounds;
        List<Vector> subset;
        public string name;

        public CellPopulation(){}

        public CellPopulation(int popSize, int popMaxSize)
        {
            CreatePopulation(popSize, popMaxSize);
        }

        public CellPopulation(int popSize, Vector dim)
        {
            //interCellEdges = new StringBuilder();
            //interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            Simulator.interCellEdges = new StringBuilder();
            Simulator.interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            populationSize = popSize;
            maxPopulationSize = (int)(dim.x * dim.y * dim.z);
            cells = new MGCell[maxPopulationSize];
            sigma = new int[maxPopulationSize];
            this.dim = dim;
            halfDim = dim / 2;

            
            for (int i = 0; i < populationSize; i++)
            {
                cells[i] = new MGCell(i);
            }
        }

        public CellPopulation(int popSize, int popMaxSize, List<Tissue> tissueList)
        {
            //interCellEdges = new StringBuilder();
            //interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            Simulator.interCellEdges = new StringBuilder();
            Simulator.interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            populationSize = popSize;
            maxPopulationSize = popMaxSize;
            cells = new MGCell[popMaxSize];
            sigma = new int[popMaxSize];
            double d = Math.Round(Math.Pow(maxPopulationSize, 1f / 3f));
            dim = new Vector(d, d, d);
            halfDim = dim / 2;

            tissues = tissueList;
            int start = 0;
            for (int i = 0; i < tissueList.Count; i++)
            {
                tissues[i].indices = new int[tissues[i].maxPopulationSize];

                start += i == 0 ? 0 : tissueList[i - 1].populationSize;
                for (int j = 0; j < tissueList[i].populationSize; j++)
                {
                    //int k = i == 0 ? 0 : i - 1;
                    //int index = i * tissueList[k].populationSize + j;

                    int index = start + j;
                    cells[index] = new MGCell(index, tissueList[i].mesh);
                    cells[index].tissueName = tissueList[i].name;
                    tissues[i].indices[j] = index;
                    tissues[i].cells[j] = cells[index];
                    sigma[index] = index;
                }
            }
        }

        public CellPopulation(int popSize, int popMaxSize, List<Tissue> tissueList, bool randomCycleTime)
        {
            //interCellEdges = new StringBuilder();
            //interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            Simulator.interCellEdges = new StringBuilder();
            Simulator.interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            populationSize = popSize;
            maxPopulationSize = popMaxSize;
            //Console.WriteLine(popMaxSize);
            cells = new MGCell[popMaxSize];
            sigma = new int[popMaxSize];
            double d = Math.Round(Math.Pow(maxPopulationSize, 1f / 3f));
            dim = new Vector(d, d, d);
            halfDim = dim / 2;

            tissues = tissueList;
            int start = 0;
            for (int i = 0; i < tissueList.Count; i++)
            {
                tissues[i].indices = new int[tissues[i].maxPopulationSize];

                start += i == 0 ? 0 : tissueList[i - 1].populationSize;
                for (int j = 0; j < tissueList[i].populationSize; j++)
                {
                    int index = start + j;

                    cells[index] = new MGCell(index, tissueList[i].mesh);
                    cells[index].tissueId = i;
                    cells[index].tissueName = tissueList[i].name;

                    tissues[i].indices[j] = index;
                    tissues[i].cells[j] = cells[index];
                    sigma[index] = index;
                    cells[index].cycleTime = new Random().Next(MGModel.cellCyclePeriod);
                }
            }
        }

        public CellPopulation(int popMaxSize, List<Tissue> tissueList)
        {
            Simulator.interCellEdges = new StringBuilder();
            Simulator.interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            populationSize = 0;
            maxPopulationSize = popMaxSize;
            cells = new MGCell[popMaxSize];
            sigma = new int[popMaxSize];
            double d = Math.Round(Math.Pow(maxPopulationSize, 1f / 3f));
            dim = new Vector(d, d, d);
            halfDim = dim / 2;

            tissues = tissueList;
            int start = 0;
            for (int i = 0; i < tissueList.Count; i++)
            {
                tissues[i].indices = new int[tissues[i].maxPopulationSize];

                start += i == 0 ? 0 : tissueList[i - 1].populationSize;
                for (int j = 0; j < tissueList[i].populationSize; j++)
                {
                    int index = start + j;
                    cells[index] = new MGCell(index, tissueList[i].cells[j]);
                    cells[index].tissueName = tissueList[i].name;
                    tissues[i].indices[j] = index;
                    tissues[i].cells[j] = cells[index];
                    sigma[index] = index;
                }
                populationSize += tissueList[i].populationSize;
            }
        }

        public CellPopulation(string folder)
        {
            //interCellEdges = new StringBuilder();
            //interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            Simulator.interCellEdges = new StringBuilder();
            Simulator.interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            CreateCellPopulation(folder);
        }
        
        public void CreatePopulation(int popSize, int popMaxSize)
        {
            //interCellEdges = new StringBuilder();
            //interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            Simulator.interCellEdges = new StringBuilder();
            Simulator.interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            populationSize = popSize;
            maxPopulationSize = popMaxSize;
            cells = new MGCell[popMaxSize];
            sigma = new int[maxPopulationSize];
            dim = new Vector(Math.Pow(maxPopulationSize, 1f / 3f), Math.Pow(maxPopulationSize, 1f / 3f), Math.Pow(maxPopulationSize, 1f / 3f));
            halfDim = dim / 2;
            
            for (int i = 0; i < populationSize; i++)
            {
                cells[i] = new MGCell(i);
            }
        }
        
        void CreateCellPopulation(string folder)
        {
            CellProvider provider = new CellProvider();
            provider.ReadParameters(folder);
            provider.ReadResults();
            provider.PlayResults(0);

            populationSize = provider.popSize;
            maxPopulationSize = provider.maxPopSize;

            cells = new MGCell[provider.maxPopSize];
            sigma = new int[provider.maxPopSize];

            int start = 0;
            
            tissues = new List<Tissue>(provider.nbOfCellTypes);
            /*
            for (int i = 0; i < provider.nbOfCellTypes; i++)
            {
                tissues.Add(new Tissue(provider.nbOfCellsPerType[i], provider.maxPopSize));
                tissues[i].indices = new int[tissues[i].maxPopulationSize];
                
                start += i == 0 ? 0 : tissues[i - 1].populationSize;
                for (int j = 0; j < tissues[i].populationSize; j++)
                {
                    int index = start + j;

                    //Console.WriteLine(index);
                    cells[index] = new MGCell(index, provider.meshPerCell[index]);
                    //Console.WriteLine(index + ", " + provider.meshPerCell[index].vertexCount());
                    cells[index].tissueId = i;

                    tissues[i].indices[j] = index;
                    tissues[i].cells[j] = cells[index];
                    sigma[index] = index;
                    cells[index].cycleTime = new Random().Next(MGModel.cellCyclePeriod);
                    tissues[i].reference += cells[index].ComputeCentreFromMesh();
                }
                tissues[i].mesh = provider.meshPerCell[start];
                tissues[i].reference /= tissues[i].populationSize; 
            }
            */

            for (int i = 0; i < provider.nbOfCellTypesPerFrame[provider.nbOfFrames-1]; i++)
            {
                //Console.WriteLine(provider.nbOfCellsPerTypePerFrame[provider.nbOfFrames - 1][i]);
                tissues.Add(new Tissue(provider.nbOfCellsPerTypePerFrame[provider.nbOfFrames - 1][i], provider.maxPopSize));
                tissues[i].indices = new int[tissues[i].maxPopulationSize];

                start += i == 0 ? 0 : tissues[i - 1].populationSize;
                for (int j = 0; j < tissues[i].populationSize; j++)
                {
                    int index = start + j;

                    //Console.WriteLine(index);
                    cells[index] = new MGCell(index, provider.meshPerCell[index]);
                    //Console.WriteLine(index + ", " + provider.meshPerCell[index].vertexCount());
                    cells[index].tissueId = i;

                    tissues[i].indices[j] = index;
                    tissues[i].cells[j] = cells[index];
                    sigma[index] = index;
                    cells[index].cycleTime = new Random().Next(MGModel.cellCyclePeriod);
                    tissues[i].reference += cells[index].ComputeCentreFromMesh();
                }
                tissues[i].mesh = provider.meshPerCell[start];
                tissues[i].reference /= tissues[i].populationSize;
            }

            //Console.WriteLine(tissues.Capacity);
            /*
            int cellTypeStartStates = 0;
            for (int i = 0; i < provider.nbOfCellTypes; i++)
            {
                if (i > 0) cellTypeStartStates += provider.nbOfCellsPerType[i - 1];
                for (int j = 0; j < provider.nbOfCellsPerType[i]; j++)
                {
                    cells[cellTypeStartStates + j] = new MGCell(cellTypeStartStates + j, provider.meshPerCell[cellTypeStartStates + j]);
                }
            }
            //*/
        }

        public void AddCell(int tissueId, MGCell cell)
        {
            if (populationSize < maxPopulationSize)
            {
                cells[populationSize] = cell;
                tissues[tissueId].AddCell(cell);
                populationSize++;
            }
        }

        public void AddTissue(Tissue tissue)
        {
            int i = tissues.Count;
            tissue.indices = new int[tissue.maxPopulationSize];
            
            int start = populationSize;
            for (int j = 0; j < tissue.populationSize; j++)
            {
                int index = start + j;
                cells[index] = new MGCell(index, tissue.mesh);
                
                tissue.indices[j] = index;
                tissue.cells[j] = cells[index];
                sigma[index] = index;
            }
            tissues.Add(tissue);
            populationSize += tissue.populationSize;
        }

        public void PositionCells(bool random)
        {
            Vector centre = new Vector(0,0,0);

            if (!random)
            {
                for (int i = 0; i < populationSize; ++i)
                {
                    double x = Helper.NCubedInN[i].x * MGModel.Rcell * cells[i].innerRadius.x * 2;
                    double y = Helper.NCubedInN[i].y * MGModel.Rcell * cells[i].innerRadius.y * 2;
                    double z = Helper.NCubedInN[i].z * MGModel.Rcell * cells[i].innerRadius.z * 2;

                    centre.x = x - halfDim.x;
                    centre.y = y - halfDim.y;
                    centre.z = z - halfDim.z;

                    cells[i].PositionCell(centre);

                    reference.y = y;
                }
            }
            else
            {
                for (int i = 0; i < populationSize; ++i)
                {
                    Random rnd = new Random();
                    double x = rnd.Next(0, (int)dim.x) * MGModel.Rcell * cells[i].innerRadius.x * 2;
                    double y = rnd.Next(0, (int)dim.y) * MGModel.Rcell * cells[i].innerRadius.y * 2;
                    double z = rnd.Next(0, (int)dim.z) * MGModel.Rcell * cells[i].innerRadius.z * 2;
                    
                    centre.x = x - halfDim.x;
                    centre.y = y - halfDim.y;
                    centre.z = z - halfDim.z;

                    cells[i].PositionCell(centre);
                }
            }
        }

        public void PositionCells(bool random, Vector refer)
        {
            Vector centre = new Vector(0, 0, 0);
            reference = refer;

            if (!random)
            {
                double x=0, y=0, z = 0;
                for (int i = 0; i < populationSize; ++i)
                {
                    x = Helper.NCubedInN[i].x * MGModel.Rcell * cells[i].innerRadius.x * 2;
                    y = Helper.NCubedInN[i].y * MGModel.Rcell * cells[i].innerRadius.y * 2;
                    z = Helper.NCubedInN[i].z * MGModel.Rcell * cells[i].innerRadius.z * 2;

                    centre.x = x - halfDim.x;
                    centre.y = y - halfDim.y;
                    centre.z = z - halfDim.z;

                    cells[i].PositionCell(centre+reference);
                }
                reference.y += y;
            }
            else
            {
                for (int i = 0; i < populationSize; ++i)
                {
                    Random rnd = new Random();
                    double x = rnd.Next(0, (int)dim.x) * MGModel.Rcell * cells[i].innerRadius.x * 2;
                    double y = rnd.Next(0, (int)dim.y) * MGModel.Rcell * cells[i].innerRadius.y * 2;
                    double z = rnd.Next(0, (int)dim.z) * MGModel.Rcell * cells[i].innerRadius.z * 2;

                    centre.x = x - halfDim.x;
                    centre.y = y - halfDim.y;
                    centre.z = z - halfDim.z;

                    cells[i].PositionCell(reference + centre);
                }
            }
        }
        
        public void PositionVEEpiTECells(int maxPopSize, int popSize)
        {
            Vector centre = new Vector(0, 0, 0);

            int[] boundary = Helper.OrderEggCylinderCells(maxPopSize, popSize);

            for (int i = 0; i < populationSize; ++i)
            {
                double x = Helper.NCubedInN[i].x * MGModel.Rcell * cells[i].innerRadius.x * 2;
                double y = Helper.NCubedInN[i].y * MGModel.Rcell * cells[i].innerRadius.y * 2;
                double z = Helper.NCubedInN[i].z * MGModel.Rcell * cells[i].innerRadius.z * 2;

                centre.x = x - halfDim.x;
                centre.y = y - halfDim.y;
                centre.z = z - halfDim.z;

                cells[i].PositionCell(reference + centre);
            }   
        }
        
        public void SetCellNeighbourhood(int cellId)
        {
            MGCell cell = cells[cellId];
            cell.neighbours.Clear();

            if (MGModel.hexagonalNeighbourhoodForCells)
            {
                cell.neighbours = Helper.GetHexagonalNeighbourhood3D(cellId, maxPopulationSize, populationSize);
                //Console.WriteLine("Hexagonal");
            }
            else if (MGModel.mooreNeighbourhoodForCells)
            {
                //cell.neighbours = Helper.GetMooreNeighbourhood3D(cellId, maxPopulationSize, populationSize);
                cell.neighbours = Helper.GetMooreNeighbourhood3D(cellId, dim, populationSize);
                //Console.WriteLine(maxPopulationSize + ", " + populationSize);
            }
            else
            {
                //Console.WriteLine("Population Size: " + populationSize);
                //Console.WriteLine("Cell ID: " + cellId);
                for (int i = 0; i < populationSize; i++)
                {
                    double distance = Vector.Distance(cell.GetPosition(), cells[i].GetPosition());

                    //Console.WriteLine(cellId + ":" + cell.GetPosition() + ", " + i + ":" + cells[i].GetPosition() + ", " + distance + ", " + MGModel.maximumNeighbourDistance);

                    if (i != cellId && distance < MGModel.maximumNeighbourDistance && cell.neighbours.Count < 26)
                    {
                        cell.neighbours.Add(i);   
                    }
                }
            }
            //Console.WriteLine(cells[cellId].neighbours.Count);
            //if(Simulator.frame== 10)
            //{
                /*
                string str = cellId + ": ";
                for (int j = 0; j < cell.neighbours.Count; j++)
                {
                    str += cell.neighbours[j] + ",";
                }
                Console.WriteLine(str);
                //*/
            //}

        }
        
        public void DynamiseSystem()
        {
            for (int i = 0; i < populationSize; i++)
            {
                cells[i].ResetCell();
                

                if (MGModel.staticNeighbourhood)
                {
                    //if (Simulator.frame < 1 || MGModel.searchNeighbours)
                    if (Simulator.frame < 1 || Simulator.frame == MGModel.nextNeighboursSearchFrame)
                    {
                        SetCellNeighbourhood(i);
                        SetInterCellEdges(i);
                        MGModel.searchNeighbours = false;
                    }

                    cells[i].ComputeForces();
                    cells[i].ComputeExternalForces();
                    cells[i].Dynamise();
                }
                else
                {
                    if (Simulator.frame < 1 || Simulator.frame == MGModel.nextNeighboursSearchFrame)
                    //if (Simulator.frame % MGModel.neighbourPeriod == 0)
                    //if (MGModel.searchNeighbours)
                    {
                        //Console.WriteLine("Searching Neighbours ...");
                        SetCellNeighbourhood(i);
                        SetInterCellEdges(i);
                        MGModel.searchNeighbours = false;
                    }
                    
                    cells[i].ComputeForces();
                    cells[i].ComputeExternalForces();
                    cells[i].Dynamise();
                }    
            }
            frame++;

            if (Simulator.frame % Simulator.logFrequency == 0)
            {
                //Console.WriteLine("Hey!");
                //statesFile = logDir + "/" + statesFile;
                
                /*
                string particleNeighboursFile = Simulator.logDir + "/InterCellEdges.mg";
                FileStream fs = new FileStream(particleNeighboursFile, FileMode.OpenOrCreate);
                fs.Write(Encoding.ASCII.GetBytes(interCellEdges.ToString()), 0, interCellEdges.ToString().Length);
                fs.Close();
                //*/
            }
        }

        public void SetInterCellEdges(int cellId)
        {
            Edge edge;
            double distance = 0;
            cells[cellId].externalEdges.clear();

            // Permutation of particles
            sigma = new int[cells[cellId].nbOfParticles];
            
            if (Simulator.frame % 10 == 0)
            {
                for (int i=0;  i < cells[cellId].nbOfParticles; i++)
                {
                    sigma[i] = i;
                }
                Randomize();
            }


            for (int i = 0; i < cells[cellId].nbOfParticles; i++)
            {
                // To make sure all external neighbours have their external forces defined
                cells[cellId].vertices[i].externalForces = new Vector();

                int countN = 0;
                for (int j = 0; j < cells[cellId].neighbours.Count; j++)
                {
                    int l = cells[cellId].neighbours[j];
                    for (int k = 0; k < cells[l].nbOfParticles; k++)
                    {
                        //distance = Vector.Distance(cells[cellId].vertices[i].GetPosition(), cells[l].vertices[k].GetPosition());
                        //if (distance < MGModel.DInteraction[cells[cellId].tissueId, cells[l].tissueId] && countN < MGModel.maxExternalNeighbours)
                        distance = Vector.Distance(cells[cellId].vertices[i].GetPosition(), cells[l].vertices[sigma[k]].GetPosition());
                        if (distance < MGModel.DInteraction[cells[cellId].tissueId, cells[l].tissueId] && countN < MGModel.maxExternalNeighbours)
                        {
                            countN++;

                            edge = new Edge(cells[cellId].vertices[i], cells[l].vertices[sigma[k]]);
                            cells[cellId].externalEdges.add(edge);
                            
                        }
                    }
                }
            }
        }

        public void SetInterCellEdgesIfNotExisting(int cellId)
        {
            Edge edge;
            double distance = 0;
            //cells[cellId].externalEdges.clear();


            for (int i = 0; i < cells[cellId].nbOfParticles; i++)
            {
                // To make sure all external neighbours have their external forces defined
                cells[cellId].vertices[i].externalForces = new Vector();

                if (cells[cellId].vertices[i].externalNeighbours.Count == 0)
                {
                    int countN = 0;
                    for (int j = 0; j < cells[cellId].neighbours.Count; j++)
                    {
                        int l = cells[cellId].neighbours[j];
                        for (int k = 0; k < cells[l].nbOfParticles; k++)
                        {
                            distance = Vector.Distance(cells[cellId].vertices[i].GetPosition(), cells[l].vertices[k].GetPosition());
                            if (distance < MGModel.DInteraction[cells[cellId].tissueId, cells[l].tissueId] && countN < MGModel.maxExternalNeighbours)
                            {
                                countN++;

                                edge = new Edge(cells[cellId].vertices[i], cells[l].vertices[k]);
                                cells[cellId].externalEdges.add(edge);
                                cells[cellId].vertices[i].externalNeighbours.Add(new int[] { l, k });
                                //Console.WriteLine("Woo hoo");

                                //interCellEdges.AppendLine( frame + "," + cellId + "," + i + "," + l + "," + k);
                            }
                        }
                    }
                }
            }
        }

        public void StochasticDynamiseSystem()
        {
            Parallel.For(0, populationSize, i =>
            {
                cells[i].ResetCell();
                SetCellNeighbourhood(i);
                SetInterCellEdges(i);
                cells[i].ChangeShape();
                cells[i].ComputeForces();
                cells[i].Dynamise();
            });
        }
        
        public void Cleavage(bool staticShape)
        {
            Vector normal = populationSize == 1 ? new Vector(1, 0, 0) : Vector.zero;
            for (int i = 0; i < populationSize; i++)
            {
                cells[i].CellCycle(appliedForces, normal, staticShape);
            }
        }

        public void Proliferation(bool staticShape)
        {
            //Console.WriteLine("Proliferate...");
            Vector normal = populationSize == 1 ? new Vector(0, 1, 0) : Vector.zero;

            for (int i = 0; i < populationSize; i++)
            {
                if(populationSize < maxPopulationSize)
                {
                    cells[i].KleinCellCycle(appliedForces, normal, staticShape);       
                }
            }
        }

        public void EPIProliferation(bool staticShape)
        {
            //Console.WriteLine("Proliferate...");
            Vector normal = new Vector(0,0,0);

            if(populationSize == 1)
            {
                normal = new Vector(0, 1, 0);
            }else if(populationSize == 2)
            {
                normal = new Vector(0, 0, 1);
            }
            else
            {
                normal = new Vector(1, 0, 0);
            }
            for (int i = 0; i < populationSize; i++)
            {
                if (populationSize < maxPopulationSize)
                {
                    cells[i].KleinCellCycle(appliedForces, normal, staticShape);
                }
            }
        }

        public void PlanarProliferation()
        {
            double x = new Random().NextDouble();
            double z = new Random().NextDouble();
            Vector normal = new Vector(x, 0, z);
            for (int i = 0; i < populationSize; i++)
            {
                cells[i].CellCycle(appliedForces, normal, false);
            }
        }

        public void Randomize()
        {
            Randomizer.Randomize<int>(sigma);
        }

        public Bounds ComputeBounds()
        {
            Vector min = new Vector();
            Vector max = new Vector();
            Vector centre = new Vector();
            for (int i = 0; i < populationSize; i++)
            {
                for (int j = 0; j < cells[i].nbOfParticles; j++)
                {
                    if (cells[i].vertices[j].v.x < min.x)
                    {
                        min.x = cells[i].vertices[j].v.x;
                    }
                    if (cells[i].vertices[j].v.y < min.y)
                    {
                        min.y = cells[i].vertices[j].v.y;
                    }
                    if (cells[i].vertices[j].v.z < min.z)
                    {
                        min.z = cells[i].vertices[j].v.z;
                    }
                    if (cells[i].vertices[j].v.x > max.x)
                    {
                        max.x = cells[i].vertices[j].v.x;
                    }
                    if (cells[i].vertices[j].v.y > max.y)
                    {
                        max.y = cells[i].vertices[j].v.y;
                    }
                    if (cells[i].vertices[j].v.z > max.z)
                    {
                        max.z = cells[i].vertices[j].v.z;
                    }
                }
                centre += cells[i].ComputeCentreFromMesh();
            }
            centre /= populationSize;

            bounds = new Bounds();
            bounds.centre = centre;
            bounds.size = max - min;
            bounds.extents = bounds.size / 2;
            bounds.max = max;
            bounds.min = min;

            Mesh[] meshes = new Mesh[populationSize];
            for (int i = 0; i < populationSize; i++)
            {
                meshes[i] = new Mesh();
                meshes[i].copy(cells[i]);
            }
            bounds.convexHull = MGConvexHull.GenerateConvexHull(meshes);


            Bounds bds = new Bounds(bounds);
            return bds;
        }

        public Bounds ComputeBounds(bool convexHull)
        {
            Vector min = new Vector();
            Vector max = new Vector();
            Vector centre = new Vector();
            for (int i = 0; i < populationSize; i++)
            {
                for (int j = 0; j < cells[i].nbOfParticles; j++)
                {
                    if (cells[i].vertices[j].v.x < min.x)
                    {
                        min.x = cells[i].vertices[j].v.x;
                    }
                    if (cells[i].vertices[j].v.y < min.y)
                    {
                        min.y = cells[i].vertices[j].v.y;
                    }
                    if (cells[i].vertices[j].v.z < min.z)
                    {
                        min.z = cells[i].vertices[j].v.z;
                    }
                    if (cells[i].vertices[j].v.x > max.x)
                    {
                        max.x = cells[i].vertices[j].v.x;
                    }
                    if (cells[i].vertices[j].v.y > max.y)
                    {
                        max.y = cells[i].vertices[j].v.y;
                    }
                    if (cells[i].vertices[j].v.z > max.z)
                    {
                        max.z = cells[i].vertices[j].v.z;
                    }
                }
                centre += cells[i].ComputeCentreFromMesh();
            }
            centre /= populationSize;

            bounds = new Bounds();
            bounds.centre = centre;
            bounds.size = max - min;
            bounds.extents = bounds.size / 2;
            bounds.max = max;
            bounds.min = min;

            Mesh[] meshes = new Mesh[populationSize];
            for (int i = 0; i < populationSize; i++)
            {
                meshes[i] = new Mesh();
                meshes[i].copy(cells[i]);
            }

            if (convexHull)
                bounds.convexHull = MGConvexHull.GenerateConvexHull(meshes);
            else
                bounds.convexHull = cells[0];


            Bounds bds = new Bounds(bounds);
            return bds;
        }

        public void Translate(Vector u)
        {
            for(int i = 0; i < populationSize; i++)
            {
                cells[i].Translate(u);
            }

            bounds.centre += u;
        }

        public double Kronecker(Vector u, Bounds bds)
        {
            double val = 0;
            if (bounds.convexHull.isPointInConvexMesh(u))
            {
                val = 1f;
            }
            return val;
        }

        public double Kronecker1(Vector u, Bounds bds)
        {
            float val = 0f;
            Vector v = new Vector((u.x - bounds.centre.x) / bounds.extents.x,
                                    (u.y - bounds.centre.y) / bounds.extents.y,
                                    (u.z - bounds.centre.z) / bounds.extents.z);

            if (v.sqrNorm() <= 1)
            {
                val = 1f;
            }
            else
            {
                val = 0f;
            }
            return val;
        }

        public void FromXZPlanarToEllipsoid(Bounds ellipsoidBounds)
        {
            Vector temp = new Vector();
            Vector radius = ellipsoidBounds.extents;
            double d = 0, R=0;

            for(int i=0; i<populationSize; i++)
            {
                for(int j=0; j<cells[i].nbOfParticles; j++)
                {
                    d = Math.Abs(ellipsoidBounds.centre.y - cells[i].vertices[j].v.y) - ellipsoidBounds.extents.y;
                    radius = ellipsoidBounds.extents + new Vector(d,d,d);

                    temp = (cells[i].vertices[j].v - ellipsoidBounds.centre);
                    R = temp.norm();

                    cells[i].vertices[j].v = new Vector(
                                                            temp.x*radius.x/R + ellipsoidBounds.centre.x,
                                                            temp.y * radius.y / R + ellipsoidBounds.centre.y,
                                                            temp.z * radius.z / R + ellipsoidBounds.centre.z
                                                       );
                }

                d = Math.Abs(ellipsoidBounds.centre.y - cells[i].nuclei[0].v.y) - ellipsoidBounds.extents.y;
                radius = ellipsoidBounds.extents + new Vector(d, d, d);

                temp = (cells[i].nuclei[0].v - ellipsoidBounds.centre);
                R = temp.norm();

                //*
                cells[i].nuclei[0].v = new Vector(
                                                        temp.x * radius.x / R + ellipsoidBounds.centre.x,
                                                        temp.y * radius.y / R + ellipsoidBounds.centre.y,
                                                        temp.z * radius.z / R + ellipsoidBounds.centre.z
                                                   );
                //*/

                //Console.WriteLine(i + ": " + cells[i].nuclei[0].v);
            }
        }

        public void FromXZSquareToXZCircle()
        {
            //double radius = Math.Sqrt(bounds.extents.x * bounds.extents.x + bounds.extents.z * bounds.extents.z);
            double radius = bounds.extents.x;

            for (int i=0; i<populationSize; i++)
            {
                for(int j=0; j<cells[i].nbOfParticles; j++)
                {
                    double y = cells[i].vertices[j].v.y;
                    Vector CP = cells[i].vertices[j].v-bounds.centre;
                    CP.y = 0;

                    if(CP.norm() > radius)
                    {
                        cells[i].vertices[j].v = radius * CP/CP.norm() + bounds.centre;
                        cells[i].vertices[j].v.y = y;
                    }
                }
            }
        }

        public void FindUpperParticlesInXZPlanarTissue()
        {    
            for (int i = 0; i < populationSize; i++)
            {
                cells[i].subset = new List<int>();
                for (int j = 0; j < cells[i].nbOfParticles; j++)
                {
                    Vector centre = cells[i].ComputeCentreFromMesh();
                    if(cells[i].vertices[j].v.y > centre.y)
                    {
                        cells[i].subset.Add(j);
                    }
                }
            }
        }

        public void FindLowerParticlesInXZPlanarTissue()
        {
            subset = new List<Vector>();
            for (int i = 0; i < populationSize; i++)
            {
                cells[i].subset = new List<int>();
                for (int j = 0; j < cells[i].nbOfParticles; j++)
                {
                    //cells[i].vertices[j].normal = new Vector(0,0,0);
                    Vector centre = cells[i].ComputeCentreFromMesh();
                    if (cells[i].vertices[j].v.y < centre.y - 0.8 * cells[i].innerRadius.y)
                    {
                        cells[i].subset.Add(j);
                        subset.Add(cells[i].vertices[j].v);
                    }
                }
                //Console.WriteLine(i + "," + cells[i].subset.Count);
            }
        }

        public void WrapTissue(Bounds bounds, Vector fictiveCentre, Vector normal, double threshold)
        {
            float scalingFactor = 10000;
            float normalScale = 0.2f;

            for(int i=0; i<populationSize; i++)
            {
                cells[i].ComputeVertexNormalsOnSubset();
                //*
                for(int j=0; j<cells[i].nbOfParticles; j++)
                {
                    cells[i].vertices[j].globalForces = new Vector();
                }
                //*/

                for (int j=0; j<cells[i].subset.Count; j++)
                {
                    Vector IJK = bounds.DichotomicSearch(cells[i].vertices[cells[i].subset[j]].v);

                    if (IJK.x >= 0 && IJK.x < bounds.resolution && IJK.y >= 0 && IJK.y < bounds.resolution && IJK.z >= 0 && IJK.z < bounds.resolution)
                    {
                        //cells[i].vertices[cells[i].subset[j]].globalForces = -bounds.gradPeronaMalik[(int)IJK.x, (int)IJK.y, (int)IJK.z];

                        normal = cells[i].vertices[cells[i].subset[j]].normal;

                        //*
                        if (Simulator.frame < 2)
                        {
                            //Vector normal = .5*(bounds.centre + fictiveCentre) - cells[i].vertices[cells[i].subset[j]].v;
                            //normal = bounds.centre - cells[i].vertices[cells[i].subset[j]].v;
                            //Vector normal = bounds.centre - cells[i].GetPosition();
                            //Console.WriteLine(i + ", " + cells[i].subset[j] + ", " + bounds.vPHI[(int)IJK.x, (int)IJK.y, (int)IJK.z]);
                            //normal = cells[i].vertices[cells[i].subset[j]].normal;
                            //normal = new Vector(0, -1, 0);
                            //normalScale = .2f;
                            //Console.WriteLine(cells[i].subset[j] + ", " + normal);
                        }

                        normal.normalize();
                        Vector globalForces = new Vector(0, 0, 0);
                        cells[i].vertices[cells[i].subset[j]].nullForces = true;
                        if (bounds.vPHI[(int)IJK.x, (int)IJK.y, (int)IJK.z] < threshold)
                        {
                            //if (Simulator.frame < 10) Console.WriteLine("Wrapping woo hooo !");

                            double dot = -bounds.gradPeronaMalik[(int)IJK.x, (int)IJK.y, (int)IJK.z] * normal;
                            globalForces = normal //* 
                                                * (normalScale * Math.Sign(dot)
                                                    + dot
                                                    );
                            //*/
                            cells[i].vertices[cells[i].subset[j]].nullForces = false;
                        }
                        cells[i].vertices[cells[i].subset[j]].globalForces = scalingFactor * globalForces;
                        //*/
                    }
                }
            }
        }

        public Vector ComputeTissueNormalsOnXZPlanarTissue(int i)
        {
            Vector u, v;
            int dimX = (int)dim.x;
            int dimY = (int)dim.y;
            int dimZ = (int)dim.z;

            if (i % dimX == 0)
            {
                u = 2 * (cells[i + 1].ComputeCentreFromMesh() - cells[i].ComputeCentreFromMesh());
            }
            else if (i % dimX == dimZ - 1)
            {
                u = 2 * (cells[i].ComputeCentreFromMesh() - cells[i - 1].ComputeCentreFromMesh());
            }
            else
            {
                u = cells[i + 1].ComputeCentreFromMesh() - cells[i - 1].ComputeCentreFromMesh();
            }

            if (i / dimZ == 0)
            {
                v = 2 * (cells[i + dimZ].ComputeCentreFromMesh() - cells[i].ComputeCentreFromMesh());
            }
            else if (i / dimZ == dimX - 1)
            {
                v = 2 * (cells[i].ComputeCentreFromMesh() - cells[i - dimZ].ComputeCentreFromMesh());
            }
            else
            {
                v = cells[i + dimZ].ComputeCentreFromMesh() - cells[i - dimZ].ComputeCentreFromMesh();
            }

            return u ^ v;
        }
        
        public void LogNumbers()
        {
            PersistantNumbers pv = 
                new PersistantNumbers(Simulator.frame, tissues.Count, populationSize);
            Simulator.numbers.Add(pv);
            for(int i = 0; i < tissues.Count; i++)
            {
                pv = new PersistantNumbers(Simulator.frame, tissues[i].id, tissues[i].populationSize);
                Simulator.numbers.Add(pv);
                for (int j=0; j < tissues[i].populationSize; j++)
                {
                    pv = new PersistantNumbers(Simulator.frame, tissues[i].cells[j].cellId, tissues[i].cells[j].nbOfParticles, i);
                    Simulator.numbers.Add(pv);
                }
            }
        }

        public void LogFinalNumbers()
        {
            PersistantNumbers pv =
                new PersistantNumbers(Simulator.frame, tissues.Count, populationSize);
            Simulator.finalNumbers.Add(pv);
            for (int i = 0; i < tissues.Count; i++)
            {
                pv = new PersistantNumbers(Simulator.frame, tissues[i].id, tissues[i].populationSize);
                Simulator.finalNumbers.Add(pv);
                for (int j = 0; j < tissues[i].populationSize; j++)
                {
                    pv = new PersistantNumbers(Simulator.frame, tissues[i].cells[j].cellId, tissues[i].cells[j].nbOfParticles, i);
                    Simulator.finalNumbers.Add(pv);
                }
            }
        }
        
        public void LogMetrics()
        {
            PersistantMetrics pm = new PersistantMetrics(Simulator.frame, Metrics.ElasticEnergy());
            Simulator.metrics.Add(pm);
        }

        public void LogNeighbourhoods()
        {
            for(int i=0; i<populationSize; i++)
            {
                for(int j=0; j<cells[i].externalEdges.getCount(); j++)
                {
                    Simulator.interCellEdges.AppendLine(Simulator.frame + ";" + 
                                                cells[i].externalEdges[j].ends[0].cellId + ";" +
                                                cells[i].externalEdges[j].ends[0].pos + ";" +
                                                cells[i].externalEdges[j].ends[1].cellId + ";" +
                                                cells[i].externalEdges[j].ends[1].pos);
                }
            }
        }

        /*
        public void LogState()
        {
            PersistantVertex pv;
            for (int i = 0; i < populationSize; i++)
            {
                for (int j = 0; j < cells[i].nbOfParticles; j++)
                {
                    pv = new PersistantVertex(Simulator.frame, cells[i].vertices[j].id, cells[i].cellId, cells[i].vertices[j].v);
                    Simulator.states.Add(pv);
                }
                pv = new PersistantVertex(Simulator.frame, cells[i].nuclei[0].id, cells[i].cellId, cells[i].nuclei[0].v);
                Simulator.states.Add(pv);
            }
        }

        public void LogFinalState()
        {
            PersistantVertex pv;
            for (int i = 0; i < populationSize; i++)
            {
                for (int j = 0; j < cells[i].nbOfParticles; j++)
                {
                    pv = new PersistantVertex(Simulator.frame, cells[i].vertices[j].id, cells[i].cellId, cells[i].vertices[j].v);
                    Simulator.finalStates.Add(pv);
                }
                pv = new PersistantVertex(Simulator.frame, cells[i].nuclei[0].id, cells[i].cellId, cells[i].nuclei[0].v);
                Simulator.finalStates.Add(pv);
            }
        }
        //*/

        //*
        public void LogState()
        {
            PersistantVertex pv;
            for (int i = 0; i < tissues.Count; i++)
            {
                for(int j=0; j<tissues[i].populationSize; j++)
                {
                    for(int k=0; k<tissues[i].cells[j].nbOfParticles; k++)
                    {
                        pv = new PersistantVertex(Simulator.frame, tissues[i].cells[j].vertices[k].id, tissues[i].cells[j].cellId, tissues[i].cells[j].vertices[k].v);
                        Simulator.states.Add(pv.Clone());
                    }
                    pv = new PersistantVertex(Simulator.frame, tissues[i].cells[j].nuclei[0].id, tissues[i].cells[j].cellId, tissues[i].cells[j].nuclei[0].v);
                    Simulator.states.Add(pv.Clone());
                }
            }
        }

        public void LogFinalState()
        {
            PersistantVertex pv;
            for (int i = 0; i < tissues.Count; i++)
            {
                for (int j = 0; j < tissues[i].populationSize; j++)
                {
                    for (int k = 0; k < tissues[i].cells[j].nbOfParticles; k++)
                    {
                        pv = new PersistantVertex(Simulator.frame, tissues[i].cells[j].vertices[k].id, tissues[i].cells[j].cellId, tissues[i].cells[j].vertices[k].v);
                        Simulator.finalStates.Add(pv.Clone());
                    }
                    pv = new PersistantVertex(Simulator.frame, tissues[i].cells[j].nuclei[0].id, tissues[i].cells[j].cellId, tissues[i].cells[j].nuclei[0].v);
                    Simulator.finalStates.Add(pv.Clone());
                }
            }
        }
        //*/
    }
}