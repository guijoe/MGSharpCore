using System;
using System.IO;
using System.Text;
using System.Linq;
using System.Collections.Generic;
using MGSharp.Core.MGCellPopulation;
using MGSharp.Core.MGModels;
using MGSharp.Core.Helpers;


namespace MGSharp
{
    class Simulator
    {
        public int nbCellTypes;
        public int nbOfSimulationSteps = 2000;
        public int popSize = 1;
        public int popMaxSize = 32;
        public static int logFrequency = 2;
        public static CellPopulation cellPopulation;
        public static string name;
        public static string logDir;
        public static List<PersistantVertex> states;
        public static List<PersistantNumbers> numbers;
        public static List<PersistantVertex> finalStates;
        public static List<PersistantNumbers> finalNumbers;
        public static List<PersistantMetrics> metrics;
        public static StringBuilder interCellEdges;
        public static bool logVTK = false;
        
        public static int frame = 0;
        public static float divisionRate;

        public virtual void SetupSimulation(){}

        public virtual void SetInitialConditions(){}

        public virtual void Update(){}

        public virtual void SetModel(){}

        public void Randomize()
        {
            cellPopulation.Randomize();
            for (int i = 0; i < cellPopulation.populationSize; i++)
            {
                cellPopulation.cells[i].Randomize();
            }
        }

        public void Log()
        {
            string statesFile = "Results.mg";
            statesFile = logDir + "/" + statesFile;
            FileStream fs = new FileStream(statesFile, FileMode.OpenOrCreate);
            
            CsvSerializer<PersistantVertex> serializer = new CsvSerializer<PersistantVertex>();
            serializer.Separator = ';';
            serializer.Serialize(fs, states);
            fs.Close();

            string finalStatesFile = "Final_States.mg";
            finalStatesFile = logDir + "/" + finalStatesFile;
            fs = new FileStream(finalStatesFile, FileMode.OpenOrCreate);

            serializer.Separator = ';';
            serializer.Serialize(fs, finalStates);
            fs.Close();


            string numbersFile = "Numbers.mg";
            numbersFile = logDir + "/" + numbersFile;
            fs = new FileStream(numbersFile, FileMode.OpenOrCreate);

            CsvSerializer<PersistantNumbers> serializer1 = new CsvSerializer<PersistantNumbers>();
            serializer1.Separator = ';';
            serializer1.Serialize(fs, numbers);
            fs.Close();

            string finalNumbersFile = "Final_Numbers.mg";
            finalNumbersFile = logDir + "/" + finalNumbersFile;
            fs = new FileStream(finalNumbersFile, FileMode.OpenOrCreate);
            serializer1.Separator = ';';
            serializer1.Serialize(fs, finalNumbers);
            fs.Close();

            //Log Metrics
            CsvSerializer<PersistantMetrics> serializer2 = new CsvSerializer<PersistantMetrics>();
            string metricsFile = "Metrics.mg";
            metricsFile = logDir + "/" + metricsFile;
            fs = new FileStream(metricsFile, FileMode.OpenOrCreate);
            serializer2.Separator = ';';
            serializer2.Serialize(fs, metrics);
            fs.Close();


            // Neighbourhoods
            string particleNeighboursFile = logDir + "/InterCellEdges.mg";
            fs = new FileStream(particleNeighboursFile, FileMode.OpenOrCreate);
            fs.Write(Encoding.ASCII.GetBytes(interCellEdges.ToString()), 0, interCellEdges.ToString().Length);
            fs.Close();

            // VTK
            if (logVTK)
            {
                LogVTK();
            }
        }

        public void LogVTK()
        {
            // VTK
            string vtk = "Results_vtk";
            vtk = logDir + "/" + vtk;

            int nbOfFrames = nbOfSimulationSteps / logFrequency;

            for (int fr = 0; fr < nbOfFrames; fr++)
            {
                string vtkFile = vtk + fr + ".vtk";

                FileStream fs = new FileStream(vtkFile, FileMode.OpenOrCreate);

                StreamWriter sw = new StreamWriter(fs);
                sw.AutoFlush = true;
                sw.WriteLine("# vtk DataFile Version 2.0");
                sw.WriteLine("MGSHARP Data");
                sw.WriteLine("ASCII");
                sw.WriteLine("DATASET POLYDATA");

                var frStates = from s in states
                                           where s.frame == logFrequency * fr
                                           select s;

                sw.WriteLine("POINTS " + (frStates.Count<PersistantVertex>()) + " float");
                foreach(var state in frStates)
                {
                    sw.WriteLine(state.v.ToString().Replace(',', ' '));
                }

                // POLYGONS
                StringBuilder polygons = new StringBuilder();
                
                var frNumbers = from n in numbers
                               where n.frame == logFrequency*fr
                               select n;

                int line = 0;
                int nTissues = frNumbers.ElementAt<PersistantNumbers>(line).n1;
                int pSize = frNumbers.ElementAt<PersistantNumbers>(line).n2;

                line++;
                int[] nCellsPerTissue = new int[nTissues];
                nCellsPerTissue[0] = frNumbers.ElementAt<PersistantNumbers>(line).n2;

                int faceCount = 0;
                int cellStart = 0;
                
                int nVerticesPerCell = 0;

                for (int i=0; i< nTissues; i++) {
                    for (int j=0; j<nCellsPerTissue[i]; j++)
                    {
                        line++;
                        nVerticesPerCell = frNumbers.ElementAt<PersistantNumbers>(line).n2;

                        int cellIndex = j;
                        if (i != 0)
                        {
                            cellIndex = nCellsPerTissue[i-1] + j;
                        }
                        
                        for (int k = 0; k < cellPopulation.cells[cellIndex].faceCount(); k++)
                        {
                            faceCount++;
                            polygons.AppendLine("3 " + (cellStart + cellPopulation.cells[cellIndex].faces[k].vertices[0].pos) + " "
                                        + (cellStart + cellPopulation.cells[cellIndex].faces[k].vertices[1].pos) + " "
                                        + (cellStart + cellPopulation.cells[cellIndex].faces[k].vertices[2].pos));
                        }

                        /*
                        for (int k=0; k<cellPopulation.tissues[i].mesh.faceCount(); k++)
                        {   
                            faceCount++;
                            polygons.AppendLine("3 " + (cellStart + cellPopulation.tissues[i].mesh.faces[k].vertices[0].pos) + " "
                                        + (cellStart + cellPopulation.tissues[i].mesh.faces[k].vertices[1].pos) + " "
                                        + (cellStart + cellPopulation.tissues[i].mesh.faces[k].vertices[2].pos));
                        }
                        */
                        cellStart += nVerticesPerCell + 1;
                    }

                    line++;
                    
                    if(i<nTissues-1)
                        nCellsPerTissue[i+1] = frNumbers.ElementAt<PersistantNumbers>(line).n2;
                }

                sw.WriteLine("POLYGONS " + faceCount + " " + (4*faceCount));
                sw.Write(polygons.ToString());

                fs.Close();
            }
        }

        public void LogParameters()
        {
            if (String.IsNullOrEmpty(logDir))
            {
                double timeStamp = (DateTime.Now.ToUniversalTime() - new DateTime(2000, 1, 1)).TotalMilliseconds;
                logDir = "../../Logs/" + name + "_" + timeStamp.ToString();

                Directory.CreateDirectory(logDir);
            }
            
            string parFile = "Parameters.mg";
            string logFile = logDir + "/" + parFile;
            
            StringBuilder parameters = new StringBuilder();
            parameters.AppendLine("# GENERAL SIMULATION PARAMETERS");
            parameters.AppendLine("NrDevSteps=" + nbOfSimulationSteps/logFrequency);

            MGModel model = new MGModel();
            var properties = model.GetType().GetProperties();
            for (int i = 0; i < properties.Length; i++)
            {
                parameters.AppendLine(properties[i].Name + "=" + properties[i].GetValue(model).ToString());
            }

            parameters.AppendLine();
            parameters.AppendLine("# CELL TYPES");
            parameters.AppendLine("NrCells=" + popSize);
            parameters.AppendLine("NrCellsMax=" + popMaxSize);
            parameters.AppendLine("NrCellTypes=" + nbCellTypes);

            for (int i = 0; i < nbCellTypes; i++)
            {
                parameters.AppendLine("NrCell" + i + "=" + cellPopulation.tissues[i].populationSize);
            }

            parameters.AppendLine();
            for (int i = 0; i < nbCellTypes; i++)
            {
                parameters.AppendLine("JCell" + i + "med=" + MGModel.J[0, i]);
            }

            for (int i = 0; i < nbCellTypes; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    parameters.AppendLine("JCell" + i + "Cell" + j + "=" + MGModel.J[i + 1, j]);
                }
            }

            parameters.AppendLine();
            parameters.AppendLine("# CELL MESHES (TRIANGULATION)");
            for (int i = 0; i < nbCellTypes; i++)
            {
                parameters.AppendLine();
                parameters.AppendLine("# CELL " + i);
                parameters.AppendLine("NrVertices=" + cellPopulation.tissues[i].mesh.vertexCount());
                parameters.AppendLine("NrTriangles=" + cellPopulation.tissues[i].mesh.faceCount());
                parameters.AppendLine();
                parameters.AppendLine("# VERTICES");
                parameters.AppendLine(cellPopulation.tissues[i].mesh.PrintVertices());
                parameters.AppendLine("# TRIANGLES");
                parameters.AppendLine(cellPopulation.tissues[i].mesh.PrintTriangles());
            }

            File.WriteAllText(logFile, parameters.ToString());
        }
    }
}