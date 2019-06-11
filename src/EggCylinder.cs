using System;
using System.Collections.Generic;
using MGSharp.Core.GeometricPrimitives;
using MGSharp.Core.MGCellPopulation;
using MGSharp.Core.MGModels;
using MGSharp.Core.Helpers;
using System.Diagnostics;

namespace MGSharp
{
    class EggCylinder: Simulator
    {
        public static void Main(string[] args)
        {
            Simulator simulator = new EggCylinder();
            Simulator.states = new List<PersistantVertex>();
            Simulator.numbers = new List<PersistantNumbers>();
            Simulator.finalStates = new List<PersistantVertex>();
            Simulator.finalNumbers = new List<PersistantNumbers>();
            Simulator.metrics = new List<PersistantMetrics>();

            MGModel.Rcell = 1;
            simulator.SetupSimulation();
            simulator.SetModel();
            simulator.SetInitialConditions();

            simulator.LogParameters();

            Stopwatch stopwatch = new Stopwatch();
            Console.WriteLine("Simulation in Progress ...");
            stopwatch.Start();
            simulator.Update();
            stopwatch.Stop();
            Console.WriteLine("Simulation ended succesfully in " + stopwatch.ElapsedMilliseconds + " milliseconds.");
            simulator.Log();

            Console.WriteLine("Simulation Logs have been written to the Log file.");
            Console.ReadKey();
        }
        
        public override void SetupSimulation()
        {
            logFrequency = 1;
            nbCellTypes = 3;
            nbOfSimulationSteps = 10;
            
            //*
            popSize = 148;
            popMaxSize = 343;

            Tissue t1 = new Tissue(25, 125, new Cylinder34(new Vector(.5f, .5f, .5f)));
            Tissue t2 = new Tissue(98, popMaxSize, new Cylinder34(new Vector(.5f, 1f, .5f)));
            Tissue t3 = new Tissue(25, 125, new Cylinder34(new Vector(.5f, .5f, .5f)));
            List<Tissue> tissueList = new List<Tissue>() { t1, t2, t3 };

            cellPopulation = new CellPopulation(popSize, popMaxSize, tissueList, true);
            //Helper.SetNCubedInN(125);
            //t1.PositionCells(false);
            //Helper.SetNCubedInN(popMaxSize);
            //t2.PositionCells(false, new Vector(0, t1.reference.y + t1.mesh.innerRadius.y + t2.mesh.innerRadius.y, 0));

            Helper.SetNCubedInN(popMaxSize);
            t2.PositionCells(false);
            Helper.SetNCubedInN(125);
            t1.PositionCells(false, new Vector(0, t2.reference.y - 4*t2.mesh.innerRadius.y - t1.mesh.innerRadius.y, 0));
            Helper.SetNCubedInN(125);
            t3.PositionCells(false, new Vector(0, t2.reference.y + t3.mesh.innerRadius.y, 0));
            //*/

            nbCellTypes = tissueList.Count;
        }

        int[] boundary25;
        public override void SetInitialConditions()
        {
            name = "Egg Cylinder";
            int[] foyer = { 8,9,10,11,12,
                            15,16,17,18,19,
                            22,23,24,25,26,
                            29,30,31,32,33,
                            36,37,38,39,40
            };

            boundary25 = new int[]{ 0,1,2,3,
                                    4,9,14,19,
                                    24,23,22,21,
                                    20,15,10,5
            };

            
            //*
            for (int i = 0; i < foyer.Length; i++)
            {
                cellPopulation.cells[foyer[i] + 25].ApicalConstriction(.5f);
                cellPopulation.cells[foyer[i] + 25].ApicalConstriction(.5f);
            }
            for (int i = 0; i < foyer.Length; i++)
            {
                cellPopulation.cells[foyer[i] + 74].polarisation = Vector.down;
                cellPopulation.cells[foyer[i] + 74].ApicalConstriction(.5f);
            }
            //*/
        }

        public override void Update()
        {
            while (frame < nbOfSimulationSteps)
            {
                
                if (frame % logFrequency == 0)
                {
                    cellPopulation.LogState();
                    cellPopulation.LogNumbers();
                    cellPopulation.LogMetrics();
                }
                cellPopulation.DynamiseSystem();
                Console.WriteLine(frame);
                frame++;
            }
            cellPopulation.LogState();
            cellPopulation.LogNumbers();
            cellPopulation.LogMetrics();
            cellPopulation.LogFinalState();
            cellPopulation.LogFinalNumbers();
        }
        
        
        public void TEFold()
        {
            //MGModel.J[1, 1] = 1f;
            int[] foyer = new int[] { 56, 57, 58,
                                      61, 62, 63,
                                      66, 67, 68 };

            for (int i = 0; i < foyer.Length; i++)
            {
                //cellPopulation.cells[foyer[i]].ApicalConstriction(1f);
                // Constrict apically
                cellPopulation.cells[foyer[i]].SetElengths(cellPopulation.cells[0]);
            }

            for (int i = 0; i < 50; i++)
            {
                cellPopulation.cells[i].SetEdgeELengths();
            }
        }


        public override void SetModel()
        {
            MGModel.dT = 0.02f;
            MGModel.Rcell = 1f;
            MGModel.maximumNeighbourDistance = 3f * MGModel.Rcell;
            MGModel.DInt = 0.02f;
            MGModel.DCol = 0.01f;
            MGModel.rho = 1.0f;
            MGModel.alpha = 2f;
            MGModel.u0 = 0.2f;
            MGModel.u1 = 0.02f;
            MGModel.springForce = 1f;
            MGModel.damping = 0.5f;
            MGModel.delta = 0.001f;
            MGModel.divisionRate = .5f;
            MGModel.cellCyclePeriod = 5000;
            MGModel.maxExternalNeighbours = 7;
            MGModel.elasticExternalSpring = false;
            MGModel.mooreNeighbourhoodForCells = false;
            MGModel.staticNeighbourhood = true;

            MGModel.J = new float[nbCellTypes + 1, nbCellTypes];
            MGModel.DInteraction = new float[nbCellTypes + 1, nbCellTypes];
            for (int i = 0; i < nbCellTypes + 1; i++)
            {
                for (int j = 0; j < nbCellTypes; j++)
                {
                    MGModel.J[i, j] = 1f;
                    MGModel.DInteraction[i, j] = 2.5f * MGModel.DInt;
                }
            }
            
            //MGModel.DInteraction[0, 0] = 2.5f * MGModel.DInt;
            //MGModel.DInteraction[1, 1] = 2.5f * MGModel.DInt;
            //MGModel.DInteraction[0, 1] = 2.5f * MGModel.DInt;
            //MGModel.DInteraction[1, 0] = 2.5f * MGModel.DInt;
        }
    }
}