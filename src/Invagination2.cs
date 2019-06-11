using System;
using System.Collections.Generic;
using MGSharp.Core.GeometricPrimitives;
using MGSharp.Core.MGCellPopulation;
using MGSharp.Core.MGModels;
using MGSharp.Core.Helpers;
using System.Diagnostics;

namespace MGSharp
{
    class Invagination2: Simulator
    {
        public static void Main(string[] args)
        {
            Simulator simulator = new Invagination2();
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
            nbCellTypes = 1;

            popSize = 25;
            popMaxSize = 125;
            nbOfSimulationSteps = 1000;

            //Cylinder34 mesh = new Cylinder34();

            //HexagonalCylinder42 mesh = new HexagonalCylinder42((float)Math.Sqrt(3) / 3);

            HexagonalCylinder42 mesh = new HexagonalCylinder42(new Vector(.945,.5,.945));

            Tissue t1 = new Tissue(popSize, popMaxSize, mesh);
            List<Tissue> tissueList = new List<Tissue>() { t1 };

            nbCellTypes = tissueList.Count;
            cellPopulation = new CellPopulation(popSize, popMaxSize, tissueList);

            Helper.SetNCubedInHexa2(popMaxSize);
            cellPopulation.PositionCells(false);
        }

        public override void SetInitialConditions()
        {
            name = "Invagination";
            int[] foyer = new int[] { 20, 21, 22, 23, 24, 29, 30, 31, 32, 33, 38, 39, 40, 41, 42, 47, 48, 49, 50, 51, 56, 57, 58, 59, 60 };
            //for (int i = 0; i < foyer.Length; i++)
            for (int i = 0; i < popSize; i++)
            {
                //cellPopulation.cells[foyer[i]].ApicalConstriction(1f);
                //cellPopulation.cells[i].ApicalConstriction(.35f);

                float r = (float) (Helper.GaussianNumber(.945, .286) - .945);
                cellPopulation.cells[i].ExpandShrinkHexagonalCell(r);
                Console.WriteLine(r);
            }
        }

        public override void SetModel()
        {
            MGModel.dT = 0.02f;
            MGModel.Rcell = 1;
            MGModel.maximumNeighbourDistance = 3f * MGModel.Rcell;
            MGModel.DInt = 0.1f;
            MGModel.DCol = 0.05f;
            MGModel.rho = 1.0f;
            MGModel.alpha = 2f;
            MGModel.u0 = 0.2f;
            MGModel.u1 = 0.02f;
            MGModel.springForce = 1f;
            MGModel.damping = 0.5f;
            MGModel.delta = 0.001f;
            MGModel.maxExternalNeighbours = 2;
            MGModel.elasticExternalSpring = false;
            MGModel.mooreNeighbourhoodForCells = false;
            MGModel.hexagonalNeighbourhoodForCells = true;
            MGModel.staticNeighbourhood = true;

            MGModel.J = new float[nbCellTypes + 1, nbCellTypes];
            for (int i = 0; i < nbCellTypes + 1; i++)
            {
                for (int j = 0; j < nbCellTypes; j++)
                {
                    MGModel.J[i, j] = .1f;
                }
            }
            MGModel.DInteraction = new float[nbCellTypes + 1, nbCellTypes];
            MGModel.DInteraction[0, 0] = 2.5f * MGModel.DInt;
            MGModel.DInteraction[1, 0] = MGModel.DInt;
        }

        public override void Update()
        {
            while (frame++ < nbOfSimulationSteps)
            {
                if (frame % logFrequency == 0) { 
                    cellPopulation.LogState();
                    cellPopulation.LogNumbers();
                    cellPopulation.LogMetrics();
                }
                cellPopulation.DynamiseSystem();
                Console.WriteLine(frame);
            }
            cellPopulation.LogState();
            cellPopulation.LogNumbers();
            cellPopulation.LogMetrics();
            cellPopulation.LogFinalState();
            cellPopulation.LogFinalNumbers();
        }
    }
}