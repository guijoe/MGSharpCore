using System;
using System.Collections.Generic;
using MGSharp.Core.GeometricPrimitives;
using MGSharp.Core.MGCellPopulation;
using MGSharp.Core.MGModels;
using MGSharp.Core.Helpers;
using System.Diagnostics;

namespace MGSharp
{
    /*
    class Invagination: Simulator
    {
        public static void Main(string[] args) { }

        public override void SetupSimulation(){}

        public override void SetInitialConditions(){}

        public override void Update(){}

        public override void SetModel(){}
    }
    */

    class Invagination: Simulator
    {
        public static void Main(string[] args)
        {
            Simulator simulator = new Invagination();
            Simulator.name = "Invagination";
            Simulator.states = new List<PersistantVertex>();
            Simulator.numbers = new List<PersistantNumbers>();
            Simulator.finalStates = new List<PersistantVertex>();
            Simulator.finalNumbers = new List<PersistantNumbers>();
            Simulator.metrics = new List<PersistantMetrics>();

            simulator.SetupSimulation();
            simulator.SetModel();
            simulator.SetInitialConditions();

            simulator.LogParameters();

            Stopwatch stopwatch = new Stopwatch();
            Console.WriteLine("Simulation in Progress ...");
            stopwatch.Start();
            simulator.Update();
            stopwatch.Stop();
            Console.WriteLine("Simulation ended succesfully in " 
                                + stopwatch.ElapsedMilliseconds 
                                    + " milliseconds.");
            simulator.Log();

            Console.WriteLine("Simulation Logs have been written to Log files.");
            Console.ReadKey();
        }
        
        public override void SetupSimulation()
        {
            logFrequency = 10;
            nbOfSimulationSteps = 10000;
            popSize = 81;
            popMaxSize = 729;

            HexagonalCylinder42 mesh = new HexagonalCylinder42((float)Math.Sqrt(3) / 3);

            Tissue t1 = new Tissue(81, popMaxSize, mesh);
            List<Tissue> tissueList = new List<Tissue>() { t1 };

            nbCellTypes = tissueList.Count;
            cellPopulation = new CellPopulation(popSize, popMaxSize, tissueList);

            Helper.SetNCubedInHexa2(popMaxSize);
            cellPopulation.PositionCells(false);
        }

        public override void SetInitialConditions()
        {
            int[] foyer = new int[] { 20, 21, 22, 23, 24, 29, 30, 31, 32, 33, 38, 39,
                                        40, 41, 42, 47, 48, 49, 50, 51, 56, 57, 58, 59, 60 };

            for (int i = 0; i < foyer.Length; i++)
            {
                cellPopulation.cells[foyer[i]].ApicalConstriction(0.5f);
            }
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
            cellPopulation.LogFinalState();
            cellPopulation.LogFinalNumbers();
        }

        public override void SetModel()
        {
            MGModel.dT = 0.02f;
            MGModel.Rcell = 1f;
            MGModel.maximumNeighbourDistance = 3f * MGModel.Rcell;
            MGModel.DInt = 0.1f;
            MGModel.DCol = 0.05f;
            MGModel.rho = 1.0f;
            MGModel.alpha = 2f;
            MGModel.u0 = 0.2f;
            MGModel.u1 = 0.02f;
            MGModel.damping = 0.5f;
            MGModel.delta = 0.001f;
            MGModel.maxExternalNeighbours = 3;
            MGModel.elasticExternalSpring = false;
            MGModel.mooreNeighbourhoodForCells = true;
            MGModel.staticNeighbourhood = true;

            MGModel.J = new float[nbCellTypes + 1, nbCellTypes];
            MGModel.DInteraction = new float[nbCellTypes + 1, nbCellTypes];
            for (int i = 0; i < nbCellTypes + 1; i++)
            {
                for (int j = 0; j < nbCellTypes; j++)
                {
                    MGModel.J[i, j] = 1f;
                    MGModel.DInteraction[i, j] = MGModel.DInt;
                }
            }
        }
    }
}