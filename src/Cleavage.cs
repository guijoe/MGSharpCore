using System;
using System.Collections.Generic;
using MGSharp.Core.GeometricPrimitives;
using MGSharp.Core.MGCellPopulation;
using MGSharp.Core.MGModels;
using MGSharp.Core.Helpers;
using System.Diagnostics;

namespace MGSharp
{
    class Cleavage: Simulator
    {
        public static void Main(string[] args)
        {
            Simulator simulator = new Cleavage();
            Simulator.name = "Cleavage";
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
            Console.WriteLine("Simulation ended succesfully in " + stopwatch.ElapsedMilliseconds + " milliseconds.");
            simulator.Log();

            Console.WriteLine("Simulation Logs have been written to the Log file.");
            Console.ReadKey();
        }
        
        public override void SetupSimulation()
        {
            logFrequency = 1;
            logVTK = true;
            nbOfSimulationSteps = 1000;
            popSize = 1;
            popMaxSize = 128;

            Tissue t1 = new Tissue(1, popMaxSize);
            List<Tissue> tissueList = new List<Tissue>() { t1 };
            nbCellTypes = tissueList.Count;
            cellPopulation = new CellPopulation(popSize, popMaxSize, tissueList);

            Helper.SetNCubedInN(popMaxSize);
            cellPopulation.PositionCells(false);
        }

        public override void SetInitialConditions()
        {
            name = "Cleavage";
            for (int i = 0; i < popSize; i++)
            {
                cellPopulation.cells[i].ScaleEdgeELengths(2f);
                //cellPopulation.cells[i].Default2Elliptic(2f, 2f, 2f, 7);
            }
        }

        public override void SetModel()
        {
            MGModel.dT = 0.02f;
            MGModel.Rcell = 1;
            MGModel.maximumNeighbourDistance = 2.5f * MGModel.Rcell;
            MGModel.DInt = 0f;
            MGModel.DCol = MGModel.DInt / 2;
            MGModel.rho = 1.0f;
            MGModel.alpha = 2f;
            MGModel.u0 = 0.2f;
            MGModel.u1 = 0.02f;
            MGModel.springForce = 1f;
            MGModel.damping = 0.5f;
            MGModel.delta = 0.001f;
            MGModel.divisionRate = 1.0f;
            MGModel.cellCyclePeriod = 200;
            MGModel.maxExternalNeighbours = 3;
            MGModel.elasticExternalSpring = false;
            MGModel.mooreNeighbourhoodForCells = true;
            MGModel.staticNeighbourhood = true;

            MGModel.J = new float[nbCellTypes + 1, nbCellTypes];
            for (int i = 0; i < nbCellTypes + 1; i++)
            {
                for (int j = 0; j < nbCellTypes; j++)
                {
                    MGModel.J[i, j] = .1f;
                }
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
                cellPopulation.Cleavage(true);

                Console.WriteLine(frame);
                frame++;
            }
            cellPopulation.LogState();
            cellPopulation.LogNumbers();
            cellPopulation.LogMetrics();
            cellPopulation.LogFinalState();
            cellPopulation.LogFinalNumbers();
        }
    }
}