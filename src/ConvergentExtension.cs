using System;
using System.Collections.Generic;
using MGSharp.Core.GeometricPrimitives;
using MGSharp.Core.MGCellPopulation;
using MGSharp.Core.MGModels;
using MGSharp.Core.Helpers;
using System.Diagnostics;

namespace MGSharp
{
    class ConvergentExtension: Simulator
    {
        public static void Main(string[] args)
        {
            Simulator simulator = new ConvergentExtension();
            Simulator.states = new List<PersistantVertex>();
            Simulator.numbers = new List<PersistantNumbers>();
            Simulator.finalStates = new List<PersistantVertex>();
            Simulator.finalNumbers = new List<PersistantNumbers>();

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
            nbCellTypes = 1;
            nbOfSimulationSteps = 1000;
            popSize = 25;
            popMaxSize = 125;

            //Cylinder34 mesh = new Cylinder34();

            //HexagonalCylinder42 mesh = new HexagonalCylinder42(1f);

            //Console.WriteLine(mesh.innerRadius);
            Tissue t1 = new Tissue(25, popMaxSize, new Cylinder34(new Vector(.5f, .5f, .5f)));
            List<Tissue> tissueList = new List<Tissue>() { t1 };

            nbCellTypes = tissueList.Count;
            cellPopulation = new CellPopulation(popSize, popMaxSize, tissueList);

            Helper.SetNCubedInN(popMaxSize);
            cellPopulation.PositionCells(false);
        }

        public override void SetInitialConditions()
        {
            name = "Convergent Extension";
            for (int i = 0; i < popSize; i++)
            {
                cellPopulation.cells[i].Cylinder34Squamous2Columnar(2f);
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
            while (frame++ < nbOfSimulationSteps)
            {
                cellPopulation.LogState();
                cellPopulation.LogNumbers();
                cellPopulation.DynamiseSystem();
                Console.WriteLine(frame);
            }
            cellPopulation.LogState();
            cellPopulation.LogNumbers();
            cellPopulation.LogFinalState();
            cellPopulation.LogFinalNumbers();
        }
    }
}