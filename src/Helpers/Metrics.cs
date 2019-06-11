using System.Collections;
using System.Collections.Generic;
using MGSharp.Core.MGCellModels;

namespace MGSharp.Core.Helpers
{
	public class Metrics
	{
        MGCell[] cells;
		public static float ElasticEnergy()
        {
            float E=0;
            for (int i=0; i<Simulator.cellPopulation.populationSize; i++)
            {
                E += Simulator.cellPopulation.cells[i].ElasticEnergyNucleusRays() 
                    + Simulator.cellPopulation.cells[i].ElasticEnergyMembraneRays();
            }
            return E;
        }

        public static float MorseEnergy()
        {
            float E = 0;
            for (int i = 0; i < Simulator.cellPopulation.populationSize; i++)
            {
                E += Simulator.cellPopulation.cells[i].MorseEnergyNucleusRays()
                    + Simulator.cellPopulation.cells[i].MorseEnergyMembraneRays();
            }
            return E;
        }

        public static float ElasticEnergyWithExternalLinks()
        {
            float E = 0;
            for (int i = 0; i < Simulator.cellPopulation.populationSize; i++)
            {
                //E += Simulator.cellPopulation.cells[i].ElasticEnergyNucleusRays() 
                //    + Simulator.cellPopulation.cells[i].ElasticEnergyMembraneRays() 
                //    + Simulator.cellPopulation.cells[i].ElasticEnergyAdhesionRays();
            }
            return E;
        }
    }
}