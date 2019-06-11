using MGSharp.Core.GeometricPrimitives;
using System;
using System.Collections.Generic;
using System.IO;

namespace MGSharp.Core.Helpers
{
    public class CellProvider
    {
        public int popSize;
        public int maxPopSize;
        public int nbOfCellTypes;
        public int nbOfFrames;
        int nbOfLinesPerFrameStates;

        public int[] nbOfCellsPerType;
        int[] nbOfVerticesPerType;
        int[][] nbOfVerticesPerTypePerCell;
        int[] nbOfVerticesPerCell;
        int[] nbOfTrianglesPerType;
        Mesh[] meshPerCellType;
        public Mesh[] meshPerCell;
        int[][] trianglesPerCellType;
        int[] cellTypePerCell;
        int cellIndex=0;

        int[] popSizePerFrame;
        int[] startOfFrame;
        public int[] nbOfCellTypesPerFrame;
        public int[][] nbOfCellsPerTypePerFrame;
        int[][][] nbOfVerticesPerTypePerCellPerFrame;
        int[][] nbOfVerticesPerCellPerFrame;
        Vector[][][][] positionsPerTypePerCellPerVerticePerFrame;

        Vector[] positions;
        //IList<PersistantVertex> states;
        Vector[] numbers;
        //int frameStartStates = 0;

        public int frame;
        public string folder;

        public void ReadParameters(string dir)
        {
            folder = dir;
            string parFile = folder + "/Parameters.mg";

            string[] parameters = File.ReadAllLines(parFile);

            int currentLine = 1;
            int.TryParse(parameters[currentLine++].Split('=')[1], out nbOfFrames);
            nbOfFrames = 1;

            while (parameters[currentLine++] != string.Empty) { }
            currentLine++;

            int.TryParse(parameters[currentLine++].Split('=')[1], out popSize);
            int.TryParse(parameters[currentLine++].Split('=')[1], out maxPopSize);
            int.TryParse(parameters[currentLine++].Split('=')[1], out nbOfCellTypes);

            nbOfCellsPerType = new int[nbOfCellTypes];
            for (int i = 0; i < nbOfCellTypes; i++)
            {
                int.TryParse(parameters[currentLine++].Split('=')[1], out nbOfCellsPerType[i]);
            }

            while (parameters[currentLine++] != string.Empty) { }

            currentLine += nbOfCellTypes * (nbOfCellTypes + 1) / 2 + nbOfCellTypes + 4;

            nbOfVerticesPerType = new int[nbOfCellTypes];
            nbOfTrianglesPerType = new int[nbOfCellTypes];
            meshPerCellType = new Mesh[nbOfCellTypes];
            trianglesPerCellType = new int[nbOfCellTypes][];

            int step = 0, pos = currentLine;
            Vector[] vs;
            int[] tris;
            for (int i = 0; i < nbOfCellTypes; i++)
            {
                int.TryParse(parameters[step + pos].Split('=')[1], out nbOfVerticesPerType[i]);
                int.TryParse(parameters[step + pos + 1].Split('=')[1], out nbOfTrianglesPerType[i]);

                currentLine += 4;
                vs = new Vector[nbOfVerticesPerType[i]];
                for (int j = 0; j < nbOfVerticesPerType[i]; j++)
                {
                    string[] xyz = parameters[currentLine++]
                        .Split(',');

                    vs[j] = new Vector(float.Parse(xyz[0]), float.Parse(xyz[1]), float.Parse(xyz[2]));
                }
                currentLine += 2;
                tris = new int[3 * nbOfTrianglesPerType[i]];
                for (int j = 0; j < nbOfTrianglesPerType[i]; j++)
                {
                    string[] xyz = parameters[currentLine++].Split(',');
                    tris[3 * j] = int.Parse(xyz[0]);
                    tris[3 * j + 1] = int.Parse(xyz[1]);
                    tris[3 * j + 2] = int.Parse(xyz[2]);

                }
                currentLine += 3;
                step += nbOfVerticesPerType[i] + nbOfTrianglesPerType[i] + 9;

                meshPerCellType[i] = new Mesh();
                //meshPerCellType[i].vertices = vs;
                //meshPerCellType[i].triangles = tris;
                trianglesPerCellType[i] = tris;
            }
        }

        public void ReadResults()
        {
            string statesFile = folder + "/Final_States.mg";
            string numbersFile = folder + "/Final_Numbers.mg";

            string[] results = File.ReadAllLines(statesFile);
            positions = new Vector[results.Length - 1];

            string[] counts = File.ReadAllLines(numbersFile);
            numbers = new Vector[counts.Length - 1];

            startOfFrame = new int[nbOfFrames];
            popSizePerFrame = new int[nbOfFrames];
            nbOfCellTypesPerFrame = new int[nbOfFrames];
            nbOfCellsPerTypePerFrame = new int[nbOfFrames][];
            nbOfVerticesPerCellPerFrame = new int[nbOfFrames][];
            nbOfVerticesPerTypePerCellPerFrame = new int[nbOfFrames][][];
            positionsPerTypePerCellPerVerticePerFrame = new Vector[nbOfFrames][][][];


            int count = 0;
            int currentLineNumbers = 1;
            int frStart = 0;
            int populationSize = 0;
            for (int fr = 0; fr < nbOfFrames; fr++)
            {
                nbOfLinesPerFrameStates = 0;
                for (int i = 0; i < populationSize; i++)
                {
                    nbOfLinesPerFrameStates += nbOfVerticesPerCell[i] + 1;
                }
                if (fr > 0)
                {
                    frStart += nbOfLinesPerFrameStates;
                    startOfFrame[fr] = frStart;
                }

                populationSize = int.Parse(counts[currentLineNumbers].Split(';')[2]);
                nbOfCellTypes = int.Parse(counts[currentLineNumbers].Split(';')[1]);
                nbOfCellsPerType = new int[nbOfCellTypes];
                nbOfVerticesPerTypePerCell = new int[nbOfCellTypes][];
                nbOfVerticesPerCell = new int[populationSize];

                popSizePerFrame[fr] = int.Parse(counts[currentLineNumbers].Split(';')[2]);
                popSize = popSizePerFrame[fr];
                cellTypePerCell = new int[popSize];
                nbOfCellTypesPerFrame[fr] = int.Parse(counts[currentLineNumbers].Split(';')[1]);
                nbOfCellsPerTypePerFrame[fr] = new int[nbOfCellTypesPerFrame[fr]];
                nbOfVerticesPerTypePerCellPerFrame[fr] = new int[nbOfCellTypesPerFrame[fr]][];
                nbOfVerticesPerCellPerFrame[fr] = new int[popSizePerFrame[fr]];
                positionsPerTypePerCellPerVerticePerFrame[fr] = new Vector[nbOfCellTypesPerFrame[fr]][][];

                currentLineNumbers++;

                int cellTypeStartStates = 1;
                int cellTypeStartNumbers = 0;
                for (int i = 0; i < nbOfCellTypes; i++)
                {
                    if (i > 0)
                    {
                        for (int j = 0; j < nbOfCellsPerType[i - 1]; j++)
                        {
                            cellTypeStartStates += nbOfVerticesPerTypePerCell[i - 1][j] + 1;
                        }
                    }
                    if (i > 0) cellTypeStartNumbers += nbOfCellsPerType[i - 1];

                    nbOfCellsPerType[i] = int.Parse(counts[currentLineNumbers].Split(';')[2]);
                    nbOfVerticesPerTypePerCell[i] = new int[nbOfCellsPerType[i]];
                    currentLineNumbers++;

                    nbOfCellsPerTypePerFrame[fr][i] = nbOfCellsPerType[i];
                    nbOfVerticesPerTypePerCellPerFrame[fr][i] = nbOfVerticesPerTypePerCell[i];
                    positionsPerTypePerCellPerVerticePerFrame[fr][i] = new Vector[nbOfCellsPerType[i]][];

                    for (int j = 0; j < nbOfCellsPerType[i]; j++)
                    {
                        string[] splitLine = counts[currentLineNumbers].Split(';');
                        if(splitLine.Length == 4)
                        {
                            cellTypePerCell[cellIndex] = int.Parse(splitLine[3]);
                        }
                        else
                        {
                            cellTypePerCell[cellIndex] = i;
                        }

                        nbOfVerticesPerTypePerCell[i][j] = int.Parse(counts[currentLineNumbers].Split(';')[2]);
                        nbOfVerticesPerCell[cellTypeStartNumbers + j] = nbOfVerticesPerTypePerCell[i][j];
                        currentLineNumbers++;

                        nbOfVerticesPerTypePerCellPerFrame[fr][i][j] = nbOfVerticesPerTypePerCell[i][j];
                        nbOfVerticesPerCellPerFrame[fr][cellTypeStartNumbers + j] = nbOfVerticesPerCell[cellTypeStartNumbers + j];
                        positionsPerTypePerCellPerVerticePerFrame[fr][i][j] = new Vector[nbOfVerticesPerTypePerCell[i][j] + 1];

                        int cellStart = frStart + cellTypeStartStates;
                        if (j > 0)
                        {
                            cellStart += nbOfVerticesPerTypePerCell[i][j - 1] + 1;
                        }
                        cellStart = count + 1;

                        for (int k = 0; k < nbOfVerticesPerTypePerCell[i][j] + 1; k++)
                        {
                            string vec = results[cellStart + k].Split(';')[3];
                            string[] xyz = vec
                                                .Split(',');
                            
                            positions[cellStart + k - 1] = new Vector(float.Parse(xyz[0]),
                                                                    float.Parse(xyz[1]),
                                                                     float.Parse(xyz[2]));

                            positionsPerTypePerCellPerVerticePerFrame[fr][i][j][k] = positions[cellStart + k - 1];
                            count++;
                        }
                        cellIndex++;
                    }
                }
            }
        }

        int lineNumber = 0;
        int cellLineNumber = 0;
        public void PlayResults(int frame)
        {
            meshPerCell = new Mesh[popSizePerFrame[frame]];
            //Console.WriteLine(popSizePerFrame[frame]);
            int frameStartStates = nbOfLinesPerFrameStates * frame;
            int cStart = 0;

            Vector[] vs;

            for(int i=0; i<popSizePerFrame[frame]; i++)
            {

            }

            for (int i = 0; i < nbOfCellTypesPerFrame[frame]; i++)
            {
                if (i > 0) cStart += nbOfCellsPerType[i - 1];
                Console.WriteLine("Nb cells paer Tissue :" + nbOfCellsPerTypePerFrame[frame][i]);
                for (int j = 0; j < nbOfCellsPerTypePerFrame[frame][i]; j++)
                {
                    //Console.WriteLine(nbOfVerticesPerTypePerCellPerFrame[frame][i][j]);
                    vs = new Vector[nbOfVerticesPerTypePerCellPerFrame[frame][i][j]];
                    cellLineNumber++;

                    meshPerCell[cStart + j] = new Mesh();

                    for (int k = 0; k < vs.Length; k++)
                    {
                        vs[k] = positionsPerTypePerCellPerVerticePerFrame[frame][i][j][k];
                        lineNumber++;

                        meshPerCell[cStart + j].AddVertex(vs[k]);
                    }

                    for (int k = 0; k < trianglesPerCellType[i].Length; k+=3)
                    {
                        meshPerCell[cStart + j].AddTriangle(trianglesPerCellType[i][k],
                                                            trianglesPerCellType[i][k+1],
                                                            trianglesPerCellType[i][k+2]);
                    }

                    /*
                    for (int k = 0; k < trianglesPerCellType[cellTypePerCell[cStart + j]].Length; k += 3)
                    {
                        meshPerCell[cStart + j].AddTriangle(trianglesPerCellType[cellTypePerCell[cStart + j]][k],
                                                            trianglesPerCellType[cellTypePerCell[cStart + j]][k + 1],
                                                            trianglesPerCellType[cellTypePerCell[cStart + j]][k + 2]);
                    }
                    */
                }
            }
        }
    }
}