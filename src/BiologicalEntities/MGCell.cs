using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using MGSharp.Core.GeometricPrimitives;
using MGSharp.Core.Helpers;
using MGSharp.Core.MGModels;

namespace MGSharp.Core.MGCellModels
{
    public class MGCell : Mesh
    {
        public int cellId;
        public int parent;
        public int tissueId;
        public int nbOfParticles;
        public Vertex[] nuclei;
        public Vector centre;
        
        public List<int> neighbours;
        public EdgeSet externalEdges;
        public EdgeSet nucleusEdges;

        public Vector[] targetVertices;
        private Vector[] potentialNucleiPos;
        public float[,] particlesDistances;
        public int[] elongationAxis;

        public bool inGrowMode = true;
        public bool inDivisionMode = false;
        public bool hasDivided = false;
        public bool isNewBorn = true;

        public int cyclePeriod = 200;
        public int cycleTime = 0;

        public int appliedForces;
        public float Rcell;

        public float[] spins;
        public int[] sigma;
        public Vector polarisation;

        public string tissueName;
        //public int tissueId;

        public MGCell()
        {
            nbOfParticles = vertexCount();
        }

        public MGCell(int i)
        {
            cellId = i;
            nuclei = new Vertex[2];
            nuclei[0] = new Vertex();
            centre = new Vector(0, 0, 0);
            
            nbOfParticles = vertexCount();
            innerRadius = new Vector(1, 1, 1);
            neighbours = new List<int>();

            sigma = new int[nbOfParticles];
            spins = new float[nbOfParticles];
            targetVertices = new Vector[nbOfParticles];
            subset = new List<int>();

            externalEdges = new EdgeSet();
            nucleusEdges = new EdgeSet();

            nuclei[0].id = nbOfParticles;
            nuclei[0].cellId = cellId;
            for (int j = 0; j < nbOfParticles; j++)
            {
                vertices[j].id = j;
                vertices[j].cellId = i;
                Edge edge = new Edge(vertices[j], nuclei[0]);
                nucleusEdges.add(edge);
                targetVertices[j] = vertices[j].Clone().v;

                int k = new Random().Next(2);
                spins[j] = (k == 0) ? MGModel.delta : -MGModel.delta;
                sigma[j] = j;
                vertices[j].externalNeighbours = new List<int[]>();
                vertices[j].globalForces = new Vector();
                subset.Add(j);
            }

            nuclei[0].v = ComputeCentreFromMesh();
            SetEdgeELengths();
            ResetCell();
        }

        public MGCell(int i, Mesh m)
        {
            cellId = i;
            nuclei = new Vertex[2];
            nuclei[0] = new Vertex();
            centre = new Vector(0, 0, 0);
            
            ShapeCell(m);
            nbOfParticles = vertexCount();
            innerRadius = m.innerRadius;
            neighbours = new List<int>();

            sigma = new int[nbOfParticles];
            spins = new float[nbOfParticles];
            targetVertices = new Vector[nbOfParticles];

            externalEdges = new EdgeSet();
            nucleusEdges = new EdgeSet();
            subset = new List<int>();

            nuclei[0].id = nbOfParticles;
            nuclei[0].cellId = cellId;
            for (int j = 0; j < nbOfParticles; j++)
            {
                vertices[j].id = j;
                vertices[j].cellId = i;
                Edge edge = new Edge(vertices[j], nuclei[0]);
                nucleusEdges.add(edge);
                targetVertices[j] = vertices[j].Clone().v;

                int k = new Random().Next(2);
                spins[j] = (k == 0) ? MGModel.delta : -MGModel.delta;
                sigma[j] = j;
                vertices[j].externalNeighbours = new List<int[]>();
                vertices[j].globalForces = new Vector();
                subset.Add(j);
            }

            nuclei[0].v = ComputeCentreFromMesh();
            SetEdgeELengths();
            ResetCell();
        }
        
        public void PositionCell(Vector c) {
            centre.Copy(c);
            nuclei[0].v.Copy(c);

            for (int i = 0; i < nbOfParticles; i++)
            {
                vertices[i].translate(c);
            }
        }

        public void Translate(Vector u)
        {
            for (int i = 0; i < nbOfParticles; i++)
            {
                vertices[i].translate(u);
            }
            nuclei[0].translate(u);
        }
        
        public void ShapeCell(Mesh m)
        {
            vertices = new VertexSet();
            edges = new EdgeSet();
            faces = new FaceSet();

            copy(m);
        }

        public void ComputeForces()
        {
            ResetCell();
            
            for (int i = 0; i < edgeCount(); i++)
            {
                //edges[i].force = MGModel.Force(edges[i].length(), MGModel.u0, edges[i].l0);
                edges[i].force = MGModel.Force(edges[i].length(), MGModel.J[tissueId, tissueId], edges[i].l0);
            }
            for (int i = 0; i < nucleusEdges.getCount(); i++)
            {
                //nucleusEdges[i].force = MGModel.Force(nucleusEdges[i].length(), MGModel.u0, nucleusEdges[i].l0);
                nucleusEdges[i].force = MGModel.Force(nucleusEdges[i].length(), MGModel.J[tissueId, tissueId], nucleusEdges[i].l0);
            }
            
            for (int i = 0; i < edgeCount(); i++)
            {
                edges[i].ends[0].internalForces += edges[i].force * edges[i].UnitVector();
                edges[i].ends[1].internalForces -= edges[i].force * edges[i].UnitVector();
            }
            for (int i = 0; i < nucleusEdges.getCount(); i++)
            {
                nucleusEdges[i].ends[0].nucleusForce0 += nucleusEdges[i].force * nucleusEdges[i].UnitVector();
                nucleusEdges[i].ends[1].force -= nucleusEdges[i].force * nucleusEdges[i].UnitVector();
            }
        }

        public void ComputeExternalForces()
        {
            if (MGModel.elasticExternalSpring)
            {
                /*
                Parallel.For(0, externalEdges.getCount(), i =>
                {
                    externalEdges[i].force = MGModel.Force(externalEdges[i].length(), MGModel.u1, externalEdges[i].l0);
                });
                */

                for (int i = 0; i < externalEdges.getCount(); i++)
                {
                    //externalEdges[i].force = MGModel.Force(externalEdges[i].length(), MGModel.u1, externalEdges[i].l0);
                    externalEdges[i].force = MGModel.Force(externalEdges[i].length(), 5f, MGModel.DCol);
                }

                for (int i = 0; i < externalEdges.getCount(); i++)
                {
                    externalEdges[i].ends[0].externalForces += externalEdges[i].force * externalEdges[i].UnitVector();
                    externalEdges[i].ends[1].externalForces -= externalEdges[i].force * externalEdges[i].UnitVector();
                }
            }
            else
            {
                for (int i = 0; i < externalEdges.getCount(); i++)
                {
                    //if (!externalEdges[i].ends[0].nullForces)
                        externalEdges[i].ends[0].externalForces += externalEdges[i].ends[1].internalForces + externalEdges[i].ends[1].nucleusForce0 + externalEdges[i].ends[1].globalForces;
                    //else
                    //    externalEdges[i].ends[0].externalForces += externalEdges[i].ends[1].internalForces + externalEdges[i].ends[1].nucleusForce0;
                }
            }
            
            for (int i = 0; i < nbOfParticles; i++)
            {
                if (!vertices[i].nullForces)
                    vertices[i].force = vertices[i].internalForces + vertices[i].nucleusForce0 + vertices[i].externalForces + vertices[i].globalForces;
                //Console.WriteLine(cellId + ", " + i + ", " + vertices[i].internalForces + "; " + vertices[i].nucleusForce0 + "; " + vertices[i].externalForces);
            }
        }

        public void Dynamise()
        {
            for(int i=0; i< nbOfParticles; i++) {
                vertices[i].Move();
            }
            nuclei[0].Move();
            cycleTime++;
        }

        public void CellCycle(int appliedForces, Vector normal, bool staticShape)
        {
            //Compute forces applied on cell
            //ChangeShape();
            //Console.WriteLine("Cycle");
            ComputeForces();

            //Cell Growth
            if (inGrowMode)
            {
                Dynamise();
            }

            //Check if cell is ready for division
            if (cycleTime == MGModel.cellCyclePeriod && !hasDivided)
            {
                inDivisionMode = true;
                inGrowMode = false;
            }

            //Cell Division
            if (inDivisionMode && Simulator.cellPopulation.populationSize < Simulator.cellPopulation.maxPopulationSize)
            {
                //Console.WriteLine("Potentially Dividing ...");
                double r = new Random().NextDouble();
                if(r < MGModel.divisionRate)
                {
                    //Console.WriteLine("Dividing ...");
                    Mitosis(normal, staticShape);
                    inDivisionMode = false;
                    MGModel.searchNeighbours = true;
                }
            }

            //Initialize new cells
            if (isNewBorn && hasDivided)
            {
                inGrowMode = true;
                isNewBorn = false;
                hasDivided = false;
            }
        }

        public void KleinCellCycle(int appliedForces, Vector normal, bool staticShape) {
            //Check if cell is ready for division
            if (cycleTime==MGModel.cellCyclePeriod && !hasDivided)
            {
                inDivisionMode = true;
                inGrowMode = false;
                //cycleTime = 0;
            }
            //Console.WriteLine("Tissue: " + tissueName + ", cycleTime: " + cycleTime);
            //Cell Division
            //if (inDivisionMode && Simulator.cellPopulation.populationSize < Simulator.cellPopulation.maxPopulationSize)
            if (inDivisionMode && Simulator.cellPopulation.tissues[tissueId].populationSize < Simulator.cellPopulation.tissues[tissueId].maxPopulationSize )
            {
                double r = new Random().NextDouble();
                if (r <= MGModel.divisionRate)
                {
                    //
                    Mitosis(normal, staticShape);
                    inDivisionMode = false;
                    MGModel.searchNeighbours = true;
                    MGModel.nextNeighboursSearchFrame = Simulator.frame + 1;
                }
            }

            //Initialize new cells
            if (isNewBorn && hasDivided)
            {
                inGrowMode = true;
                isNewBorn = false;
                hasDivided = false;
            }
        }

        public void Mitosis(Vector normal, bool staticShape)
        {
            //Console.WriteLine("Dividing Cell: " + cellId);
            //Console.WriteLine("Face Count: " + faceCount());

            Plane splitPlane = CreateCutPlane(normal);
            MeshSplitter splitter = new MeshSplitter(this, splitPlane);
            splitter.MeshInitialize();
            splitter.MeshSplit();
            Mesh lowerMesh = splitter.CreateMeshLower();
            Mesh upperMesh = splitter.CreateMeshUpper();
            
            MGCell newCell = new MGCell(Simulator.cellPopulation.populationSize, lowerMesh);
            if (newCell != null)
            {
                newCell.cellId = Simulator.cellPopulation.populationSize;
                newCell.appliedForces = appliedForces;
                newCell.Rcell = Rcell;

                newCell.inDivisionMode = false;
                newCell.inGrowMode = true;
                newCell.isNewBorn = true;
                newCell.parent = cellId;
                newCell.cyclePeriod = 200;
                newCell.cycleTime = 0;
                newCell.nuclei[0].v = newCell.ComputeCentreFromMesh();
                newCell.centre = newCell.ComputeCentreFromMesh();
                newCell.spins = spins;
                newCell.polarisation = polarisation;
                subset = new List<int>();

                for(int i = 0; i<newCell.nbOfParticles; i++)
                {
                    newCell.vertices[i].globalForces = new Vector();
                    newCell.subset.Add(i);
                    subset.Add(i);
                }
                Simulator.cellPopulation.AddCell(tissueId, newCell);
            }

            //Booleans
            inDivisionMode = false;
            hasDivided = true;
            isNewBorn = true;
            inGrowMode = true;
            cycleTime = 0;
            //cyclePeriod = 100;

            if (staticShape)
            {
                SetEdgeELengths();
            }
            else
            {
                //SetEdgeELengths();
                if(polarisation == Vector.down)
                {
                    SetElengths(MGModel.modelMeshes[1]);
                    newCell.SetElengths(MGModel.modelMeshes[1]);
                }
                else
                {
                    SetElengths(MGModel.modelMeshes[0]);
                    newCell.SetElengths(MGModel.modelMeshes[0]);
                }
                
            }

            copy(upperMesh);
            for (int i = 0; i < nbOfParticles; i++)
            {
                vertices[i].globalForces = new Vector();
                subset.Add(i);
            }
            //SetEdgeELengths();
            nuclei[0].v = ComputeCentreFromMesh();
            centre = ComputeCentreFromMesh();
            for (int i = 0; i < nbOfParticles; i++)
            {
                targetVertices[i] = vertices[i].Clone().v;
            }

            if (!staticShape)
            {
                if (polarisation == Vector.down)
                {
                    SetElengths(MGModel.modelMeshes[1]);
                    newCell.SetElengths(MGModel.modelMeshes[1]);
                }
                else
                {
                    SetElengths(MGModel.modelMeshes[0]);
                    newCell.SetElengths(MGModel.modelMeshes[0]);
                }
            }
        }

        public void ChangeShape()
        {
            double dist = 0;
            float[] minDist = new float[nbOfParticles];
            int jAfter = 0;
            int jCurrent = 0;
            Vector meanPoint = centre;
            Vector norm = Vector.zero;
            Vector tang = Vector.zero;
            Vector tP;

            Vector spinSiteCoordinatesInCurrentConfig, spinSiteCoordinatesAfterFlip;
            double energyContributionInCurrentConfig = 100f, energyContributionAfterFlip = 100f;
            Vector[] currentPos = new Vector[nbOfParticles];

            int spinSite = new Random().Next(nbOfParticles);

            for (int i = 0; i < spinSite; i++)
            {
                norm = vertices[i].GetPosition() - meanPoint;
                currentPos[i] = (nucleusEdges[i].l0 + spins[i]) * (norm / norm.norm()) + meanPoint;
            }
            for (int i = spinSite + 1; i < nbOfParticles; i++)
            {
                norm = vertices[i].GetPosition() - meanPoint;
                currentPos[i] = (nucleusEdges[i].l0 + spins[i]) * (norm / norm.norm()) + meanPoint;
            }

            norm = vertices[spinSite].GetPosition() - meanPoint;
            spinSiteCoordinatesInCurrentConfig = (nucleusEdges[spinSite].l0 + spins[spinSite]) * (norm / norm.norm()) + meanPoint;
            spinSiteCoordinatesAfterFlip = (nucleusEdges[spinSite].l0 - spins[spinSite]) * (norm / norm.norm()) + meanPoint;

            for (int j = 0; j < nbOfParticles; j++)
            {
                dist = Vector.Distance(spinSiteCoordinatesInCurrentConfig, targetVertices[j]);
                if (dist < energyContributionInCurrentConfig)
                {
                    energyContributionInCurrentConfig = dist;
                    jCurrent = j;
                }
            }

            for (int j = 0; j < nbOfParticles; j++)
            {
                dist = Vector.Distance(spinSiteCoordinatesAfterFlip, targetVertices[j]);
                if (dist < energyContributionAfterFlip)
                {
                    energyContributionAfterFlip = dist;
                    jAfter = j;
                }
            }

            double energyGap = energyContributionAfterFlip - energyContributionInCurrentConfig;

            if (energyGap < 0)
            {
                spins[spinSite] = -spins[spinSite];
                currentPos[spinSite] = spinSiteCoordinatesAfterFlip;
                tP = targetVertices[jAfter];
            }
            else
            {
                currentPos[spinSite] = spinSiteCoordinatesInCurrentConfig;
                tP = targetVertices[jCurrent];
                if (MGModel.T > 0)
                {
                    currentPos[spinSite] = spinSiteCoordinatesInCurrentConfig;
                    double p = Math.Exp(-energyGap / MGModel.T);
                    double r = new Random().NextDouble();
                    if (r <= p)
                    {
                        spins[spinSite] = -spins[spinSite];
                        currentPos[spinSite] = spinSiteCoordinatesAfterFlip;
                        tP = targetVertices[jAfter];
                    }
                }
            }

            for (int i = 0; i < nbOfParticles; i++)
            {
                nucleusEdges[i].l0 += spins[i];
            }

            /*
            norm = currentPos[spinSite] - meanPoint;
            Plane plane = new Plane(currentPos[spinSite], norm);
            Vector direction = (tP - currentPos[spinSite]).normalized;
            tang = plane.PointOrthogonalProjection(direction);

            float coef = 0, max = 0; int jMax = 0;
            for (int j = 0; j < particleNeighbourhoods[spinSite].Count; j++)
            {
                coef = Vector.Dot(tang, vertices[particleNeighbourhoods[spinSite][j]].GetPosition() - currentPos[spinSite]);
                if (Math.Abs(coef) > Math.Abs(max))
                {
                    max = coef;
                    jMax = j;
                }
            }
            //if(Mathf.Abs(max) - MGModel.delta > MGModel.epsilon)
            //{
            //cellReq[spinSite][jMax] -= Mathf.Sign(max) * 0.01f;
            //Debug.Log(cellReq[spinSite][jMax]);
            //}
            */
        }

        public void Reshape(Mesh mesh)
        {
            float[,] costs = new float[vertexCount(), mesh.vertexCount()];
            for(int i=0; i<vertexCount(); i++)
            {
                for(int j=0; j<mesh.vertexCount(); j++)
                {
                    costs[i,j] = (float)Vector.Distance(vertices[i].v, mesh.vertices[j].v);
                }
            }

            int[] mapping = HungarianAlgorithm.FindAssignments(costs);

            for(int i=0; i<mapping.Length; i++)
            {
                vertices[i].v = mesh.vertices[mapping[i]].v;
            }

            //Console.WriteLine(mapping.Length);
        }

        public void ApicalConstriction(float r)
        {
            Vector meanPoint = Vector.zero;
            for (int i = 0; i < nbOfParticles; i++)
            {
                meanPoint += vertices[i].v;
            }
            meanPoint /= nbOfParticles;

            double coef;
            Vector O2;
            for (int i = 0; i < nbOfParticles; i++)
            {
                targetVertices[i] = vertices[i].v - meanPoint;

                O2 = new Vector(0, targetVertices[i].y);

                double xz = Math.Sqrt(targetVertices[i].x * targetVertices[i].x + targetVertices[i].z * targetVertices[i].z);
                //double MM2 = (r * innerRadius.x) * (targetVertices[i].y + innerRadius.y);
                double MM2 = (r * xz) * (targetVertices[i].y + innerRadius.y);
                Vector MO2 = O2 - targetVertices[i];
                coef = MO2.norm() < MGModel.epsilon ? 1 : 1 - MM2 / MO2.norm();

                Vector M2O2 = coef * MO2;
                Vector M2 = O2 - M2O2;
                targetVertices[i] = M2;    
            }

            if (polarisation == Vector.down)
            {
                for (int i = 0; i < nbOfParticles; i++)
                {
                    targetVertices[i].y = -targetVertices[i].y;
                    vertices[i].v = 2 * meanPoint - vertices[i].v;
                }
            }

            for (int i = 0; i < nbOfParticles; i++)
            {
                targetVertices[i] += meanPoint;
                //vertices[i].v = targetVertices[i];
                nucleusEdges[i].l0 = (targetVertices[i] - meanPoint).norm();
            }

            for (int j = 0; j < edgeCount(); j++)
            {
                edges[j].l0 = Vector.Distance(targetVertices[edges[j].ends[0].pos], targetVertices[edges[j].ends[1].pos]);
            }
        }

        public void SetElengths(Mesh cell)
        {
            Vector cellCentre = cell.ComputeCentreFromMesh();
            for (int i = 0; i < nbOfParticles; i++)
            {
                nucleusEdges[i].l0 = Vector.Distance(cellCentre, cell.vertices[i].v);
            }

            for (int j = 0; j < edgeCount(); j++)
            {
                edges[j].l0 = cell.edges[j].l0;
            }
        }

        public void ApicalConstriction1(float r)
        {
            //this.appliedForces = forces;
            //Debug.Log("Apical Constriction");
            Vector meanPoint = Vector.zero;
            for (int i = 0; i < nbOfParticles; i++)
            {
                meanPoint += vertices[i].GetPosition();
            }
            meanPoint /= nbOfParticles;

            float R = (float) innerRadius.x;
            float H = (float) innerRadius.y * 2f;

            //Console.WriteLine(R + ", " + r + ", " + H);

            float delta = 3 * (3 * R - r) * (R + r);
            float t1 = (float) (Math.Sqrt(delta) - (3 * R - r)) / 2f;
            //Console.WriteLine(t1);

            //float t2 = (-Mathf.Sqrt(delta) - (3 * R - r)) / 2f;
            //float zero = t1 * t1 + (3 * R - r) * t1 + R * R + (2 * R * R - 3 * R * r + r * r - 3 * R * R);
            //Debug.Log("zero: " + zero);
            //float V0 = 4 * R *R * H;
            //float V = (4f * H / (3f*(t1+r))) * ((R + t1) * (R + t1) * (R + t1) - (R - r) * (R - r) * (R - r));

            Vector O2, B, I, PI, IF, IB, F, MP, O2M, M;
            for (int i = 0; i < nbOfParticles; i++)
            {
                targetVertices[i] = vertices[i].GetPosition() - meanPoint;
                double xz = Math.Sqrt(targetVertices[i].x * targetVertices[i].x + targetVertices[i].z * targetVertices[i].z);

                //O2 = new Vector(0, targetVertices[i].y, targetVertices[i].z);
                O2 = new Vector(0, targetVertices[i].y);
                Vector O2P = targetVertices[i] - O2;
                B = new Vector(targetVertices[i].x, -innerRadius.y, targetVertices[i].z);
                F = new Vector(targetVertices[i].x, innerRadius.y, targetVertices[i].z);
                I = B + new Vector(0f, t1 / (t1 + r) * H, 0f);
                PI = I - targetVertices[i];
                //Console.WriteLine(i + ", I: " + I + ", PI: " + PI);

                Vector O2Pnorm = (Vector)O2P.Clone(); O2Pnorm.normalize();
                //Console.WriteLine(i + ", " + O2P + ", " + O2Pnorm);
                if (targetVertices[i].y > I.y)
                {
                    IF = F - I;
                    MP = (2*xz * r * PI.norm() / IF.norm()) * O2Pnorm;
                    O2M = O2P - MP;
                }
                else
                {
                    IB = B - I;
                    MP = -(2*xz * t1 * PI.norm() / IB.norm()) * O2Pnorm;
                    O2M = O2P - MP;
                }

                M = O2M + O2;
                targetVertices[i] = M;// +meanPoint;
                //Console.WriteLine(i + ", " +M+", "+ targetVertices[i]);
            }

            if (polarisation == Vector.down)
            {
                for (int i = 0; i < nbOfParticles; i++)
                {
                    targetVertices[i].y = -targetVertices[i].y;
                    vertices[i].v = 2 * meanPoint - vertices[i].v;
                }
            }

            //Console.WriteLine(nuclei[0].v + "," + meanPoint);
            for (int i = 0; i < nbOfParticles; i++)
            {
                targetVertices[i] += meanPoint;
                //Console.WriteLine(i + ", " + targetVertices[i]);
                nucleusEdges[i].l0 = (targetVertices[i] - meanPoint).norm();
                //if(cellId==5 || cellId==73)
                    //Console.WriteLine(cellId + ", " + i + ", " + nucleusEdges[i].l0);
            }

            for (int j = 0; j < edgeCount(); j++)
            {
                edges[j].l0 = Vector.Distance(targetVertices[edges[j].ends[0].pos], targetVertices[edges[j].ends[1].pos]);
            }
        }

        public void ExpandShrinkHexagonalCell(float deltaR)
        {
            Vector meanPoint = new Vector();
            for(int i=0; i<nbOfParticles; i++)
            {
                meanPoint += vertices[i].v;
            }
            meanPoint /= nbOfParticles;

            for(int i=0; i<nbOfParticles; i++)
            {
                double x = vertices[i].v.x - meanPoint.x;
                double y = vertices[i].v.y - meanPoint.y;
                double z = vertices[i].v.z - meanPoint.z;

                double R = Math.Sqrt(x * x + z * z);

                Vector direction = new Vector(x, 0, z);

                direction = (R / innerRadius.x) * (innerRadius.x + deltaR) * direction / direction.norm();

                targetVertices[i] = new Vector(direction.x,y,direction.z) + meanPoint;
            }

            for (int i = 0; i < nbOfParticles; i++)
            {
                nucleusEdges[i].l0 = (targetVertices[i] - meanPoint).norm();
            }

            for (int j = 0; j < edgeCount(); j++)
            {
                edges[j].l0 = Vector.Distance(targetVertices[edges[j].ends[0].pos], targetVertices[edges[j].ends[1].pos]);
            }
        }

        public void Cylinder34Squamous2Columnar(float scaleH)
        {
            //int gridSize = 2;
            //int ySize = 3;
            //innerRadius = new Vector(.5f, .5f, .5f);
            //innerRadius0 = new Vector(r, 1f, r);

            Mesh m = new Cylinder34(new Vector(.5f, .5f, .5f));
            targetVertices = m.GetVertexPositions();
            innerRadius.y *= scaleH;

            for (int i = 0; i < vertices.getCount(); i++)
            {
                //Vector vec = new Vector();
                //vec.x = 2 * innerRadius.x * vec.x / gridSize - innerRadius.x;
                targetVertices[i].y = targetVertices[i].y * scaleH;
                //vec.z = 2 * innerRadius.z * vec.z / gridSize - innerRadius.z;
            }

            for (int i = 0; i < nbOfParticles; i++)
            {
                nucleusEdges[i].l0 = (targetVertices[i]).norm();
            }

            for (int j = 0; j < edgeCount(); j++)
            {
                edges[j].l0 = Vector.Distance(targetVertices[edges[j].ends[0].pos], targetVertices[edges[j].ends[1].pos]);
            }
        }

        public void Cylinder34Squamous2Columnar2(float scaleH)
        {
            //int gridSize = 2;
            //int ySize = 3;
            //innerRadius = new Vector(.5f, .5f, .5f);
            //innerRadius0 = new Vector(r, 1f, r);

            Mesh m = new Cylinder34(new Vector(.5f, .5f, .5f));
            targetVertices = m.GetVertexPositions();
            innerRadius.y *= scaleH;

            Vector nucleus = new Vector(0, -.5f, 0);
            for (int i = 0; i < vertices.getCount(); i++)
            {
                //Vector vec = new Vector();
                //vec.x = 2 * innerRadius.x * vec.x / gridSize - innerRadius.x;
                targetVertices[i].y = targetVertices[i].y * scaleH;
                //vec.z = 2 * innerRadius.z * vec.z / gridSize - innerRadius.z;
            }

            for (int i = 0; i < nbOfParticles; i++)
            {
                nucleusEdges[i].l0 = (targetVertices[i] - nucleus).norm();
            }

            for (int j = 0; j < edgeCount(); j++)
            {
                edges[j].l0 = Vector.Distance(targetVertices[edges[j].ends[0].pos], targetVertices[edges[j].ends[1].pos]);
            }
        }

        public void Default2Elliptic(float a, float b, float c, int appliedForces)
        {
            this.appliedForces = appliedForces;
            for (int i = 0; i < nbOfParticles; i++)
            {
                targetVertices[i].x = a * targetVertices[i].x + (1 - a) * nuclei[0].GetPosition().x;
                targetVertices[i].y = b * targetVertices[i].y + (1 - b) * nuclei[0].GetPosition().y;
                targetVertices[i].z = c * targetVertices[i].z + (1 - c) * nuclei[0].GetPosition().z;
            }
        }

        public void Default2Elliptic2(Vector radius)
        {
            for (int i = 0; i < nbOfParticles; i++)
            {
                vertices[i].v.x = radius.x * vertices[i].v.x + (1 - radius.x) * nuclei[0].GetPosition().x;
                vertices[i].v.y = radius.y * vertices[i].v.y + (1 - radius.y) * nuclei[0].GetPosition().y;
                vertices[i].v.z = radius.z * vertices[i].v.z + (1 - radius.z) * nuclei[0].GetPosition().z;
            }
            innerRadius %= radius;
            SetEdgeELengths();
        }

        private Plane CreateCutPlane(Vector normal)
        {
            GetElongationAxis();
            int min = (int)elongationAxis[0];
            int max = (int)elongationAxis[1];
            //Console.WriteLine(Vector.Distance(vertices[min].v, vertices[max].v));
            Vector centre = Vector.Lerp(vertices[min].GetPosition(), vertices[max].GetPosition(), .5f);

            potentialNucleiPos = new Vector[2];
            potentialNucleiPos[0] = Vector.Lerp(vertices[max].GetPosition(), centre, .5f);
            potentialNucleiPos[1] = Vector.Lerp(vertices[min].GetPosition(), centre, .5f);

            if (normal == Vector.zero)
            {
                normal = (vertices[max].GetPosition() - vertices[min].GetPosition());
                normal.normalize();
            }

            Vector newNormal = (vertices[max].GetPosition() - vertices[min].GetPosition()) ^ (ComputeCentreFromMesh() - vertices[max].GetPosition());

            return new Plane(centre, normal);
            //return new Plane(centre, newNormal);
        }
        
        public void GetElongationAxis()
        {
            ComputeParticleDistances();
            elongationAxis = Helper.DetermineElongationAxis(particlesDistances, nbOfParticles);
        }

        public void ComputeParticleDistances()
        {
            particlesDistances = new float[nbOfParticles, nbOfParticles];
            for (int i = 0; i < nbOfParticles; i++)
            {
                int j = 0;
                while (j++ < i)
                {
                    particlesDistances[i, j] = (float)Vector.Distance(vertices[j].GetPosition(), vertices[i].GetPosition());
                    particlesDistances[j, i] = particlesDistances[i, j];
                }
            }
        }
        
        public void ResetCell()
        {
            //externalEdges.clear();
            nuclei[0].force = new Vector();
            Parallel.For(0, nbOfParticles, i =>
            {
                vertices[i].internalForces = new Vector();
                vertices[i].externalForces = new Vector();
                vertices[i].nucleusForce0 = new Vector();
                vertices[i].force = new Vector();
            });
        }
        
        public Vector GetPosition() {
            //return centre;
            return nuclei[0].GetPosition();
            //return ComputeCentreFromMesh();
        }
        
        public void SetEdgeELengths()
        {
            for(int i=0; i < edgeCount(); i++)
            {
                edges[i].l0 = edges[i].length();
            }
            for (int i = 0; i < nucleusEdges.getCount(); i++)
            {
                nucleusEdges[i].l0 = nucleusEdges[i].length();
            }
            for (int i = 0; i < externalEdges.getCount(); i++)
            {
                externalEdges[i].l0 = MGModel.DCol;
            }
        }

        public void ScaleEdgeELengths(float scale)
        {
            for (int i = 0; i < edgeCount(); i++)
            {
                edges[i].l0 *= scale;
                //Console.WriteLine(edges[i].l0);
            }
            for (int i = 0; i < nucleusEdges.getCount(); i++)
            {
                nucleusEdges[i].l0 *= scale;
                //Console.WriteLine(nucleusEdges[i].l0);
            }
        }

        public float ElasticEnergyMembraneRays()
        {
            double E=0;
            for(int i=0; i<edgeCount(); i++)
            {
                E += (edges[i].l0 - edges[i].length())* (edges[i].l0 - edges[i].length());
            }
            return (float)E;
        }

        public float ElasticEnergyNucleusRays()
        {
            double E = 0;
            for (int i = 0; i < nucleusEdges.getCount(); i++)
            {
                E += (nucleusEdges[i].l0 - nucleusEdges[i].length()) * (nucleusEdges[i].l0 - nucleusEdges[i].length());
            }
            return (float)E;
        }

        public float MorseEnergyMembraneRays()
        {
            double E = 0;
            for (int i = 0; i < edgeCount(); i++)
            {
                double exponent = MGModel.rho * (edges[i].l0 - edges[i].length());
                E += MGModel.J[tissueId, tissueId] * (Math.Exp(2 * exponent) - MGModel.alpha * Math.Exp(exponent));
            }
            return (float)E;
        }

        public float MorseEnergyNucleusRays()
        {
            double E = 0;
            for (int i = 0; i < nucleusEdges.getCount(); i++)
            {
                double exponent = MGModel.rho * (nucleusEdges[i].l0 - nucleusEdges[i].length());
                E += MGModel.J[tissueId, tissueId] * (Math.Exp(2 * exponent) - MGModel.alpha * Math.Exp(exponent));
            }
            return (float)E;
        }

        public void Randomize()
        {
            Randomizer.Randomize<int>(sigma);
        }
    }
}