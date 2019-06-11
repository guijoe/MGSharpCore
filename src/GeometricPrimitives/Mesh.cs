using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using MGSharp.Core.LinearAlgebra;

namespace MGSharp.Core.GeometricPrimitives
{
    public class Mesh
    {
        public const int F0DEL = 1 << 16;
        public const int F1DEL = 2 << 16;
        public const int SINGLEFACE = 3 << 16;
        public const int DELETED = 4 << 16;
        public const int SPLIT = 8 << 16;
        public const int BOUNDARY = 16 << 16;
        public const int LABELMASK = 65535;
        public const int BOUNDARY_RECTANGLE = 1;
        public const int BOUNDARY_CIRCLE = 2;
        public const int BOUNDARY_CHORD_RECTANGLE = 3;
        public const int BOUNDARY_CHORD_CIRCLE = 4;
        public const int BOUNDARY_PRESET = 5;
        public const int GEOMETRY_IMAGE = 1;
        public const int NORMAL_MAP = 2;

        public bool test(int i)
        {
            if (i == 0) return false;
            else return true;
        }

        public void copy(Mesh M)
        {
            int i;
            destroy();
            for (i = 0; i < M.vertices.getCount(); i++)
            {
                M.vertices[i].pos = AddVertex(M.vertices[i].v);
                vertices[M.vertices[i].pos].tex = (Vector)M.vertices[i].tex.Clone();
                vertices[M.vertices[i].pos].normal = (Vector)M.vertices[i].normal.Clone();
                vertices[M.vertices[i].pos].n = (Vector)M.vertices[i].n.Clone();
            }
            for (i = 0; i < M.faces.getCount(); i++)
            {
                AddTriangle(M.faces[i].vertices[0].pos,
                    M.faces[i].vertices[1].pos,
                    M.faces[i].vertices[2].pos);
            }
        }

        public Mesh Clone()
        {
            Mesh m = new Mesh();
            m.copy(this);
            return m;
        }

        public Mesh(Mesh M)
        {
            vertices = new VertexSet();
            edges = new EdgeSet();
            faces = new FaceSet();
            copy(M);
        }

        public Mesh()
        {
            vertices = new VertexSet();
            edges = new EdgeSet();
            faces = new FaceSet();
        }

        public void destroy()
        {
            vertices.clear();
            edges.clear();
            faces.clear();
        }

        //will always return vertices in sequence, unless vertices are deleted
        public int AddVertex(Vector v)
        {
            Vertex vert = new Vertex();
            vert.v = (Vector)v.Clone();
            vert.n = new Vector(0, 1, 0);
            vert.normal = new Vector(0, 1, 0);
            vert.tex = new Vector(0, 0, 0);
            //return vertices.addWithoutCheck(vert);
            return vertices.add(vert);
        }

        public void SetTex(int v, Vector tex)
        {
            if (v < 0) return;
            if (v >= vertices.getCount()) return;
            vertices[v].tex = (Vector)tex.Clone();
        }

        public void SetVertex(int v, Vector vert)
        {
            if (v < 0) return;
            if (v >= vertices.getCount()) return;
            vertices[v].v = (Vector)vert.Clone();
        }

        public void SetNormal(int v, Vector n)
        {
            if (v < 0) return;
            if (v >= vertices.getCount()) return;
            vertices[v].n = (Vector)n.Clone();
            vertices[v].normal = (Vector)n.Clone();
        }

        public Vertex getVertex(int i)
        {
            if (i < 0) return null;
            if (i >= vertices.getCount()) return null;
            return vertices[i];
        }

        public void Recenter()
        {
            int i;
            Vector min, max, c;
            if (vertices.getCount() == 0) return;
            min = (Vector)vertices[0].v.Clone();
            max = (Vector)vertices[0].v.Clone();
            for (i = 1; i < vertices.getCount(); i++)
            {
                if (vertices[i].v.x < min.x)
                    min.x = vertices[i].v.x;
                if (vertices[i].v.y < min.y)
                    min.y = vertices[i].v.y;
                if (vertices[i].v.z < min.z)
                    min.z = vertices[i].v.z;
                if (vertices[i].v.x > max.x)
                    max.x = vertices[i].v.x;
                if (vertices[i].v.y > max.y)
                    max.y = vertices[i].v.y;
                if (vertices[i].v.z > max.z)
                    max.z = vertices[i].v.z;
            }
            c = (max + min) / 2.0;
            for (i = 0; i < vertices.getCount(); i++)
            {
                vertices[i].v -= c;
            }
        }

        public void Rescale()
        {
            int i;
            Vector min, max;
            double maxw;
            if (vertices.getCount() == 0) return;
            min = (Vector)vertices[0].v.Clone();
            max = (Vector)vertices[0].v.Clone();
            for (i = 1; i < vertices.getCount(); i++)
            {
                if (vertices[i].v.x < min.x)
                    min.x = vertices[i].v.x;
                if (vertices[i].v.y < min.y)
                    min.y = vertices[i].v.y;
                if (vertices[i].v.z < min.z)
                    min.z = vertices[i].v.z;
                if (vertices[i].v.x > max.x)
                    max.x = vertices[i].v.x;
                if (vertices[i].v.y > max.y)
                    max.y = vertices[i].v.y;
                if (vertices[i].v.z > max.z)
                    max.z = vertices[i].v.z;
            }
            maxw = max.x - min.x;
            if (max.y - min.y > maxw) maxw = max.y - min.y;
            if (max.z - min.z > maxw) maxw = max.z - min.z;
            for (i = 0; i < vertices.getCount(); i++)
            {
                vertices[i].v /= maxw;
            }
        }

        public int AddTriangle(int v1, int v2, int v3)
        {
            int i;
            Vertex v;
            if ((v1 == v2) || (v1 == v3) || (v2 == v3)) return -1;
            Edge e1 = AddEdge(v1, v2);
            Edge e2 = AddEdge(v2, v3);
            Edge e3 = AddEdge(v3, v1);
            Face f = new Face();
            f.vertices.add(vertices[v1]);
            f.vertices.add(vertices[v2]);
            f.vertices.add(vertices[v3]);

            f.vertices[0].pos = v1;
            f.vertices[1].pos = v2;
            f.vertices[2].pos = v3;

            f.edges.add(e1);
            f.edges.add(e2);
            f.edges.add(e3);
            //check for duplicate triangles
            //if the face exists, all vertices will have the face
            //   in their adjacency lists
            v = vertices[v1];
            for (i = 0; i < v.faces.getCount(); i++)
            {
                if (f.Equals(v.faces[i]))
                {
                    return faces.index(v.faces[i]);
                }
            }
            vertices[v1].add(f);
            vertices[v2].add(f);
            vertices[v3].add(f);
            e1.add(f);
            e2.add(f);
            e3.add(f);
            return faces.addWithoutCheck(f);
        }

        //take a triangle already in the mesh and read it
        public void AddTriangle(int v1, int v2, int v3, Face f)
        {
            int i;
            Vertex v;
            Edge e1 = AddEdge(v1, v2);
            Edge e2 = AddEdge(v2, v3);
            Edge e3 = AddEdge(v3, v1);
            f.vertices.add(vertices[v1]);
            f.vertices.add(vertices[v2]);
            f.vertices.add(vertices[v3]);
            f.edges.add(e1);
            f.edges.add(e2);
            f.edges.add(e3);
            //check for duplicate triangles
            //if the face exists, all vertices will have the face
            //   in their adjacency lists
            v = vertices[v1];
            for (i = 0; i < v.faces.getCount(); i++)
            {
                if (f.Equals(v.faces[i]))
                {
                    return;
                }
            }
            vertices[v1].add(f);
            vertices[v2].add(f);
            vertices[v3].add(f);
            e1.add(f);
            e2.add(f);
            e3.add(f);
        }

        public void ReverseNormals()
        {
            int i;

            for (i = 0; i < vertices.getCount(); i++)
            {
                vertices[i].normal = -vertices[i].normal;
            }
            for (i = 0; i < faces.getCount(); i++)
            {
                faces[i].normal = -faces[i].normal;
            }
        }

        public int Genus()
        {
            int i, j;
            int N, F, E, G, B;
            bool loop = true;
            //Euler-Poincare formula: N+F-E=2-2G-B
            //N=no. vertices
            //F=no. faces
            //E=no. edges
            //G=genus
            //B=no. of boundary loops
            //thus G=(2-B-N-F+E)/2;

            N = vertices.getCount();
            F = faces.getCount();
            E = edges.getCount();
            for (i = 0; i < E; i++)
            {
                edges[i].label = 0;
            }
            B = 0;
            //visit boundaries
            while (loop)
            {
                loop = false;
                for (i = 0; i < E; i++)
                {
                    if ((edges[i].label == 0) && (edges[i].boundary()))
                    {
                        loop = true;
                        B++;
                        Edge e = edges[i];
                        while (e != null)
                        {
                            e.label = 1;
                            Edge e2 = null;
                            for (j = 0; j < e.ends[0].edges.getCount(); j++)
                            {
                                if ((e.ends[0].edges[j].label == 0) && (e.ends[0].edges[j].boundary()))
                                {
                                    e2 = e.ends[0].edges[j];
                                }
                            }
                            if (e2 == null)
                                for (j = 0; j < e.ends[1].edges.getCount(); j++)
                                {
                                    if ((e.ends[1].edges[j].label == 0) && (e.ends[1].edges[j].boundary()))
                                    {
                                        e2 = e.ends[1].edges[j];
                                    }
                                }
                            e = e2;
                        }
                    }
                    else edges[i].label = 1;
                }
            }

            G = (2 - B - N - F + E);
            return G / 2;
        }

        public void ComputeFaceNormals()
        {
            int i;
            for (i = 0; i < faces.getCount(); i++)
            {
                faces[i].ComputeNormal();
            }
        }

        public void ComputeVertexNormals()
        {
            int i;

            ComputeFaceNormals();
            for (i = 0; i < vertices.getCount(); i++)
            {
                vertices[i].ComputeNormal();
            }
        }

        public bool Connected()
        {
            if (!ConnectedFaces()) return false;
            if (!ConnectedEdges()) return false;
            return true;
        }

        public bool ConnectedFaces()
        {
            int i;

            if (faces.getCount() == 0) return true;
            for (i = 0; i < faces.getCount(); i++)
            {
                faces[i].label = 0;
            }

            CheckConnected(faces[0]);
            for (i = 0; i < faces.getCount(); i++)
            {
                if (faces[i].label == 0) return false;
            }
            return true;
        }

        public bool ConnectedEdges()
        {
            int i;

            if (edges.getCount() == 0) return true;

            for (i = 0; i < edges.getCount(); i++)
            {
                edges[i].label = 0;
            }

            for (i = 0; i < vertices.getCount(); i++)
            {
                vertices[i].label = 0;
            }

            CheckConnected(edges[0]);
            for (i = 0; i < edges.getCount(); i++)
            {
                if (edges[i].label == 0) return false;
            }
            for (i = 0; i < vertices.getCount(); i++)
            {
                if (vertices[i].label == 0) return false;
            }
            return true;
        }

        //only 2-manifold meshes can be cut and parameterised.
        public bool Manifold()
        {
            int i;
            for (i = 0; i < edges.getCount(); i++)
            {
                if (edges[i].faces.getCount() > 2) return false;
            }
            return true;
        }

        //Check that no vertex is connected to two or more
        //	triangle fans, if so, split the fans
        //This happens if the mesh touches itself at exactly
        // one point, a vertex (i.e. not an edge or face). 
        public void CheckFans()
        {
            int i, j, k;
            bool recheck = false;
            for (i = 0; i < vertices.getCount(); i++)
            {
                Vertex v = vertices[i];
                for (j = 0; j < v.faces.getCount(); j++)
                {
                    v.faces[j].label = 0;
                }
                if (v.faces.getCount() == 0) continue;
                v.faces[0].label = 1;
                bool changed = true;
                while (changed)
                {
                    changed = false;
                    for (j = 0; j < v.faces.getCount(); j++)
                    {
                        for (k = j + 1; k < v.faces.getCount(); k++)
                        {
                            if ((v.faces[j].label == 0) || (v.faces[k].label == 0))
                            {
                                Edge e = v.faces[j].adjacent(v.faces[k]);
                                if (e != null)
                                {
                                    if (v.faces[j].label > 0)
                                    {
                                        v.faces[k].label = v.faces[j].label;
                                        changed = true;
                                    }
                                    else if (v.faces[k].label > 0)
                                    {
                                        v.faces[j].label = v.faces[k].label;
                                        changed = true;
                                    }
                                }
                            }
                        }
                    }
                }

                bool touch = false;
                for (j = 0; j < v.faces.getCount(); j++)
                {
                    if (v.faces[j].label == 0)
                    {
                        touch = true;
                    }
                }
                if (touch == true)
                {
                    Console.Write("Found touching vertex\n");
                    Vertex vnew = new Vertex();
                    vnew.dist = v.dist;
                    vnew.label = v.label;
                    vnew.n = (Vector)v.n.Clone();
                    vnew.normal = (Vector)v.normal.Clone();
                    vnew.tex = (Vector)v.tex.Clone();
                    vnew.v = (Vector)v.v.Clone();
                    FaceSet vfaces = v.faces.Clone();
                    vertices.addWithoutCheck(vnew);
                    recheck = true;
                    for (j = 0; j < vfaces.getCount(); j++)
                    {
                        if (vfaces[j].label == 1)
                        {
                            vfaces[j].replace(v, vnew);
                            for (k = 0; k < vfaces[j].edges.getCount(); k++)
                            {
                                vfaces[j].edges[k].replace(v, vnew);
                            }
                        }
                    }
                }
            }

            if (recheck) CheckFans();
        }

        //cut mesh to be equivalent to a disc
        //  only works if the mesh does not
        //  "cross" itself. That is, inside and outside is always clear.
        public void DiscCut()
        {
            int i, j;
            EdgeSet cutpath = new EdgeSet();
            EdgeSet candidates = new EdgeSet();
            Edge selected;
            double min;
            int valence, x;
            Vertex v;
            bool changed;
            bool boundary = false;

            //mark boundary edges
            for (i = 0; i < edges.getCount(); i++)
            {
                edges[i].label = 0;
                //mark boundaries, marked as split as well, so that we don't try to split them
                if (edges[i].faces.getCount() < 2)
                {
                    edges[i].label |= F1DEL | BOUNDARY | SPLIT;
                    boundary = true;
                }
            }

            //no faces are deleted
            for (i = 0; i < faces.getCount(); i++)
            {
                faces[i].label = 0;
            }

            //for computation of geodesic distance
            for (i = 0; i < vertices.getCount(); i++)
            {
                vertices[i].dist = 1e15;
            }
            //remove a seed triangle
            int fs = 0;
            DiscDeleteFace(faces[fs], candidates);
            //update geodesic distance for each 
            // of the seed triangle vertices
            VertexSet border = new VertexSet();
            border.add(faces[fs].vertices[0]);
            border.add(faces[fs].vertices[1]);
            border.add(faces[fs].vertices[2]);
            Console.Write("Distance...\n");
            CalculateGeodesicDistance(border);

            //successively remove edges that are only adjacent 
            // to one triangle, and remove the edge and triangle
            Console.Write("Removing triangles...\n");
            selected = candidates[0];
            while (selected != null)
            {
                selected = null;
                min = 1e16;
                min = -1.0;
                for (i = 0; i < candidates.getCount(); i++)
                {
                    if ((!(test(candidates[i].label & BOUNDARY))) &&
                        ((candidates[i].label & SINGLEFACE) > 0))
                    {
                        if (!((test(candidates[i].label & F0DEL)) && (test(candidates[i].label & F1DEL))))
                        {
                            if (candidates[i].ends[0].dist > min)
                            {
                                selected = candidates[i];
                                min = selected.ends[0].dist;
                            }
                            if (candidates[i].ends[1].dist > min)
                            {
                                selected = candidates[i];
                                min = selected.ends[1].dist;
                            }
                        }
                    }
                }
                if (selected != null)
                {
                    //one of the two must be deleted,
                    // one is deleted already
                    //  DiscDeleteFace checks if the face is already deleted
                    // The corresponding edge is also deleted.
                    candidates.remove(selected);
                    DiscDeleteEdge(selected);
                    DiscDeleteFace(selected.faces[0], candidates);
                    DiscDeleteFace(selected.faces[1], candidates);
                }
            }

            //the cut path can now be built
            for (i = 0; i < edges.getCount(); i++)
            {
                if ((!(test(edges[i].label & DELETED))) &&
                    (!(test(edges[i].label & BOUNDARY))))
                    cutpath.add(edges[i]);
            }

            Console.Write("Pruning...{0}\n", cutpath.getCount());
            //prune vertices that
            // are connected to only one edge, delete the edge also.
            // Repeat until only loops remain
            changed = true;
            while (changed)
            {
                changed = false;
                for (i = cutpath.getCount() - 1; i >= 0; i--)
                {
                    valence = 0;
                    v = cutpath[i].ends[0];
                    for (j = 0; j < v.edges.getCount(); j++)
                    {
                        if (!(test(v.edges[j].label & DELETED))) valence++;
                    }
                    //store the valence in the label
                    v.label = valence;
                    //if the valence is one, prune the edge
                    if (valence == 1)
                    {
                        cutpath[i].label |= DELETED;
                        cutpath.remove(i);
                        changed = true;
                        continue;
                    }
                    //Recompute valence by disregarding boundary edges, since they are not part of the cut
                    //  homever, if a boundary is found, then the valence should be increased by one, since
                    //  a boundary can be regarded as a cut that has already been made.
                    //  A vertex can be joined to a maximum of two boundary edges, any more and the graph is
                    //	not connected nicely. Should this be detected and cut elsewhere? 
                    valence = 0;
                    x = 0;
                    for (j = 0; j < v.edges.getCount(); j++)
                    {
                        if (test(v.edges[j].label & BOUNDARY)) x = 1;
                        if ((!(test(v.edges[j].label & DELETED))) && (!(test(v.edges[j].label & BOUNDARY)))) valence++;
                    }
                    //store the valence in the label
                    v.label = valence + x;

                    valence = 0;
                    v = cutpath[i].ends[1];
                    for (j = 0; j < v.edges.getCount(); j++)
                    {
                        if (!(test(v.edges[j].label & DELETED))) valence++;
                    }
                    //store the valence in the label
                    v.label = valence;
                    //if the valence is one, prune the edge
                    if (valence == 1)
                    {
                        cutpath[i].label |= DELETED;
                        cutpath.remove(i);
                        changed = true;
                        continue;
                    }
                    //recompute valence by disregarding boundary edges, since they are not part of the cut
                    valence = 0;
                    x = 0;
                    for (j = 0; j < v.edges.getCount(); j++)
                    {
                        if (test(v.edges[j].label & BOUNDARY)) x = 1;
                        if ((!(test(v.edges[j].label & DELETED))) && (!(test(v.edges[j].label & BOUNDARY)))) valence++;
                    }
                    //store the valence in the label
                    v.label = valence + x;
                }
            }

            //now try to remove serrations in the path
            // we do this by checxkxing the triangles adjacent to
            // each edge in the path, if another edge on the
            // path is part of the triangle, then we check 
            // the third edge. If the length of the third edge
            // is shorter than the first two, then we replace 
            // the two edges by the third edge.
            Console.Write("Smoothing cut...\n");
            while (SmoothCut(cutpath)) ;
            Random rnd = new Random();

            //Test for genus 0 surface
            if ((cutpath.getCount() == 0) && (!boundary))
            {
                //If we get a genus 0 surface, we should rather use
                // an alternative parameterisation and avoid the
                //	disc cut

                //create a random connected path to cut along	
                VertexSet cutvert = new VertexSet();
                i = rnd.Next() % edges.getCount();
                cutpath.add(edges[i]);
                min = edges.getCount();
                //min=sqrt(min)*2;
                min = 8;
                v = edges[i].ends[1];
                cutvert.add(edges[i].ends[0]);
                cutvert.add(edges[i].ends[1]);
                //remark edges
                for (i = 0; i < edges.getCount(); i++)
                {
                    edges[i].label = DELETED;
                }
                for (i = 0; i < min; i++)
                {
                    selected = v.edges[rnd.Next() % (v.edges.getCount())];
                    while (cutpath.contains(selected))
                    {
                        selected = v.edges[rnd.Next() % (v.edges.getCount())];
                    }
                    Vertex nv;
                    if (selected.ends[0] == v) nv = selected.ends[1];
                    else nv = selected.ends[0];
                    if (!cutvert.contains(nv))
                    {
                        cutvert.add(nv);
                        cutpath.add(selected);
                        v = nv;
                    }
                }
                //remark edges
                for (i = 0; i < cutpath.getCount(); i++)
                {
                    cutpath[i].label = 0;
                }
                //compute valences for path
                for (i = 0; i < cutvert.getCount(); i++)
                {
                    valence = 0;
                    v = cutvert[i];
                    for (j = 0; j < v.edges.getCount(); j++)
                    {
                        if (!(test(v.edges[j].label & DELETED))) valence++;
                    }
                    //store the valence in the label
                    v.label = valence;
                }
                //while (SmoothCut(cutpath));
            }
            if ((cutpath.getCount() == 0) && (!boundary))
            {
                //failed to cut the surface!
                //we should hopefully never get here
            }
            Console.Write("Performing cut...{0}\n", cutpath.getCount());
            CutMesh(cutpath);
            Console.Write("Cut done.\n");
        }

        //Cut for geometry image using the technique described in Gu's thesis 
        //  on geometry images
        //	C# doesn't support default parameters, so two 
        //  methods are required
        public void GuCut() { GuCut(false); }
        public void GuCut(bool flipedge)
        {
            EdgeSet eboundary = new EdgeSet();
            VertexSet vboundary = new VertexSet();
            Edge e, ce;
            Vertex v, vp;
            Face f;
            int i, j, c;
            double maxstretch, stretch, Eold, Enew;
            EdgeSet cutpath = new EdgeSet();
            Mesh opt = null;
            int iter = 1000;

            //first cut into a disc
            Console.Write("Cutting\n");
            DiscCut();

            //then compute Floater parameterisation
            Console.Write("Floater param to circle\n");
            ParamFloater(BOUNDARY_CIRCLE);
            Enew = AllStretch();
            Eold = Enew + 1.0;
            Console.Write("Stretch before: {0}\n", Enew);
            //refine until stretch increases
            while ((Eold > Enew) && (--iter > 0))
            {
                //save curent mesh
                opt = this.Clone();
                if (!flipedge)
                {
                    //If any border triangle is degenerate in parameter space,
                    // that is, three vertices are on the border, split all edges not on the border.
                    // This creates  4 new triangles (usually), and 2 are discarded 
                    cutpath.clear();
                    for (i = 0; i < faces.getCount(); i++)
                    {
                        c = 0;
                        for (j = 0; j < 3; j++)
                            if (faces[i].vertices[j].boundary())
                                c++;
                        if (c == 3)
                        {
                            for (j = 0; j < 3; j++)
                                if (!(faces[i].edges[j].boundary()))
                                {
                                    cutpath.add(faces[i].edges[j]);
                                }
                        }
                    }
                    if (cutpath.getCount() > 0)
                    {
                        for (i = 0; i < cutpath.getCount(); i++)
                        {
                            SplitEdge(cutpath[i]);
                        }
                        //reparameterise
                        Console.Write("Parameterise...\n");
                        ParamFloater(BOUNDARY_CIRCLE);
                        //check if stretch is minimised
                        Eold = Enew;
                        Enew = AllStretch();
                        Eold = Enew + 1.0;
                        Console.Write("Stretch new: {0}\n", Enew);
                        continue;
                    }
                }
                else
                {
                    //flip the offending edge instead of introducing new vertices
                    for (i = 0; i < vertices.getCount(); i++)
                    {
                        vertices[i].pos = i;
                    }
                    for (i = 0; i < faces.getCount(); i++)
                    {
                        c = 0;
                        for (j = 0; j < 3; j++)
                            if (faces[i].vertices[j].boundary())
                                c++;
                        if (c == 3)
                        {
                            Edge et = null;
                            for (j = 0; j < 3; j++)
                            {
                                if (!faces[i].edges[j].boundary())
                                {
                                    et = faces[i].edges[j];
                                }
                            }
                            if (et != null)
                            {
                                Vertex v1, v2, v3 = null, v4 = null;
                                v1 = et.ends[0];
                                v2 = et.ends[1];
                                for (j = 0; j < 3; j++)
                                {
                                    if ((et.faces[0].vertices[j] != v1) &&
                                        (et.faces[0].vertices[j] != v2))
                                        v3 = et.faces[0].vertices[j];
                                    if ((et.faces[1].vertices[j] != v1) &&
                                        (et.faces[1].vertices[j] != v2))
                                        v4 = et.faces[1].vertices[j];
                                }
                                RemoveFace(et.faces[0]);
                                RemoveFace(et.faces[1]);
                                edges.remove(et);
                                AddTriangle(v1.pos, v3.pos, v4.pos);
                                AddTriangle(v2.pos, v3.pos, v4.pos);
                            }
                        }
                    }

                    //reparameterise
                    Console.Write("Parameterise...\n");
                    ParamFloater(BOUNDARY_CIRCLE);
                    //check if stretch is minimised
                    Eold = Enew;
                    Enew = AllStretch();
                    Eold = Enew + 1.0;
                    Console.Write("Stretch new: {0}\n", Enew);
                    continue;
                }
                maxstretch = -1;
                f = null;
                for (i = 0; i < faces.getCount(); i++)
                {
                    //don't consider triangles on the boundary
                    if (faces[i].edges[0].boundary()) continue;
                    if (faces[i].edges[1].boundary()) continue;
                    if (faces[i].edges[2].boundary()) continue;
                    stretch = faces[i].stretch();
                    if (maxstretch < stretch)
                    {
                        maxstretch = stretch;
                        f = faces[i];
                    }
                }
                if (f == null) Console.Write("No maximal stretch!\n");
                //Check for high stretch.
                //if found, join vertex of high stretch to current border 
                //  using shortest path	
                eboundary.clear();
                vboundary.clear();
                ComputeBoundary(eboundary, vboundary);

                Console.Write("Calculating distances...\n");
                CalculateGeodesicDistance(vboundary);

                v = f.vertices[0];
                cutpath.clear();
                //follow path from vertex to boundary
                Console.Write("Building cutpath...\n");
                while (v.next != null)
                {
                    vp = v.next;
                    ce = new Edge();
                    ce.ends[0] = v;
                    ce.ends[1] = vp;
                    //find matching edge
                    e = v.edges.match(ce);
                    if (e == null) Console.Write("Bad edge! v1={0} v2={1}\n", v, vp);
                    //add to cutpath
                    cutpath.add(e);
                    if (e.boundary()) Console.Write("Boundary edge!\n");
                    //advance to next node
                    v = vp;
                }

                Console.Write("Cutting...\n");
                //cut the mesh
                CutMesh(cutpath);

                //reparameterise
                Console.Write("Parameterise...\n");
                ParamFloater(BOUNDARY_CIRCLE);
                //check if stretch is minimised
                Eold = Enew;
                Enew = AllStretch();
                Console.Write("Stretch new: {0}\n", Enew);
            }
            //revert to previous mesh (stretch was better)
            copy(opt);
            Enew = AllStretch();
            Console.Write("Stretch after: {0}\n", Enew);
        }

        public void SplitEdge(Edge e)
        {
            Face f;
            Vector pos;
            int v1, v2, v3, vpos, i, j;

            pos = (e.ends[0].v + e.ends[1].v) / 2.0;
            vpos = AddVertex(pos);

            //remove edge
            e.ends[0].edges.remove(e);
            e.ends[1].edges.remove(e);
            edges.remove(e);

            for (i = 0; i < 2; i++)
            {
                f = e.faces[i];
                if (f != null)
                {
                    for (j = 0; j < 3; j++)
                    {
                        f.vertices[j].faces.remove(f);
                        f.edges[j].remove(f);
                    }
                    for (j = 0; j < 3; j++)
                    {
                        if ((f.vertices[j] != e.ends[0]) && (f.vertices[j] != e.ends[1]))
                        {
                            Vertex tmp = f.vertices[0];
                            f.vertices[0] = f.vertices[j];
                            f.vertices[j] = tmp;
                        }
                    }
                    v1 = vertices.index(f.vertices[0]);
                    v2 = vertices.index(f.vertices[1]);
                    v3 = vertices.index(f.vertices[2]);

                    faces.remove(f);
                    AddTriangle(v1, v2, vpos);
                    AddTriangle(v1, vpos, v3);
                }
            }

        }

        public void CutMesh(EdgeSet cutpath)
        {
            int i, j, k, valence, x;
            VertexSet cutvert = new VertexSet();
            Vertex[] newvert;
            int count;
            Vertex v;
            Edge adj;
            bool changed;
            Edge e;

            Console.Write("Cutting {0}/{1} edges\n", cutpath.getCount(), edges.getCount());
            //label edges as if they are not in the cutpath
            for (i = 0; i < edges.getCount(); i++)
            {
                if (edges[i].boundary()) edges[i].label = BOUNDARY | SPLIT;
                else edges[i].label = DELETED;
            }

            //label edges in cutpath as not deleted
            for (i = 0; i < cutpath.getCount(); i++)
            {
                //label edge as not deleted
                if (!cutpath[i].boundary())
                    cutpath[i].label = 0;
                else
                {
                    Console.Write("Boundary in cutpath!\n");
                    cutpath.remove(i);
                    i--;
                }
            }

            //label vertices according to valence (undeleted edges)
            for (i = 0; i < cutpath.getCount(); i++)
            {
                valence = 0;
                v = cutpath[i].ends[0];
                x = 0;
                for (j = 0; j < v.edges.getCount(); j++)
                {
                    if (test(v.edges[j].label & BOUNDARY)) x = 1;
                    if ((!(test(v.edges[j].label & DELETED))) && (!(test(v.edges[j].label & BOUNDARY)))) valence++;
                }
                //store the valence in the label
                v.label = valence + x;

                valence = 0;
                v = cutpath[i].ends[1];
                x = 0;
                for (j = 0; j < v.edges.getCount(); j++)
                {
                    if (test(v.edges[j].label & BOUNDARY)) x = 1;
                    if ((!(test(v.edges[j].label & DELETED))) && (!(test(v.edges[j].label & BOUNDARY)))) valence++;
                }
                //store the valence in the label
                v.label = valence + x;
            }

            //find the vertices in the cut
            for (i = 0; i < cutpath.getCount(); i++)
            {
                //check valence
                //for genus 0 surfaces, the cut will not necessarily be a loop
                if (cutpath[i].ends[0].label >= 2)
                    cutvert.add(cutpath[i].ends[0]);
                if (cutpath[i].ends[1].label >= 2)
                    cutvert.add(cutpath[i].ends[1]);
            }

            for (i = 0; i < cutvert.getCount(); i++)
            {
                count = 0;
                v = cutvert[i];
                //create a new vertex for each edge
                // entering this vertex (each edge on a cutpath)
                newvert = new Vertex[v.label];
                newvert[0] = v;
                for (j = 1; j < v.label; j++)
                {
                    newvert[j] = new Vertex();
                    newvert[j].v = (Vector)v.v.Clone();
                    newvert[j].tex = (Vector)v.tex.Clone();
                    newvert[j].normal = (Vector)v.normal.Clone();
                    newvert[j].n = (Vector)v.n.Clone();
                    newvert[j].dist = v.dist;
                    vertices.add(newvert[j]);
                }
                //reset labels for cut
                for (j = 0; j < v.faces.getCount(); j++)
                {
                    v.faces[j].label = -1;
                }

                changed = true;
                //label the faces
                while (changed)
                {
                    changed = false;
                    //first label a triangle next to a cut edge
                    for (j = 0; (!changed) && (j < v.faces.getCount()); j++)
                    {
                        if (v.faces[j].label < 0)
                            for (k = 0; (!changed) && (k < 3); k++)
                            {
                                if ((!(test(v.faces[j].edges[k].label & DELETED)))
                                    && (!(test(v.faces[j].edges[k].label & BOUNDARY)))
                                //	&&(!(v.faces[j].edges[k].label&SPLIT)) 
                                )
                                {
                                    v.faces[j].label = count;
                                    changed = true;
                                    count++;
                                }
                            }
                    }

                    if (changed == false) break;
                    //now propogate the information to the surrounding faces
                    while (changed)
                    {
                        changed = false;
                        for (j = 0; j < v.faces.getCount(); j++)
                        {
                            for (k = j + 1; k < v.faces.getCount(); k++)
                            {
                                if ((v.faces[j].label < 0) || (v.faces[k].label < 0))
                                {
                                    adj = v.faces[j].adjacent(v.faces[k]);
                                    if (adj != null)
                                    {
                                        if (test(adj.label & DELETED))
                                        {
                                            //faces should have same label
                                            if (v.faces[j].label >= 0)
                                            {
                                                v.faces[k].label = v.faces[j].label;
                                                changed = true;
                                            }
                                            else if (v.faces[k].label >= 0)
                                            {
                                                v.faces[j].label = v.faces[k].label;
                                                changed = true;
                                            }
                                        }
                                        else
                                        {
                                            //faces should have different labels
                                            //not sure what yet though
                                        }
                                    }
                                }
                            }
                        }
                    }
                    //go to next iteration, to choose a new cut face to label
                    changed = true;
                }

                for (j = 0; j < v.faces.getCount(); j++)
                {
                    if (v.faces[j].label < 0)
                    {
                        Console.Write("Critical error, unlabeled face\n");
                    }
                    if (v.faces[j].label >= v.label)
                    {
                        Console.Write("Critical error, incorrectly labeled face v.label={0} face={1}\n", v.label, v.faces[j].label);
                    }
                }

                for (j = 0; j < v.edges.getCount(); j++)
                {
                    if (!(test(v.edges[j].label & DELETED)))
                    {
                        if (v.edges[j].faces.getCount() >= 2)
                        {
                            if (v.edges[j].faces[0].label == v.edges[j].faces[1].label)
                            {
                                Console.Write("Critical error, faces have same label but shouldn't\n");
                                Console.Write("label={0}\n", v.edges[j].faces[0].label);
                            }
                        }
                    }
                    else
                    {
                        if (v.edges[j].faces.getCount() >= 2)
                        {
                            if (v.edges[j].faces[0].label != v.edges[j].faces[1].label)
                            {
                                Console.Write("Critical error, faces should have the same label\n");
                            }
                        }
                    }
                }

                //now, update the edges to use the correct vertex
                //if the edge was a cut edge, split the edge and use the
                //   face labels to determine which edge belongs
                //   to which face. These new edges are no longer cut edges.

                //Console.Write("Splitting edges\n");
                // go backwards, because we will be removing edges while worxkxing	
                for (j = v.edges.getCount() - 1; j >= 0; j--)
                {
                    if ((test(v.edges[j].label & DELETED)) || (test(v.edges[j].label & SPLIT)))
                    {
                        //edge not part of cutpath: just change the vertex
                        //Console.Write("Splitting non-cutpath edge\n");	
                        e = v.edges[j];
                        if (e.faces[0] != null)
                        {
                            if (e.ends[0] == v)
                            {
                                v.edges.remove(e);
                                e.ends[0] = newvert[e.faces[0].label];
                                e.ends[0].edges.add(e);
                            }
                            else
                            if (e.ends[1] == v)
                            {
                                v.edges.remove(e);
                                e.ends[1] = newvert[e.faces[0].label];
                                e.ends[1].edges.add(e);
                            }
                        }
                        else
                        {
                            if (e.ends[0] == v)
                            {
                                v.edges.remove(e);
                                e.ends[0] = newvert[e.faces[1].label];
                                e.ends[0].edges.add(e);
                            }
                            else
                            if (e.ends[1] == v)
                            {
                                v.edges.remove(e);
                                e.ends[1] = newvert[e.faces[1].label];
                                e.ends[1].edges.add(e);
                            }
                        }
                        //Console.Write("Split done\n");	
                    }
                    else
                    {
                        //change the existing edge to use the correct vertex
                        //Console.Write("Splitting cutpath edge face0=%x face1=%x\n", v.edges[j].faces[0], v.edges[j].faces[1]);	
                        Edge newedge;
                        newedge = v.edges[j].Clone();
                        e = v.edges[j];
                        if (e.ends[0] == v)
                        {
                            v.edges.remove(e);
                            e.ends[0] = newvert[e.faces[0].label];
                            e.ends[0].edges.add(e);
                        }
                        else
                        if (e.ends[1] == v)
                        {
                            v.edges.remove(e);
                            e.ends[1] = newvert[e.faces[0].label];
                            e.ends[1].edges.add(e);
                        }
                        //Console.Write("Actual split\n");	
                        //edge part of cutpath: split edge
                        //Console.Write("Replacing ends\n");
                        if (newedge.ends[0] == v)
                        {
                            newedge.ends[0] = newvert[newedge.faces[1].label];
                        }
                        else
                        if (newedge.ends[1] == v)
                        {
                            newedge.ends[1] = newvert[newedge.faces[1].label];
                        }
                        //Console.Write("Adding edge to vertices\n");	
                        newedge.ends[0].add(newedge);
                        newedge.ends[1].add(newedge);
                        //Console.Write("Selecting faces\n");	
                        //newedge.faces[0]=null;
                        newedge.faces.remove(0);
                        //don't have to set the face to null, since replace
                        // will do it for us
                        //e.faces[1]=null;
                        //Console.Write("Replacing edges in triangle %x\n", newedge.faces[1]);	
                        newedge.faces[0].replace(e, newedge);
                        newedge.label |= SPLIT;
                        e.label |= SPLIT;
                        //Console.Write("Adding edge to set of edges\n");	
                        edges.add(newedge);
                        //Console.Write("Split done\n");	
                    }
                }

                //now update the triangles
                for (j = v.faces.getCount() - 1; j >= 0; j--)
                {
                    v.faces[j].replace(v, newvert[v.faces[j].label]);
                }
            }
        }

        public bool SmoothCut(EdgeSet cutpath)
        {
            int i, j, k;
            Edge other;
            Face tri;
            Face f = new Face();
            int[] tag = new int[3] { 0, 0, 0 };
            Edge e1;
            Edge e2;
            Vertex v;
            int valence, x;

            for (i = 0; i < cutpath.getCount(); i++)
            {
                //check for later connecting edges
                for (j = i + 1; j < cutpath.getCount(); j++)
                {
                    //test if they share an endpoint, but there must be only two edges 
                    // sharing this vertex, thus the valence (stored in the label) must be 2
                    tri = null;
                    if (cutpath[i].ends[0].label == 2)
                    {
                        if (cutpath[j].ends[0] == cutpath[i].ends[0])
                        {
                            f.vertices[0] = cutpath[i].ends[0];
                            f.vertices[1] = cutpath[i].ends[1];
                            f.vertices[2] = cutpath[j].ends[1];
                            tri = faces.match(f);
                        }
                        if (cutpath[j].ends[1] == cutpath[i].ends[0])
                        {
                            f.vertices[0] = cutpath[i].ends[0];
                            f.vertices[1] = cutpath[i].ends[1];
                            f.vertices[2] = cutpath[j].ends[0];
                            tri = faces.match(f);
                        }
                    }
                    if (cutpath[i].ends[1].label == 2)
                    {
                        if (cutpath[j].ends[0] == cutpath[i].ends[1])
                        {
                            f.vertices[0] = cutpath[i].ends[0];
                            f.vertices[1] = cutpath[i].ends[1];
                            f.vertices[2] = cutpath[j].ends[1];
                            tri = faces.match(f);
                        }
                        if (cutpath[j].ends[1] == cutpath[i].ends[1])
                        {
                            f.vertices[0] = cutpath[i].ends[0];
                            f.vertices[1] = cutpath[i].ends[1];
                            f.vertices[2] = cutpath[j].ends[0];
                            tri = faces.match(f);
                        }
                    }
                    //if we find such a triangle, find the other edge
                    //  and see if it is shorter
                    if (tri != null)
                    {
                        for (k = 0; k < 3; k++) tag[k] = 0;
                        for (k = 0; k < 3; k++)
                        {
                            if (tri.edges[k] == cutpath[j]) tag[k] = 1;
                            if (tri.edges[k] == cutpath[i]) tag[k] = 1;
                        }
                        other = null;
                        for (k = 0; k < 3; k++)
                        {
                            if (tag[k] == 0) other = tri.edges[k];
                        }
                        if (other.length() < cutpath[i].length() + cutpath[j].length())
                        {
                            //the edge is shorter, replace the other two by this edge
                            e1 = cutpath[i];
                            e2 = cutpath[j];
                            cutpath.add(other);
                            cutpath.remove(i);
                            cutpath.remove(j);
                            e1.label |= DELETED;
                            e2.label |= DELETED;
                            other.label &= ~DELETED;

                            //Recompute valence by disregarding boundary edges, since they are not part of the cut
                            //  homever, if a boundary is found, then the valence should be increased by one, since
                            //  a boundary can be regarded as a cut that has already been made.
                            //  A vertex can be joined to a maximum of two boundary edges, any more and the graph is
                            //	not connected nicely. Should this be detected and cut elsewhere? 
                            v = other.ends[0];
                            valence = 0;
                            x = 0;
                            for (j = 0; j < v.edges.getCount(); j++)
                            {
                                if (test(v.edges[j].label & BOUNDARY)) x = 1;
                                if ((!(test(v.edges[j].label & DELETED))) && (!(test(v.edges[j].label & BOUNDARY)))) valence++;
                            }
                            //store the valence in the label
                            v.label = valence + x;

                            v = other.ends[1];
                            valence = 0;
                            x = 0;
                            for (j = 0; j < v.edges.getCount(); j++)
                            {
                                if (test(v.edges[j].label & BOUNDARY)) x = 1;
                                if ((!(test(v.edges[j].label & DELETED))) && (!(test(v.edges[j].label & BOUNDARY)))) valence++;
                            }
                            //store the valence in the label
                            v.label = valence + x;

                            return true;
                        }
                    }
                }
            }
            return false;
        }

        //calculate the geodesic distance, the shortest path
        // on the surface, between this vertex and all other vertices
        //  uses Dxixxjxkstra's algorithm
        // The labels of vertices are set to point to the next in the path to v
        public void CalculateGeodesicDistance(Vertex v)
        {
            VertexSet set = new VertexSet();
            set.add(v);
            Dxixxjxkstra(set);
        }

        public void CalculateGeodesicDistance(VertexSet set)
        {
            Dxixxjxkstra(set);
        }

        public void Dxixxjxkstra(VertexSet known)
        {
            int i, j;
            VertexSet S = new VertexSet();
            VertexSet W = new VertexSet();
            Vertex x, v;
            double min, d;

            for (j = 0; j < vertices.getCount(); j++)
            {
                vertices[j].label = -1;
                vertices[j].dist = 1e15;
                vertices[j].next = null;
            }
            //apply Dxixxjxkstra's algorithm
            for (j = 0; j < known.getCount(); j++)
            {
                v = known[j];
                v.dist = 0;
                v.label = 0;
            }
            for (j = 0; j < known.getCount(); j++)
            {
                v = known[j];
                //set distance for neighbouring edges
                for (i = 0; i < v.edges.getCount(); i++)
                {
                    d = v.edges[i].length();
                    if (v.edges[i].ends[0] == v)
                    {
                        if (v.edges[i].ends[1].dist > d)
                        {
                            v.edges[i].ends[1].dist = d;
                            v.edges[i].ends[1].next = v;
                        }
                    }
                    else
                    {
                        if (v.edges[i].ends[0].dist > d)
                        {
                            v.edges[i].ends[0].dist = d;
                            v.edges[i].ends[0].next = v;
                        }
                    }
                }

                S.add(v);

                //compute W, the set of all vertices connected to S, but not in S
                for (i = 0; i < v.edges.getCount(); i++)
                {
                    if (v.edges[i].ends[0] == v)
                    {
                        if (!S.contains(v.edges[i].ends[1]))
                            W.add(v.edges[i].ends[1]);
                    }
                    else
                    {
                        if (!S.contains(v.edges[i].ends[0]))
                            W.add(v.edges[i].ends[0]);
                    }
                }
            }

            while ((S.getCount() != vertices.getCount()) && (W.getCount() > 0))
            {
                //find the vertex X in W for which the distance is minimum
                //there will be at least one element in W, because we stop when
                //   W is empty (and S contains all vertices);
                x = W[0];
                min = x.dist;
                for (i = 1; i < W.getCount(); i++)
                {
                    if (W[i].dist < min)
                    {
                        x = W[i];
                        min = x.dist;
                    }
                }
                //Add X to S
                S.add(x);
                //update W, the set of all vertices connected to S, but not in S
                W.remove(x);
                //update closest distance for all vertices outside S (but only if they are connected to X), and update W
                for (i = 0; i < x.edges.getCount(); i++)
                {
                    if (!(S.contains(x.edges[i].ends[0])))
                    {
                        d = x.edges[i].length();
                        if (d + x.dist < x.edges[i].ends[0].dist)
                        {
                            x.edges[i].ends[0].dist = d + x.dist;
                            //label prior vertex in sequence
                            x.edges[i].ends[0].next = x;
                        }
                        W.add(x.edges[i].ends[0]);
                    }
                    if (!(S.contains(x.edges[i].ends[1])))
                    {
                        d = x.edges[i].length();
                        if (d + x.dist < x.edges[i].ends[1].dist)
                        {
                            x.edges[i].ends[1].dist = d + x.dist;
                            //label prior vertex in sequence
                            x.edges[i].ends[1].next = x;
                        }
                        W.add(x.edges[i].ends[1]);
                    }
                }
                //next iteration
            }
        }

        public void DiscDeleteEdge(Edge e)
        {
            if (test(e.label & BOUNDARY)) return;
            e.label |= DELETED;
        }

        public void DiscDeleteFace(Face f, EdgeSet candidates)
        {
            int i;

            if (test(f.label & DELETED)) return;
            f.label |= DELETED;
            for (i = 0; i < 3; i++)
            {
                if (f.edges[i].faces[0] == f)
                {
                    f.edges[i].label |= F0DEL;
                }
                else
                {
                    f.edges[i].label |= F1DEL;
                }
                if ((test(f.edges[i].label & F0DEL)) && (test(f.edges[i].label & F1DEL)))
                {
                    //DiscDeleteEdge(f.edges[i]);
                }
                else
                {
                    candidates.add(f.edges[i]);
                }
            }
        }

        public Edge AddEdge(int v1, int v2)
        {
            int i;
            Vertex v;
            Edge e = new Edge();
            e.ends[0] = vertices[v1];
            e.ends[1] = vertices[v2];
            //if the edge exists, both vertices will have the edge
            //   in their adjacency lists
            v = vertices[v1];
            for (i = 0; i < v.edges.getCount(); i++)
            {
                if (e.Equals(v.edges[i]))
                {
                    return v.edges[i];
                }
            }
            edges.addWithoutCheck(e);
            vertices[v1].add(e);
            vertices[v2].add(e);
            return e;
        }

        public VertexSet ComputeRing(Vertex v, int level)
        {
            VertexSet ring = new VertexSet();
            ring.add(v);
            if (level <= 0) return ring;
            return ComputeRing(ring, level);
        }

        public VertexSet ComputeRing(VertexSet set, int level)
        {
            VertexSet ring = new VertexSet();
            Vertex v;
            int i, j;

            for (i = 0; i < set.getCount(); i++)
            {
                v = set[i];
                ring.add(v);
            }

            if (level <= 0) return ring;

            for (i = 0; i < set.getCount(); i++)
            {
                v = set[i];
                for (j = 0; j < v.edges.getCount(); j++)
                {
                    ring.add(v.edges[j].ends[1]);
                    ring.add(v.edges[j].ends[0]);
                }
            }
            if (level == 1) return ring;
            return ComputeRing(ring, level - 1);
        }

        public VertexSet ComputeExclusiveRing(Vertex v, int level)
        {
            VertexSet ring = new VertexSet();
            if (level <= 0) return ring;
            ring.add(v);
            return ComputeExclusiveRing(ring, level);
        }

        public VertexSet ComputeExclusiveRing(Face t, int level)
        {
            VertexSet ring = new VertexSet();
            if (level <= 0) return ring;
            ring.add(t.vertices[0]);
            ring.add(t.vertices[1]);
            ring.add(t.vertices[2]);
            return ComputeExclusiveRing(ring, level);
        }

        public VertexSet ComputeExclusiveRing(VertexSet set, int level)
        {
            VertexSet ring, inner;
            Vertex v;
            int i;

            inner = ComputeRing(set, level - 1);
            ring = ComputeRing(set, level);
            for (i = 0; i < inner.getCount(); i++)
            {
                v = inner[i];
                ring.remove(v);
            }
            return ring;
        }

        public void CircleBoundary(VertexSet vboundary)
        {
            int i, n;
            Vector p = new Vector();
            //Circle boundary
            n = vboundary.getCount();
            for (i = 0; i < n; i++)
            {
                p.z = 0;
                p.x = Math.Cos(i * 2.0 * Math.PI / n) * 0.5 + 0.5;
                p.y = Math.Sin(i * 2.0 * Math.PI / n) * 0.5 + 0.5;
                vboundary[i].tex = (Vector)p.Clone();
            }
        }

        public void CircleChordBoundary(VertexSet vboundary)
        {
            int i, n;
            double dist, cd;
            Vector p = new Vector();
            //Circle boundary
            n = vboundary.getCount();
            dist = 0;
            for (i = 1; i < n; i++)
            {
                dist += (vboundary[i].v - vboundary[i - 1].v).norm();
            }
            dist += (vboundary[n - 1].v - vboundary[0].v).norm();
            cd = 0;
            for (i = 0; i < n; i++)
            {
                p.z = 0;
                p.x = Math.Cos(cd / dist * 2.0 * Math.PI) * 0.5 + 0.5;
                p.y = Math.Sin(cd / dist * 2.0 * Math.PI) * 0.5 + 0.5;
                vboundary[i].tex = (Vector)p.Clone();
                if (i < n - 1)
                    cd += (vboundary[i + 1].v - vboundary[i].v).norm();
            }
        }

        public void RectangleBoundary(VertexSet vboundary)
        {
            int i, n, basei, corner;
            double min;
            Vector p = new Vector();
            //parameterise in a straight line
            n = vboundary.getCount();
            if (n == 0) return;
            p.x = p.y = p.z = 0.0;
            for (i = 0; i < n; i++)
            {
                p.x = i / ((double)n);
                vboundary[i].tex = (Vector)p.Clone();
            }
            //find point closest to 1/4 distance, make it the corner
            basei = 0;
            corner = 0;
            min = Math.Abs(vboundary[basei].tex.x - 1.0 / 4.0);
            for (i = basei; i < n; i++)
            {
                if (Math.Abs(vboundary[i].tex.x - 1.0 / 4.0) < min)
                {
                    min = Math.Abs(vboundary[i].tex.x - 1.0 / 4.0);
                    corner = i;
                }
            }
            for (i = basei; i < corner; i++)
            {
                p.x = p.z = 0.0;
                p.y = i / ((double)corner);
                vboundary[i].tex = (Vector)p.Clone();
            }

            //find point closest to 2/4 distance, make it the corner
            basei = corner;
            corner = basei;
            min = Math.Abs(vboundary[basei].tex.x - 2.0 / 4.0);
            for (i = basei; i < n; i++)
            {
                if (Math.Abs(vboundary[i].tex.x - 2.0 / 4.0) < min)
                {
                    min = Math.Abs(vboundary[i].tex.x - 2.0 / 4.0);
                    corner = i;
                }
            }
            for (i = basei; i < corner; i++)
            {
                p.z = 0.0;
                p.y = 1.0;
                p.x = (i - basei) / ((double)(corner - basei));
                vboundary[i].tex = (Vector)p.Clone();
            }

            //find point closest to 3/4 distance, make it the corner
            basei = corner;
            corner = basei;
            min = Math.Abs(vboundary[basei].tex.x - 3.0 / 4.0);
            for (i = basei; i < n; i++)
            {
                if (Math.Abs(vboundary[i].tex.x - 3.0 / 4.0) < min)
                {
                    min = Math.Abs(vboundary[i].tex.x - 3.0 / 4.0);
                    corner = i;
                }
            }
            for (i = basei; i < corner; i++)
            {
                p.z = 0.0;
                p.x = 1.0;
                p.y = 1.0 - (i - basei) / ((double)(corner - basei));
                vboundary[i].tex = (Vector)p.Clone();
            }

            //parameterise last stretch of boundary
            basei = corner;
            for (i = basei; i < n; i++)
            {
                p.z = 0.0;
                p.y = 0.0;
                p.x = 1.0 - (i - basei) / ((double)(n - basei));
                vboundary[i].tex = (Vector)p.Clone();
            }
        }

        public void RectangleChordBoundary(VertexSet vboundary)
        {
            int i, j, n, corner;
            double min;
            double dist, cd;
            Vector p = new Vector();
            //parameterise in a straight line
            n = vboundary.getCount();
            if (n == 0) return;
            p.x = p.y = p.z = 0.0;

            dist = 0;
            for (i = 1; i < n; i++)
            {
                dist += (vboundary[i].v - vboundary[i - 1].v).norm();
            }
            dist += (vboundary[n - 1].v - vboundary[0].v).norm();
            cd = 0;
            for (i = 0; i < n; i++)
            {
                p.x = cd / dist;
                vboundary[i].tex = (Vector)p.Clone();
                if (i < n - 1)
                    cd += (vboundary[i + 1].v - vboundary[i].v).norm();
            }

            //find corners
            //first corner is already at 0,0
            //last corner is the first corner
            for (j = 1; j < 4; j++)
            {
                //find closest parameter
                min = Math.Abs(vboundary[0].tex.x - j / 4.0);
                corner = 0;
                for (i = 1; i < n; i++)
                {
                    if (Math.Abs(vboundary[i].tex.x - j / 4.0) < min)
                    {
                        min = Math.Abs(vboundary[i].tex.x - j / 4.0);
                        corner = i;
                    }
                }
                //set the parameter value
                vboundary[corner].tex.x = j / 4.0;
            }
            //move each parameter to the correct segment
            for (i = 0; i < n; i++)
            {
                p = (Vector)vboundary[i].tex.Clone();
                //first segment
                if (p.x < 1.0 / 4.0)
                {
                    p.y = p.x / 0.25;
                    p.x = p.z = 0.0;
                }
                else if (p.x < 2.0 / 4.0)
                {
                    p.x = (p.x - 1.0 / 4.0) / 0.25;
                    p.y = 1.0;
                    p.z = 0.0;
                }
                else if (p.x < 3.0 / 4.0)
                {
                    p.y = 1.0 - (p.x - 2.0 / 4.0) / 0.25;
                    p.x = 1.0;
                    p.z = 0.0;
                }
                else
                {
                    p.x = 1.0 - (p.x - 3.0 / 4.0) / 0.25;
                    p.y = 0.0;
                    p.z = 0.0;
                }
                vboundary[i].tex = (Vector)p.Clone();
            }
        }

        public void ComputeBoundary(EdgeSet eboundary)
        {
            int i;

            //mark boundary edges
            for (i = 0; i < edges.getCount(); i++)
            {
                if (edges[i].faces.getCount() < 2)
                {
                    eboundary.add(edges[i]);
                }
            }
        }

        //Only works if there is precisely one boundary loop!!
        public void ComputeBoundary(EdgeSet eboundary, VertexSet vboundary)
        {
            int i;
            Edge e;
            Vertex v;

            //mark boundary edges
            Console.Write("Total edges: {0}\n", edges.getCount());
            for (i = 0; i < edges.getCount(); i++)
            {
                edges[i].label = 0;
                //mark boundaries
                if (edges[i].faces.getCount() < 2)
                {
                    edges[i].label |= BOUNDARY;
                    eboundary.add(edges[i]);
                }
            }
            Console.Write("Boundary edges: {0}\n", eboundary.getCount());

            //remove edges one for one to get boundary vertices
            //srand(time(null));
            //e=eboundary[rand()%eboundary.getCount()];
            if (eboundary.getCount() == 0) Console.Write("Fatal: bounray count=0\n");
            e = eboundary[0];
            v = e.ends[1];
            vboundary.add(e.ends[0]);
            vboundary.add(e.ends[1]);
            eboundary.remove(0);

            while (eboundary.getCount() > 0)
            {
                for (i = 0; i < eboundary.getCount(); i++)
                {
                    e = eboundary[i];
                    if (e.ends[0] == v)
                    {
                        vboundary.add(e.ends[1]);
                        eboundary.remove(i);
                        v = e.ends[1];
                        e = null;
                        break;
                    }
                    else if (e.ends[1] == v)
                    {
                        vboundary.add(e.ends[0]);
                        eboundary.remove(i);
                        v = e.ends[0];
                        e = null;
                        break;
                    }
                }
                if (e != null)
                    Console.Write("Broken loop! {0} edges left\n", eboundary.getCount());
                if (e != null)
                    break;
            }

        }

        //Tutte
        public void ParamBarycentric(int boundary_type)
        {
            int i, j, n, c;
            double lambda;
            EdgeSet eboundary = new EdgeSet();
            VertexSet vboundary = new VertexSet();
            Vertex v, v2;
            SparseMatrix M = new SparseMatrix();
            double[] bu;
            double[] bv;
            //for second chance to solve
            double[] scu;
            double[] scv;

            ComputeBoundary(eboundary, vboundary);
            if (vboundary.getCount() == 0)
            {
                Console.Write("No boundary, perhaps the mesh is too small?\n");
                return;
            }
            switch (boundary_type)
            {
                case BOUNDARY_CIRCLE:
                    CircleBoundary(vboundary);
                    break;
                case BOUNDARY_RECTANGLE:
                    RectangleBoundary(vboundary);
                    break;
                case BOUNDARY_CHORD_CIRCLE:
                    CircleChordBoundary(vboundary);
                    break;
                case BOUNDARY_CHORD_RECTANGLE:
                    RectangleChordBoundary(vboundary);
                    break;
                case BOUNDARY_PRESET:
                    break;
                default:
                    RectangleBoundary(vboundary);
                    break;
            }

            //label non-border vertices according to equation variable number
            n = 0;
            for (i = 0; i < vertices.getCount(); i++)
            {
                if (!vboundary.contains(vertices[i]))
                {
                    vertices[i].label = n;
                    n++;
                }
                else vertices[i].label = -1;
            }

            M.Resize(n, n);
            bu = new double[n];
            for (i = 0; i < n; i++) bu[i] = 0.0;
            bv = new double[n];
            for (i = 0; i < n; i++) bv[i] = 0.0;
            scu = new double[n];
            scv = new double[n];
            //border is paramerised. Set up system of equations setting each vertex to the
            //   barycentre of its neighbours (Tutte)
            for (i = 0; i < vertices.getCount(); i++)
            {
                //check for a boundary
                v = vertices[i];
                if (v.label != -1)
                {
                    c = v.label;
                    M[c, c] = 1.0;
                    //lambda=1.0/v.valence();
                    //Compute the 1 ring
                    VertexSet ring;
                    ring = ComputeExclusiveRing(v, 1);
                    //From Floater: lambda=w_xixj / sum(w_xixj)
                    // for barycentric w_xixj=1
                    lambda = 1.0 / (ring.getCount());

                    for (j = 0; j < ring.getCount(); j++)
                    {
                        v2 = ring[j];
                        if (v2 == v) continue;
                        //boundary goes into b
                        if (v2.label == -1)
                        {
                            bu[c] += lambda * v2.tex.x;
                            bv[c] += lambda * v2.tex.y;
                        }
                        else
                        {
                            M[c, v2.label] = -lambda;
                        }
                    }
                }
            }

            for (i = 0; i < n; i++) scu[i] = bu[i];
            for (i = 0; i < n; i++) scv[i] = bv[i];
            c = M.GaussSeidel(bu, 100000, 1e-7);
            if (c < 0)
            {
                Console.Write("Failed to solve\n");
                return;
            }
            if (c == 0)
            {
                //failed to obtain accuracy, try again
                c = M.GaussSeidel(scu, 100000, 1e-7, bu);
                if (c == 0)
                    Console.Write("Failed to solve to desired accuracy\n");
                for (i = 0; i < n; i++) bu[i] = scu[i];
            }
            c = M.GaussSeidel(bv, 100000, 1e-7);
            if (c < 0)
            {
                Console.Write("Failed to solve\n");
                return;
            }
            if (c == 0)
            {
                //failed to obtain accuracy, try again
                c = M.GaussSeidel(scv, 100000, 1e-7, bv);
                if (c == 0)
                    Console.Write("Failed to solve to desired accuracy\n");
                for (i = 0; i < n; i++) bv[i] = scv[i];
            }

            for (i = 0; i < vertices.getCount(); i++)
            {
                //check for a boundary
                v = vertices[i];
                if (v.label != -1)
                {
                    c = v.label;
                    v.tex.x = bu[c];
                    v.tex.y = bv[c];
                    v.tex.z = 0.0;
                }
            }
        }

        public double angle(Vector v1, Vector v2)
        {
            double proj, ang;
            v1 = (Vector)v1.Clone();
            v2 = (Vector)v2.Clone();
            v2.normalize();
            v1.normalize();
            proj = v1 * v2;
            //perp=(v1-proj*v2).norm();
            ang = Math.Acos(proj);
            return ang;
        }

        //Stretch for a vertex
        public double Sigma(Vertex v)
        {
            int i;
            double sumarea = 0.0, sumstretch = 0.0;
            double area, stretch;

            for (i = 0; i < v.faces.getCount(); i++)
            {
                area = v.faces[i].area();
                sumarea += area;
                stretch = v.faces[i].stretch();
                sumstretch += area * stretch * stretch;
            }

            return Math.Sqrt(sumstretch / sumarea);
        }

        //Stretch for entire mesh
        public double AllStretch()
        {
            int i;
            double sumarea = 0.0, sumstretch = 0.0;
            double area, stretch;

            for (i = 0; i < faces.getCount(); i++)
            {
                area = faces[i].area();
                sumarea += area;
                stretch = faces[i].stretch();
                sumstretch += area * stretch * stretch;
            }

            return Math.Sqrt(sumstretch / sumarea);
        }

        //determine if the vector from v1 to v2 and onwards crosses the vector from v3 to v4
        // 
        public bool cross(Vector v1, Vector v2, Vector v3, Vector v4)
        {
            Vector dir, e1, e2, p1, p2;
            dir = v2 - v1;
            dir.normalize();
            e1 = v3 - v2;
            e2 = v4 - v2;
            //check if either of the points are behind
            if (e1 * dir < 0) return false;
            if (e2 * dir < 0) return false;
            p1 = e1 - (e1 * dir) * dir;
            p2 = e2 - (e2 * dir) * dir;
            //don't care about positive crossings of v1--v2 and v3--v4, because we won't test these
            //check if the perpendicular vectors are on opposite sides of the direction
            if (p1 * p2 <= 0.0) return true;
            return false;
        }

        //area of triangle
        public double area(Vector v1, Vector v2, Vector v3)
        {
            return ((v2 - v1) ^ (v3 - v1)).norm() / 2.0;
        }

        public void ParamFloater(int boundary_type)
        {
            int i, j, k, n, c;
            double lambda, sumtheta, theta, totalw;
            EdgeSet eboundary = new EdgeSet();
            VertexSet vboundary = new VertexSet();
            Vertex v, v2 = null, v3 = null;
            SparseMatrix M = new SparseMatrix();
            double[] bu;
            double[] bv;
            //for second chance to solve
            double[] scu;
            double[] scv;
            bool trianglesleft = true;

            Console.Write("Boundary\n");
            ComputeBoundary(eboundary, vboundary);
            if (vboundary.getCount() == 0)
            {
                Console.Write("No boundary, perhaps the mesh is too small?\n");
                return;
            }
            Console.Write("Boundary param\n");
            switch (boundary_type)
            {
                case BOUNDARY_CIRCLE:
                    CircleBoundary(vboundary);
                    break;
                case BOUNDARY_RECTANGLE:
                    RectangleBoundary(vboundary);
                    break;
                case BOUNDARY_CHORD_CIRCLE:
                    CircleChordBoundary(vboundary);
                    break;
                case BOUNDARY_CHORD_RECTANGLE:
                    RectangleChordBoundary(vboundary);
                    break;
                case BOUNDARY_PRESET:
                    break;
                default:
                    RectangleBoundary(vboundary);
                    break;
            }

            //label non-border vertices according to equation variable number
            n = 0;
            for (i = 0; i < vertices.getCount(); i++)
            {
                if (!vboundary.contains(vertices[i]))
                {
                    vertices[i].label = n;
                    n++;
                }
                else vertices[i].label = -1;
            }

            Console.Write("Main param...\n");
            M.Resize(n, n);
            bu = new double[n];
            for (i = 0; i < n; i++) bu[i] = 0.0;
            bv = new double[n];
            for (i = 0; i < n; i++) bv[i] = 0.0;
            scu = new double[n];
            scv = new double[n];
            //border is paramerised. Set up system of equations setting each vertex to the
            //   barycentre of its neighbours (Tutte)
            for (i = 0; i < vertices.getCount(); i++)
            {
                //check for a boundary
                v = vertices[i];
                //corresponds to edges
                v.wij = new double[v.valence()];
                for (j = 0; j < v.valence(); j++)
                {
                    v.wij[j] = 0.0;
                }

                if (v.label != -1)
                {
                    VertexSet ring;
                    ring = ComputeExclusiveRing(v, 1);
                    //save texture coordinates of boundary in normals
                    for (j = 0; j < ring.getCount(); j++)
                    {
                        ring[j].n = (Vector)ring[j].tex.Clone();
                    }
                    c = v.label;
                    M[c, c] = 1.0;
                    //Now do local parameterisation according to Floater
                    //set p in the middle
                    v.tex = new Vector(0, 0, 0);
                    //determine sum of angles 
                    sumtheta = 0.0;
                    for (j = 0; j < v.faces.getCount(); j++)
                    {
                        if (v == v.faces[j].vertices[0])
                        {
                            v2 = v.faces[j].vertices[1];
                            v3 = v.faces[j].vertices[2];
                        }
                        else if (v == v.faces[j].vertices[1])
                        {
                            v2 = v.faces[j].vertices[0];
                            v3 = v.faces[j].vertices[2];
                        }
                        else if (v == v.faces[j].vertices[2])
                        {
                            v2 = v.faces[j].vertices[0];
                            v3 = v.faces[j].vertices[1];
                        }
                        sumtheta += Math.Abs(angle(v2.v - v.v, v3.v - v.v));
                        //this face has not been parameterised yet
                        v.faces[j].label = 0;
                    }

                    //place first adjacent vertex horizontally, preserving distance
                    j = 0;
                    if (v == v.faces[j].vertices[0])
                    {
                        v2 = v.faces[j].vertices[1];
                        v3 = v.faces[j].vertices[2];
                    }
                    else if (v == v.faces[j].vertices[1])
                    {
                        v2 = v.faces[j].vertices[0];
                        v3 = v.faces[j].vertices[2];
                    }
                    else if (v == v.faces[j].vertices[2])
                    {
                        v2 = v.faces[j].vertices[0];
                        v3 = v.faces[j].vertices[1];
                    }
                    theta = 0.0;
                    theta += 2 * Math.PI * Math.Abs(angle(v2.v - v.v, v3.v - v.v)) / sumtheta;
                    v2.tex.x = (v2.v - v.v).norm();
                    v2.tex.y = v2.tex.z = 0.0;
                    v3.tex.x = Math.Cos(theta);
                    v3.tex.y = Math.Sin(theta);
                    v3.tex.z = 0.0;
                    v3.tex *= (v3.v - v.v).norm();
                    v2 = v3;

                    //place vertices at shape preserving angles around the vertex
                    //this creates a local shape preserving parameterization
                    trianglesleft = true;
                    while (trianglesleft)
                    {
                        //find a triangle adjacent to the current one, that 
                        //has not been computed 
                        Edge e = new Edge();
                        e.ends[0] = v;
                        e.ends[1] = v2;

                        trianglesleft = false;
                        for (j = 0; j < v.faces.getCount(); j++)
                        {
                            if (v.faces[j].label == 0)
                            {
                                //check to see if common edge is shared
                                if (v.faces[j].edges[0].Equals(e)) break;
                                if (v.faces[j].edges[1].Equals(e)) break;
                                if (v.faces[j].edges[2].Equals(e)) break;
                            }
                        }
                        //found an adjacent triangle
                        if (j < v.faces.getCount())
                        {
                            int[] tag = new int[3] { 0, 0, 0 };
                            for (k = 0; k < 3; k++)
                            {
                                if (v == v.faces[j].vertices[k]) tag[k] = 1;
                                if (v2 == v.faces[j].vertices[k]) tag[k] = 1;
                            }
                            for (k = 0; k < 3; k++)
                            {
                                if (tag[k] == 0) v3 = v.faces[j].vertices[k];
                            }
                            //the undetermined vertex is placed according to
                            //its angle to the vertex, and the parametric distance
                            //   is set to be the same as the real distance
                            theta += 2 * Math.PI * Math.Abs(angle(v2.v - v.v, v3.v - v.v)) / sumtheta;
                            v3.tex.x = Math.Cos(theta);
                            v3.tex.y = Math.Sin(theta);
                            v3.tex.z = 0.0;
                            v3.tex *= (v3.v - v.v).norm();
                            v2 = v3;
                            trianglesleft = true;
                            v.faces[j].label = 1;
                        }
                    }

                    //for each vertex determine a triangle that  
                    //  v lies in, and then adjust the weights
                    //  according to the barycentric coordinates
                    //Note: ring.getCount()==v.valence()
                    for (k = 0; k < ring.getCount(); k++)
                    {
                        for (j = 0; j < v.faces.getCount(); j++)
                        {
                            if (v == v.faces[j].vertices[0])
                            {
                                v2 = v.faces[j].vertices[1];
                                v3 = v.faces[j].vertices[2];
                            }
                            else if (v == v.faces[j].vertices[1])
                            {
                                v2 = v.faces[j].vertices[0];
                                v3 = v.faces[j].vertices[2];
                            }
                            else if (v == v.faces[j].vertices[2])
                            {
                                v2 = v.faces[j].vertices[0];
                                v3 = v.faces[j].vertices[1];
                            }
                            //determine if this is the triangle
                            if (cross(ring[k].tex, v.tex, v2.tex, v3.tex))
                            {
                                double l1, l2, l3, a;
                                int i2, i3;
                                i2 = ring.index(v2);
                                i3 = ring.index(v3);
                                if (i2 < 0) Console.Write("Bad index: should not happen!\n");
                                if (i3 < 0) Console.Write("Bad index: should not happen!\n");
                                a = area(ring[k].tex, v2.tex, v3.tex);
                                l1 = area(v.tex, v2.tex, v3.tex) / a;
                                l2 = area(ring[k].tex, v.tex, v3.tex) / a;
                                l3 = area(ring[k].tex, v2.tex, v.tex) / a;
                                v.wij[k] += l1;
                                v.wij[i2] += l2;
                                v.wij[i3] += l3;
                            }
                        }
                    }

                    totalw = 0.0;
                    for (j = 0; j < ring.getCount(); j++)
                    {
                        totalw += v.wij[j];
                        //if it is a border vertex, put the normal
                        // back in the texture
                        if (ring[j].label == -1)
                        {
                            ring[j].tex = (Vector)ring[j].n.Clone();
                        }
                    }
                    //if (totalw==0.0) totalw=1.0;

                    for (j = 0; j < ring.getCount(); j++)
                    {
                        v2 = ring[j];
                        lambda = v.wij[j] / totalw;
                        if (v2 == v) continue;
                        //boundary goes into b
                        if (v2.label == -1)
                        {
                            bu[c] += lambda * v2.tex.x;
                            bv[c] += lambda * v2.tex.y;
                        }
                        else
                        {
                            M[c, v2.label] = -lambda;
                        }
                    }
                }
            }

            for (i = 0; i < n; i++) scu[i] = bu[i];
            for (i = 0; i < n; i++) scv[i] = bv[i];
            Console.Write("Solving\n");
            c = M.GaussSeidel(bu, 100000, 1e-7);
            //c=M.ConjugateGradient(bu, 100000, 1e-7);
            if (c < 0)
            {
                Console.Write("Failed to solve\n");
                return;
            }
            if (c == 0)
            {
                //failed to obtain accuracy, try again
                c = M.GaussSeidel(scu, 100000, 1e-7, bu);
                //c=M.ConjugateGradient(scu, 100000, 1e-7, bu);
                if (c == 0)
                    Console.Write("Failed to solve to desired accuracy\n");
                for (i = 0; i < n; i++) bu[i] = scu[i];
            }
            c = M.GaussSeidel(bv, 100000, 1e-7);
            //c=M.ConjugateGradient(bv, 100000, 1e-7);
            if (c < 0)
            {
                Console.Write("Failed to solve\n");
                return;
            }
            if (c == 0)
            {
                //failed to obtain accuracy, try again
                c = M.GaussSeidel(scv, 100000, 1e-7, bv);
                //c=M.ConjugateGradient(scv, 100000, 1e-7, bv);
                if (c == 0)
                    Console.Write("Failed to solve to desired accuracy\n");
                for (i = 0; i < n; i++) bv[i] = scv[i];
            }

            for (i = 0; i < vertices.getCount(); i++)
            {
                //check for a boundary
                v = vertices[i];
                if (v.label != -1)
                {
                    c = v.label;
                    v.tex.x = bu[c];
                    v.tex.y = bv[c];
                    v.tex.z = 0.0;
                }
            }
        }

        //Mean value parameterisation, also due to Floater
        public void ParamMeanValue(int boundary_type)
        {
            int i, j, k, n, c;
            double lambda, sumtheta, theta, totalw;
            EdgeSet eboundary = new EdgeSet();
            VertexSet vboundary = new VertexSet();
            Vertex v, v2 = null, v3 = null;
            SparseMatrix M = new SparseMatrix();
            double[] bu;
            double[] bv;
            //for second chance to solve
            double[] scu;
            double[] scv;
            bool trianglesleft = true;

            ComputeBoundary(eboundary, vboundary);
            if (vboundary.getCount() == 0)
            {
                Console.Write("No boundary, perhaps the mesh is too small?\n");
                return;
            }
            switch (boundary_type)
            {
                case BOUNDARY_CIRCLE:
                    CircleBoundary(vboundary);
                    break;
                case BOUNDARY_RECTANGLE:
                    RectangleBoundary(vboundary);
                    break;
                case BOUNDARY_CHORD_CIRCLE:
                    CircleChordBoundary(vboundary);
                    break;
                case BOUNDARY_CHORD_RECTANGLE:
                    RectangleChordBoundary(vboundary);
                    break;
                case BOUNDARY_PRESET:
                    break;
                default:
                    RectangleBoundary(vboundary);
                    break;
            }

            //label non-border vertices according to equation variable number
            n = 0;
            for (i = 0; i < vertices.getCount(); i++)
            {
                if (!vboundary.contains(vertices[i]))
                {
                    vertices[i].label = n;
                    n++;
                }
                else vertices[i].label = -1;
            }

            M.Resize(n, n);
            bu = new double[n];
            for (i = 0; i < n; i++) bu[i] = 0.0;
            bv = new double[n];
            for (i = 0; i < n; i++) bv[i] = 0.0;
            scu = new double[n];
            scv = new double[n];
            //border is paramerised. Set up system of equations setting each vertex to the
            //   barycentre of its neighbours (Tutte)
            for (i = 0; i < vertices.getCount(); i++)
            {
                //check for a boundary
                v = vertices[i];
                //corresponds to edges
                v.wij = new double[v.valence()];
                for (j = 0; j < v.valence(); j++)
                {
                    v.wij[j] = 0.0;
                }

                if (v.label != -1)
                {
                    VertexSet ring;
                    ring = ComputeExclusiveRing(v, 1);
                    //save texture coordinates of boundary in normals
                    for (j = 0; j < ring.getCount(); j++)
                    {
                        ring[j].normal = (Vector)ring[j].tex.Clone();
                    }
                    c = v.label;
                    M[c, c] = 1.0;
                    //Now do local parameterisation according to Floater
                    //set p in the middle
                    v.tex = new Vector(0, 0, 0);
                    //determine sum of angles 
                    sumtheta = 0.0;
                    for (j = 0; j < v.faces.getCount(); j++)
                    {
                        if (v == v.faces[j].vertices[0])
                        {
                            v2 = v.faces[j].vertices[1];
                            v3 = v.faces[j].vertices[2];
                        }
                        else if (v == v.faces[j].vertices[1])
                        {
                            v2 = v.faces[j].vertices[0];
                            v3 = v.faces[j].vertices[2];
                        }
                        else if (v == v.faces[j].vertices[2])
                        {
                            v2 = v.faces[j].vertices[0];
                            v3 = v.faces[j].vertices[1];
                        }
                        sumtheta += Math.Abs(angle(v2.v - v.v, v3.v - v.v));
                        //this face has not been parameterised yet
                        v.faces[j].label = 0;
                    }

                    //place first adjacent vertex horizontally, preserving distance
                    j = 0;
                    if (v == v.faces[j].vertices[0])
                    {
                        v2 = v.faces[j].vertices[1];
                        v3 = v.faces[j].vertices[2];
                    }
                    else if (v == v.faces[j].vertices[1])
                    {
                        v2 = v.faces[j].vertices[0];
                        v3 = v.faces[j].vertices[2];
                    }
                    else if (v == v.faces[j].vertices[2])
                    {
                        v2 = v.faces[j].vertices[0];
                        v3 = v.faces[j].vertices[1];
                    }
                    theta = 0.0;
                    theta += 2 * Math.PI * Math.Abs(angle(v2.v - v.v, v3.v - v.v)) / sumtheta;
                    v2.tex.x = (v2.v - v.v).norm();
                    v2.tex.y = v2.tex.z = 0.0;
                    v3.tex.x = Math.Cos(theta);
                    v3.tex.y = Math.Sin(theta);
                    v3.tex.z = 0.0;
                    v3.tex *= (v3.v - v.v).norm();
                    v2 = v3;

                    //place vertices at shape preserving angles around the vertex
                    //this creates a local shape preserving parameterization
                    trianglesleft = true;
                    while (trianglesleft)
                    {
                        //find a triangle adjacent to the current one, that 
                        //has not been computed 
                        Edge e = new Edge();
                        e.ends[0] = v;
                        e.ends[1] = v2;

                        trianglesleft = false;
                        for (j = 0; j < v.faces.getCount(); j++)
                        {
                            if (v.faces[j].label == 0)
                            {
                                //check to see if common edge is shared
                                if (v.faces[j].edges[0].Equals(e)) break;
                                if (v.faces[j].edges[1].Equals(e)) break;
                                if (v.faces[j].edges[2].Equals(e)) break;
                            }
                        }
                        //found an adjacent triangle
                        if (j < v.faces.getCount())
                        {
                            int[] tag = new int[3] { 0, 0, 0 };
                            for (k = 0; k < 3; k++)
                            {
                                if (v == v.faces[j].vertices[k]) tag[k] = 1;
                                if (v2 == v.faces[j].vertices[k]) tag[k] = 1;
                            }
                            for (k = 0; k < 3; k++)
                            {
                                if (tag[k] == 0) v3 = v.faces[j].vertices[k];
                            }
                            //the undetermined vertex is placed according to
                            //its angle to the vertex, and the parametric distance
                            //   is set to be the same as the real distance
                            theta += 2 * Math.PI * Math.Abs(angle(v2.v - v.v, v3.v - v.v)) / sumtheta;
                            v3.tex.x = Math.Cos(theta);
                            v3.tex.y = Math.Sin(theta);
                            v3.tex.z = 0.0;
                            v3.tex *= (v3.v - v.v).norm();
                            v2 = v3;
                            trianglesleft = true;
                            v.faces[j].label = 1;
                        }
                    }

                    for (k = 0; k < ring.getCount(); k++)
                    {
                        double alpha, beta;
                        Edge e;
                        Edge e2 = new Edge();

                        e2.ends[0] = v;
                        e2.ends[1] = v2;
                        v2 = ring[k];
                        e = v.edges.match(e2);
                        //find angles adjacent to edge
                        Face f = e.faces[0];
                        if (v == f.vertices[0])
                        {
                            v2 = f.vertices[1];
                            v3 = f.vertices[2];
                        }
                        else if (v == f.vertices[1])
                        {
                            v2 = f.vertices[0];
                            v3 = f.vertices[2];
                        }
                        else if (v == f.vertices[2])
                        {
                            v2 = f.vertices[0];
                            v3 = f.vertices[1];
                        }
                        alpha = angle(v3.v - v.v, v2.v - v.v);

                        f = e.faces[1];
                        if (v == f.vertices[0])
                        {
                            v2 = f.vertices[1];
                            v3 = f.vertices[2];
                        }
                        else if (v == f.vertices[1])
                        {
                            v2 = f.vertices[0];
                            v3 = f.vertices[2];
                        }
                        else if (v == f.vertices[2])
                        {
                            v2 = f.vertices[0];
                            v3 = f.vertices[1];
                        }
                        beta = angle(v3.v - v.v, v2.v - v.v);

                        //Formula from mean-value paper by Floater
                        v.wij[k] = (Math.Tan(alpha / 2.0) + Math.Tan(beta / 2.0)) / (v2.v - v.v).norm();
                    }

                    totalw = 0.0;
                    for (j = 0; j < ring.getCount(); j++)
                    {
                        totalw += v.wij[j];
                        //if it is a border vertex, put the normal
                        // back in the texture
                        if (ring[j].label == -1)
                        {
                            ring[j].tex = (Vector)ring[j].normal.Clone();
                        }
                    }
                    //if (totalw==0.0) totalw=1.0;

                    for (j = 0; j < ring.getCount(); j++)
                    {
                        v2 = ring[j];
                        lambda = v.wij[j] / totalw;
                        if (v2 == v) continue;
                        //boundary goes into b
                        if (v2.label == -1)
                        {
                            bu[c] += lambda * v2.tex.x;
                            bv[c] += lambda * v2.tex.y;
                        }
                        else
                        {
                            M[c, v2.label] = -lambda;
                        }
                    }
                }
            }

            for (i = 0; i < n; i++) scu[i] = bu[i];
            for (i = 0; i < n; i++) scv[i] = bv[i];
            c = M.GaussSeidel(bu, 100000, 1e-7);
            if (c < 0)
            {
                Console.Write("Failed to solve\n");
                return;
            }
            if (c == 0)
            {
                //failed to obtain accuracy, try again
                c = M.GaussSeidel(scu, 100000, 1e-7, bu);
                if (c == 0)
                    Console.Write("Failed to solve to desired accuracy\n");
                for (i = 0; i < n; i++) bu[i] = scu[i];
            }
            c = M.GaussSeidel(bv, 100000, 1e-7);
            if (c < 0)
            {
                Console.Write("Failed to solve\n");
                return;
            }
            if (c == 0)
            {
                //failed to obtain accuracy, try again
                c = M.GaussSeidel(scv, 100000, 1e-7, bv);
                if (c == 0)
                    Console.Write("Failed to solve to desired accuracy\n");
                for (i = 0; i < n; i++) bv[i] = scv[i];
            }

            for (i = 0; i < vertices.getCount(); i++)
            {
                //check for a boundary
                v = vertices[i];
                if (v.label != -1)
                {
                    c = v.label;
                    v.tex.x = bu[c];
                    v.tex.y = bv[c];
                    v.tex.z = 0.0;
                }
            }
        }
        //Use Floater parameterisation, followed by stretch minimisation of Yoshizawa, Belyaev and Seidel (2004)
        public void ParamStretchmin(int boundary_type, double eta)
        {
            int i, j, n, c;
            double lambda, sigma, totalw, Enew, Eold;
            EdgeSet eboundary = new EdgeSet();
            VertexSet vboundary = new VertexSet();
            Vertex v, v2;
            SparseMatrix M = new SparseMatrix();
            double[] bu;
            double[] bv;
            double[] su;
            double[] sv;
            //for second chance to solve
            double[] scu;
            double[] scv;
            int iter = 500;
            Mesh opt = null;

            ParamMeanValue(boundary_type);
            //ParamBarycentric(boundary_type);
            ComputeBoundary(eboundary, vboundary);
            if (vboundary.getCount() == 0)
            {
                Console.Write("No boundary, perhaps the mesh is too small?\n");
                return;
            }
            switch (boundary_type)
            {
                case BOUNDARY_CIRCLE:
                    CircleBoundary(vboundary);
                    break;
                case BOUNDARY_RECTANGLE:
                    RectangleBoundary(vboundary);
                    break;
                case BOUNDARY_CHORD_CIRCLE:
                    CircleChordBoundary(vboundary);
                    break;
                case BOUNDARY_CHORD_RECTANGLE:
                    RectangleChordBoundary(vboundary);
                    break;
                case BOUNDARY_PRESET:
                    break;
                default:
                    RectangleBoundary(vboundary);
                    break;
            }

            //label non-border vertices according to equation variable number
            n = 0;
            for (i = 0; i < vertices.getCount(); i++)
            {
                if (!vboundary.contains(vertices[i]))
                {
                    vertices[i].label = n;
                    n++;
                }
                else vertices[i].label = -1;
            }
            Console.Write("variables={0}\n", n);

            M.Resize(n, n);
            bu = new double[n];
            for (i = 0; i < n; i++) bu[i] = 0.0;
            bv = new double[n];
            for (i = 0; i < n; i++) bv[i] = 0.0;
            su = new double[n];
            sv = new double[n];
            scu = new double[n];
            scv = new double[n];

            Eold = AllStretch() + 1.0;
            while (--iter > 0)
            {
                Enew = AllStretch();
                Console.Write("Eold={0} Enew={1}\n", Eold, Enew);
                if ((Enew >= Eold) || (Double.IsNaN(Enew))) break;
                Eold = Enew;
                opt = this.Clone();
                Console.Write("Stretchmin iter {0}\n", iter);

                for (i = 0; i < n; i++) bu[i] = 0.0;
                for (i = 0; i < n; i++) bv[i] = 0.0;
                for (i = 0; i < vertices.getCount(); i++)
                {
                    //check for a boundary
                    v = vertices[i];

                    if (v.label != -1)
                    {
                        c = v.label;
                        M[c, c] = 1.0;
                        VertexSet ring;
                        ring = ComputeExclusiveRing(v, 1);
                        for (j = 0; j < ring.getCount(); j++)
                        {
                            sigma = Sigma(ring[j]);
                            sigma = Math.Pow(sigma, eta);
                            if (sigma != 0.0)
                                v.wij[j] /= sigma;
                            else
                                Console.Write("bad sigma\n");
                        }

                        totalw = 0.0;
                        for (j = 0; j < ring.getCount(); j++)
                        {
                            totalw += v.wij[j];
                        }
                        //if (totalw==0.0) {Console.Write("smallw\n"); totalw=1.0;}

                        for (j = 0; j < ring.getCount(); j++)
                        {
                            v2 = ring[j];
                            lambda = v.wij[j] / totalw;
                            if (lambda < 0) Console.Write("-lambda\n");
                            if (v2 == v) continue;
                            //boundary goes into b
                            if (v2.label == -1)
                            {
                                if (v2.tex.x < 0.0) Console.Write("smallx\n");
                                if (v2.tex.x > 1.0) Console.Write("largex\n");
                                if (v2.tex.y < 0.0) Console.Write("smally\n");
                                if (v2.tex.y > 1.0) Console.Write("largey\n");
                                bu[c] += lambda * v2.tex.x;
                                bv[c] += lambda * v2.tex.y;
                            }
                            else
                            {
                                M[c, v2.label] = -lambda;
                            }
                        }
                    }
                }

                for (i = 0; i < vertices.getCount(); i++)
                {
                    //check for a boundary
                    v = vertices[i];
                    if (v.label != -1)
                    {
                        c = v.label;
                        su[c] = v.tex.x;
                        sv[c] = v.tex.y;
                    }
                }

                for (i = 0; i < n; i++) scu[i] = bu[i];
                for (i = 0; i < n; i++) scv[i] = bv[i];
                //c=M.GaussSeidel(bu, 100000, 1e-7, su);
                c = M.GaussSeidel(bu, 100000, 1e-7);
                if (c < 0)
                {
                    Console.Write("Failed to solve\n");
                    copy(opt);
                    return;
                }
                if (c == 0)
                {
                    //failed to obtain accuracy, try again
                    c = M.GaussSeidel(scu, 100000, 1e-7, bu);
                    if (c == 0)
                    { Console.Write("Failed to solve to desired accuracy\n"); eta /= 2.0; continue; }
                    for (i = 0; i < n; i++) bu[i] = scu[i];
                }
                //c=M.GaussSeidel(bv, 100000, 1e-7, sv);
                c = M.GaussSeidel(bv, 100000, 1e-7);
                if (c < 0)
                {
                    Console.Write("Failed to solve\n");
                    copy(opt);
                    return;
                }
                if (c == 0)
                {
                    //failed to obtain accuracy, try again
                    c = M.GaussSeidel(scv, 100000, 1e-7, bv);
                    if (c == 0)
                    { Console.Write("Failed to solve to desired accuracy\n"); eta /= 2.0; continue; }
                    for (i = 0; i < n; i++) bv[i] = scv[i];
                }
                for (i = 0; i < vertices.getCount(); i++)
                {
                    //check for a boundary
                    v = vertices[i];
                    if (v.label != -1)
                    {
                        c = v.label;
                        v.tex.x = bu[c];
                        v.tex.y = bv[c];
                        v.tex.z = 0.0;
                    }
                }
            }
            Console.Write("Stretchmin done\n");
            copy(opt);
            Console.Write("opt assigned\n");
        }

        //Loop subdivision only works on triangle meshes
        //Loop subdivision only works on 2-manifold meshes
        public Mesh SubdivideLoop()
        {
            int i, j, k;
            Vector nv;
            Mesh nm = new Mesh();

            for (i = 0; i < vertices.getCount(); i++)
            {
                nv = new Vector(0.0, 0.0, 0.0);
                if (vertices[i].boundary())
                {
                    //find two boundary edges
                    // each opposite end contributes 1/8
                    // while this vertex contributes 6/8
                    for (j = 0; j < vertices[i].edges.getCount(); j++)
                    {
                        Edge e = vertices[i].edges[j];
                        if (e.boundary())
                        {
                            nv += 1.0 / 8.0 * e.ends[0].v;
                            nv += 1.0 / 8.0 * e.ends[1].v;
                        }
                    }
                    //this vertex has been added with weight 2*1.0/8.0
                    //still need 4.0/8.0
                    nv += 0.5 * vertices[i].v;
                }
                else
                {
                    VertexSet ring = ComputeExclusiveRing(vertices[i], 1);
                    //the valence
                    int n = vertices[i].valence(); //==ring.getCount()
                                                   //Loop's method
                                                   //double num=(3.0+2.0*Math.Cos(2.0*Math.PI/n));
                                                   //double beta=1.0/n*(5.0/8.0-num*num/64.0);
                                                   //Warren and Weimer's method			                   
                    double beta = 3.0 / (n * (n + 2.0));
                    nv = vertices[i].v * (1 - n * beta);
                    for (j = 0; j < n; j++)
                    {
                        nv += beta * ring[j].v;
                    }
                }
                vertices[i].label = nm.AddVertex(nv);
            }
            for (i = 0; i < edges.getCount(); i++)
            {
                nv = new Vector();
                Edge e = edges[i];
                if (e.boundary())
                {
                    nv = 0.5 * e.ends[0].v;
                    nv += 0.5 * e.ends[1].v;
                }
                else
                {
                    nv = 3.0 / 8.0 * e.ends[0].v;
                    nv += 3.0 / 8.0 * e.ends[1].v;
                    //for each face find the vertex
                    //not on the edge
                    for (k = 0; k < 2; k++)
                    {
                        nv += 1.0 / 8.0 * e.faces[k].OffEdge(e).v;
                    }
                }
                e.label = nm.AddVertex(nv);
            }

            //vertices are created, now add the triangles
            for (i = 0; i < faces.getCount(); i++)
            {
                Face f = faces[i];
                //center triangle
                nm.AddTriangle(f.edges[0].label, f.edges[1].label, f.edges[2].label);
                //triangles connected to existing vertices
                for (j = 0; j < 3; j++)
                {
                    for (k = j + 1; k < 3; k++)
                    {
                        Vertex v = null;
                        //find common vertex
                        if (f.edges[j].ends[0] == f.edges[k].ends[0])
                            v = f.edges[j].ends[0];
                        if (f.edges[j].ends[0] == f.edges[k].ends[1])
                            v = f.edges[j].ends[0];
                        if (f.edges[j].ends[1] == f.edges[k].ends[1])
                            v = f.edges[j].ends[1];
                        if (f.edges[j].ends[1] == f.edges[k].ends[0])
                            v = f.edges[j].ends[1];
                        nm.AddTriangle(f.edges[j].label, v.label, f.edges[k].label);
                    }
                }
            }

            return nm;
        }

        //this is an interpolating subdivision scheme
        // so existing vertices are left unaltered
        public Mesh SubdivideButterfly()
        {
            int i, j, k;
            Vector nv;
            Mesh nm = new Mesh();

            //add existing vertices with no changes
            for (i = 0; i < vertices.getCount(); i++)
            {
                vertices[i].label = nm.AddVertex(vertices[i].v);
            }
            //create new vertices for each edge
            for (i = 0; i < edges.getCount(); i++)
            {
                nv = new Vector();
                Edge e = edges[i];
                if (e.boundary())
                {
                    nv = 9.0 / 16.0 * e.ends[0].v;
                    nv += 9.0 / 16.0 * e.ends[1].v;
                    //find adjacent boundary edges
                    for (k = 0; k < 2; k++)
                    {
                        for (j = 0; j < e.ends[k].edges.getCount(); j++)
                        {
                            Edge eb = e.ends[k].edges[j];
                            if ((eb != e) && (eb.boundary()))
                            {
                                if (eb.ends[0] == e.ends[k])
                                    nv -= 1.0 / 16.0 * eb.ends[1].v;
                                else
                                    nv -= 1.0 / 16.0 * eb.ends[0].v;
                            }
                        }
                    }
                }
                else
                {
                    //number of irregular settings encountered
                    // to compute average
                    int irreg = 0;
                    //check for regular setting
                    if ((e.ends[0].valence() == 6)
                        && (e.ends[1].valence() == 6))
                    {
                        nv += 0.5 * e.ends[0].v;
                        nv += 0.5 * e.ends[1].v;
                        //for each face find the vertex
                        //not on the edge
                        for (k = 0; k < 2; k++)
                        {
                            nv += 1.0 / 8.0 * e.faces[k].OffEdge(e).v;
                        }
                        //Now find the other faces adjacent to the
                        // edge faces, and find the point not on the
                        // edge faces (4 of them).
                        for (k = 0; k < 2; k++)
                        {
                            for (j = 0; j < 3; j++)
                            {
                                if (e.faces[k].edges[j] != e)
                                {
                                    Face adj = null;
                                    if (e.faces[k] == e.faces[k].edges[j].faces[0])
                                        adj = e.faces[k].edges[j].faces[1];
                                    else
                                        adj = e.faces[k].edges[j].faces[0];
                                    if (adj != null)
                                        nv -= 1.0 / 16.0 * adj.OffEdge(e.faces[k].edges[j]).v;
                                    else
                                    {
                                        //use a reflection
                                    }
                                }
                            }
                        }

                        //since we will divide set irreg to 
                        // something not zero
                        irreg = 1;
                    }
                    //and now the irregular setting
                    for (k = 0; k < 2; k++)
                    {
                        if (e.ends[k].valence() != 6)
                        {
                            double sumwi = 0.0;
                            VertexSet ring = ComputeExclusiveRing(e.ends[k], 1);
                            ring = OrderRing(e.ends[k], e.ends[1 - k], ring);
                            if (e.ends[k].valence() == 3)
                            {
                                nv += 5.0 / 12.0 * ring[0].v; sumwi += 5.0 / 12.0;
                                nv -= 1.0 / 12.0 * ring[1].v; sumwi += -1.0 / 12.0;
                                nv -= 1.0 / 12.0 * ring[2].v; sumwi += -1.0 / 12.0;
                            }
                            else
                            if (e.ends[k].valence() == 4)
                            {
                                nv += 3.0 / 8.0 * ring[0].v; sumwi += 3.0 / 8.0;
                                //ring[1] weight is zero
                                nv -= 1.0 / 8.0 * ring[2].v; sumwi += -1.0 / 8.0;
                                //ring[3] weight is zero
                            }
                            else
                            {
                                int n = ring.getCount();
                                for (j = 0; j < ring.getCount(); j++)
                                {
                                    double wi = 0.25 + Math.Cos(2.0 * Math.PI * j / n) +
                                        0.5 * Math.Cos(4.0 * Math.PI * j / n);
                                    wi /= n;
                                    sumwi += wi;
                                    nv += wi * ring[j].v;
                                }
                            }
                            nv += (1.0 - sumwi) * e.ends[k].v;
                            irreg++;
                        }
                    }
                    nv /= irreg;
                }
                e.label = nm.AddVertex(nv);
            }

            //vertices are created, now add the triangles
            for (i = 0; i < faces.getCount(); i++)
            {
                Face f = faces[i];
                //center triangle
                nm.AddTriangle(f.edges[0].label, f.edges[1].label, f.edges[2].label);
                //triangles connected to existing vertices
                for (j = 0; j < 3; j++)
                {
                    for (k = j + 1; k < 3; k++)
                    {
                        Vertex v = null;
                        //find common vertex
                        if (f.edges[j].ends[0] == f.edges[k].ends[0])
                            v = f.edges[j].ends[0];
                        if (f.edges[j].ends[0] == f.edges[k].ends[1])
                            v = f.edges[j].ends[0];
                        if (f.edges[j].ends[1] == f.edges[k].ends[1])
                            v = f.edges[j].ends[1];
                        if (f.edges[j].ends[1] == f.edges[k].ends[0])
                            v = f.edges[j].ends[1];
                        nm.AddTriangle(f.edges[j].label, v.label, f.edges[k].label);
                    }
                }
            }

            return nm;
        }

        //m is the vertex for which ring is the 1-ring
        //start is the first vertex in the ordered ring
        public VertexSet OrderRing(Vertex m, Vertex start, VertexSet ring)
        {
            int i, j;
            bool found;
            int miss = 0;
            Vertex v = start;
            VertexSet r = new VertexSet();
            ring.remove(start);
            r.add(start);

            while (ring.getCount() > 0)
            {
                found = false;
                for (i = 0; i < m.faces.getCount(); i++)
                {
                    for (j = 0; j < m.faces[i].edges.getCount(); j++)
                    {
                        Edge e = m.faces[i].edges[j];
                        if ((e.ends[0] == v) && (ring.contains(e.ends[1])))
                        {
                            v = e.ends[1];
                            r.add(v);
                            ring.remove(v);
                            found = true;
                            break;
                        }
                        else if ((e.ends[1] == v) && (ring.contains(e.ends[0])))
                        {
                            v = e.ends[0];
                            r.add(v);
                            ring.remove(v);
                            found = true;
                            break;
                        }
                    }
                    if (found) break;
                }
                //no edge found, this could be a problem
                if (!found)
                    miss++;
                //the ring of vertices is probably broken
                // by a boundary, we should use reflection
                // to create virtual interpolation points.
                // instead we just look for the other edge of the boundary
                //it is not sufficient that this is a boundary point,
                // it must also only have one face in common with the
                // vertex for which the ring has been computed
                if (miss > 5)
                {
                    miss = 0;
                    for (i = 0; i < m.edges.getCount(); i++)
                    {
                        if (m.edges[i].boundary())
                        {
                            if (m.edges[i].ends[0] == m)
                                v = m.edges[i].ends[1];
                            else
                                v = m.edges[i].ends[0];
                            if (ring.contains(v))
                            {
                                r.add(v);
                                ring.remove(v);
                                found = true;
                                break;
                            }
                        }
                    }
                }
            }
            return r;
        }

        //boundary specifies whether the boundary edge must be subdivided
        //this must only happen every second iteration
        public Mesh SubdivideSqrt3(bool boundary)
        {
            int i, j, k;
            int v1, v2, v3, v4, v5;
            Mesh nm = new Mesh();
            //edge flipping happens as a second pass
            //so we have two meshes
            Mesh nm2 = new Mesh();
            Vector nv;
            //smooth exisitng vertices
            for (i = 0; i < vertices.getCount(); i++)
            {
                VertexSet ring;
                nv = new Vector();
                if (vertices[i].boundary())
                {
                    if (boundary)
                    {
                        //find two boundary edges
                        // each opposite end contributes 4/27
                        // while this vertex contributes 19/27
                        for (j = 0; j < vertices[i].edges.getCount(); j++)
                        {
                            Edge e = vertices[i].edges[j];
                            if (e.boundary())
                            {
                                nv += 4.0 / 27.0 * e.ends[0].v;
                                nv += 4.0 / 27.0 * e.ends[1].v;
                            }
                        }
                        //this vertex has been added with weight 2*4.0/27.0
                        //still need 11.0/27.0
                        nv += 11.0 / 27.0 * vertices[i].v;
                    }
                    else
                    {
                        nv = (Vector)vertices[i].v.Clone();
                    }
                }
                else
                {
                    int n = vertices[i].valence();
                    double alpha = 4.0 - 2.0 * Math.Cos(2.0 * Math.PI / n);
                    alpha /= 9.0;
                    nv = (1.0 - alpha) * vertices[i].v;
                    ring = ComputeExclusiveRing(vertices[i], 1);
                    for (j = 0; j < ring.getCount(); j++)
                    {
                        nv += alpha / n * ring[j].v;
                    }
                }
                vertices[i].label = nm.AddVertex(nv);
                //label this as : not a new vertex
                nm.vertices[vertices[i].label].label = 0;
            }
            //add a new vertex for each face
            for (i = 0; i < faces.getCount(); i++)
            {
                nv = new Vector();
                if ((faces[i].boundary()) && (boundary))
                {
                    Edge e = null;
                    //there should be precisely one boundary edge
                    //in the first pass any triangle with two boundary
                    // edges is subdivided so that each new triangle
                    // has only one boundary edge. Since this code
                    // should only be invoked every second iteration
                    // only one bounadry edge should be present.
                    for (j = 0; j < faces[i].edges.getCount(); j++)
                    {
                        if (faces[i].edges[j].boundary())
                            e = faces[i].edges[j];
                    }
                    //split the boundary edge
                    //find two boundary edges
                    nv = new Vector();
                    for (k = 0; k < e.ends[0].edges.getCount(); k++)
                    {
                        Edge e2 = e.ends[0].edges[k];
                        if (e2.boundary() && (e2 != e))
                        {
                            if (e2.ends[0] != e.ends[0])
                            {
                                nv = 1.0 / 27.0 * (e2.ends[0].v +
                                             16.0 * e.ends[0].v +
                                             10.0 * e.ends[1].v);
                            }
                            else
                            {
                                nv = 1.0 / 27.0 * (e2.ends[1].v +
                                             16.0 * e.ends[0].v +
                                             10.0 * e.ends[1].v);
                            }
                        }
                    }
                    //nv=e.ends[0].v*2.0/3.0+e.ends[1].v*1.0/3.0;
                    v4 = nm.AddVertex(nv);

                    nv = new Vector();
                    for (k = 0; k < e.ends[1].edges.getCount(); k++)
                    {
                        Edge e2 = e.ends[1].edges[k];
                        if (e2.boundary() && (e2 != e))
                        {
                            if (e2.ends[0] != e.ends[1])
                            {
                                nv = 1.0 / 27.0 * (e2.ends[0].v +
                                             16.0 * e.ends[1].v +
                                             10.0 * e.ends[0].v);
                            }
                            else
                            {
                                nv = 1.0 / 27.0 * (e2.ends[1].v +
                                             16.0 * e.ends[1].v +
                                             10.0 * e.ends[0].v);
                            }
                        }
                    }
                    //nv=e.ends[0].v*1.0/3.0+e.ends[1].v*2.0/3.0;
                    v5 = nm.AddVertex(nv);
                    //since vertices are added sequntially we have v5=v4+1
                    faces[i].label = v4;
                    nm.vertices[v4].label = 1;
                    nm.vertices[v5].label = 1;
                }
                else
                {
                    nv = faces[i].vertices[0].v / 3.0;
                    nv += faces[i].vertices[1].v / 3.0;
                    nv += faces[i].vertices[2].v / 3.0;
                    faces[i].label = nm.AddVertex(nv);
                    //label this as : a new vertex
                    nm.vertices[faces[i].label].label = 1;
                }
            }
            //now create the new triangles
            //edge flips are done later
            for (i = 0; i < faces.getCount(); i++)
            {
                if ((faces[i].boundary()) && (boundary))
                {
                    //boundary edges are only subdivided
                    // every second iteration and are not subject
                    // to edge flips.
                    Edge e = null;
                    //there should be precisely one boundary edge
                    //in the first pass any triangle with two boundary
                    // edges is subdivided so that each new triangle
                    // has only one boundary edge. Since this code
                    // should only be invoked every second iteration
                    // only one bounadry edge should be present.
                    for (j = 0; j < faces[i].edges.getCount(); j++)
                    {
                        if (faces[i].edges[j].boundary())
                            e = faces[i].edges[j];
                    }
                    v1 = faces[i].OffEdge(e).label;
                    v2 = e.ends[0].label;
                    v3 = e.ends[1].label;
                    v4 = faces[i].label;
                    v5 = v4 + 1;
                    nm.AddTriangle(v1, v2, v4);
                    nm.AddTriangle(v1, v4, v5);
                    nm.AddTriangle(v1, v5, v3);
                }
                else
                {
                    v4 = faces[i].label;
                    v1 = faces[i].vertices[0].label;
                    v2 = faces[i].vertices[1].label;
                    v3 = faces[i].vertices[2].label;
                    nm.AddTriangle(v1, v2, v4);
                    nm.AddTriangle(v2, v3, v4);
                    nm.AddTriangle(v3, v1, v4);
                }
            }

            //now do edge flipping, create a new mesh with the
            // same vertices
            for (i = 0; i < nm.vertices.getCount(); i++)
            {
                //add the vertex and store the index
                nm.vertices[i].pos = nm2.AddVertex(nm.vertices[i].v);
            }
            //now add the triangles in such a way that the edges
            // are flipped
            for (i = 0; i < nm.faces.getCount(); i++)
            {
                //not added yet
                nm.faces[i].label = 0;
            }
            for (i = 0; i < nm.edges.getCount(); i++)
            {
                //only non boundary edges need to be flipped
                Edge e = nm.edges[i];
                if (!e.boundary())
                {
                    //only flip if both ends of the edge are "old vertices"
                    if ((e.ends[0].label == 0) &&
                        (e.ends[1].label == 0))
                    {
                        //and if oppposite ends are "new vertices"
                        if ((e.faces[0].OffEdge(e).label == 1)
                            && (e.faces[1].OffEdge(e).label == 1))
                        {
                            v1 = e.ends[0].pos;
                            v2 = e.ends[1].pos;
                            v3 = e.faces[0].OffEdge(e).pos;
                            v4 = e.faces[1].OffEdge(e).pos;
                            nm2.AddTriangle(v1, v3, v4);
                            nm2.AddTriangle(v3, v2, v4);
                            //these triangles are added
                            e.faces[0].label = 1;
                            e.faces[1].label = 1;
                        }
                    }
                }
            }
            for (i = 0; i < nm.faces.getCount(); i++)
            {
                if (nm.faces[i].label == 0)
                {
                    v1 = nm.faces[i].vertices[0].pos;
                    v2 = nm.faces[i].vertices[1].pos;
                    v3 = nm.faces[i].vertices[2].pos;
                    nm2.AddTriangle(v1, v2, v3);
                }
            }
            return nm2;
        }


        public void AlignAdjacentFaces(FaceSet queue)
        {
            int i;
            Face f;

            while (queue.getCount() > 0)
            {
                f = queue[0];
                queue.remove(0);
                if (f.label == 0) return;
                for (i = 0; i < f.edges.getCount(); i++)
                {
                    Edge e = f.edges[i];
                    Face adj;
                    if (!e.boundary())
                    {
                        if (e.faces[0] == f) adj = e.faces[1];
                        else adj = e.faces[0];
                        //might cause acute angles to face inward!
                        if ((adj.normal * f.normal < 0.0) && (adj.label == 0))
                        {
                            //flip adj to have same orientation
                            // as this face;
                            Vertex tmp;
                            tmp = adj.vertices[0];
                            adj.vertices[0] = adj.vertices[2];
                            adj.vertices[2] = tmp;
                            adj.normal = -adj.normal;
                        }
                        if (adj.label == 0)
                        {
                            adj.label = 1;
                            queue.add(adj);
                        }
                    }
                }
            }
        }
        //try to fix problems in alignment
        // of normals caused by incorrect
        // triangle vertex orders
        //Only works on triangles
        public void AlignNormals()
        {
            int i;
            ComputeFaceNormals();
            for (i = 0; i < faces.getCount(); i++)
            {
                faces[i].label = 0;
            }
            faces[0].label = 1;
            FaceSet queue = new FaceSet();
            queue.add(faces[0]);
            AlignAdjacentFaces(queue);
            ComputeVertexNormals();
        }

        public void SaveTexInNormal()
        {
            int i;
            for (i = 0; i < vertices.getCount(); i++)
            {
                vertices[i].normal = vertices[i].tex;
            }
        }

        public void WriteOFF(string filename)
        {
            int nv, nf;
            int i, v1, v2, v3;
            StreamWriter f = new StreamWriter(filename);
            if (f == null) return;

            f.Write("OFF\n");
            nv = vertices.getCount();
            nf = faces.getCount();
            f.Write("{0} {1} 0\n", nv, nf);
            for (i = 0; i < nv; i++)
            {
                f.Write("{0} {1} {2}\n", vertices[i].v.x, vertices[i].v.y, vertices[i].v.z);
                vertices[i].label = i;
            }
            for (i = 0; i < nf; i++)
            {
                v1 = faces[i].vertices[0].label;
                v2 = faces[i].vertices[1].label;
                v3 = faces[i].vertices[2].label;
                f.Write("3 {0} {1} {2}\n", v1, v2, v3);
            }
            f.Close();
        }

        //wavefront obj
        public void WriteOBJ(string filename)
        {
            int nv, nf;
            int i, v1, v2, v3;
            StreamWriter f = new StreamWriter(filename);
            if (f == null) return;

            nv = vertices.getCount();
            nf = faces.getCount();
            for (i = 0; i < nv; i++)
            {
                f.Write("v {0} {1} {2}\n", vertices[i].v.x, vertices[i].v.y, vertices[i].v.z);
                vertices[i].label = i + 1;
            }
            for (i = 0; i < nf; i++)
            {
                v1 = faces[i].vertices[0].label;
                v2 = faces[i].vertices[1].label;
                v3 = faces[i].vertices[2].label;
                f.Write("f {0} {1} {2}\n", v1, v2, v3);
            }
            f.Close();
        }

        public void WriteParamOBJ(string filename)
        {
            int nv, nf;
            int i, v1, v2, v3;
            StreamWriter f = new StreamWriter(filename);
            if (f == null) return;

            nv = vertices.getCount();
            nf = faces.getCount();
            for (i = 0; i < nv; i++)
            {
                f.Write("v {0} {1} {2}\n", vertices[i].tex.x, vertices[i].tex.y, vertices[i].tex.z);
                vertices[i].label = i + 1;
            }
            for (i = 0; i < nf; i++)
            {
                v1 = faces[i].vertices[0].label;
                v2 = faces[i].vertices[1].label;
                v3 = faces[i].vertices[2].label;
                f.Write("f {0} {1} {2}\n", v1, v2, v3);
            }
            f.Close();
        }

        public void WriteParamOFF(string filename)
        {
            int nv, nf;
            int i, v1, v2, v3;
            StreamWriter f = new StreamWriter(filename);
            if (f == null) return;

            f.Write("OFF\n");
            nv = vertices.getCount();
            nf = faces.getCount();
            f.Write("{0} {1} 0\n", nv, nf);
            for (i = 0; i < nv; i++)
            {
                f.Write("{0} {1} {2}\n", vertices[i].tex.x, vertices[i].tex.y, vertices[i].tex.z);
                vertices[i].label = i;
            }
            for (i = 0; i < nf; i++)
            {
                v1 = faces[i].vertices[0].label;
                v2 = faces[i].vertices[1].label;
                v3 = faces[i].vertices[2].label;
                f.Write("3 {0} {1} {2}\n", v1, v2, v3);
            }
            f.Close();
        }

        //size is the recommended size, but it may have to be larger to avoid
        //	vertices being placed in the same cell
        public Vector[,] CreateGeometryImage(int size, bool resize)
        {
            return CreateGeometryImage(size, resize, GEOMETRY_IMAGE);
        }
        public Vector[,] CreateGeometryImage(int size, bool resize, int type)
        {
            int[,] occupied;
            Vector[,] gim = null;
            int i, j, x, y, newsize, closest;
            bool restart = true;
            double minx, maxx, miny, maxy, tx, ty;
            double alpha, beta, gamma;
            double basearea, min;
            Face f;

            newsize = size;
            while (restart)
            {
                restart = false;
                size = newsize;
                //no vertices have been placed
                occupied = new int[size, size];
                for (y = 0; y < size; y++)
                {
                    for (x = 0; x < size; x++)
                    {
                        occupied[y, x] = -1;
                    }
                }

                //place the vertices in cells				
                for (i = 0; i < vertices.getCount(); i++)
                {
                    x = (int)(vertices[i].tex.x * (size - 1) + 0.4);
                    y = (int)(vertices[i].tex.y * (size - 1) + 0.4);
                    //if not occupied, then use the vertex
                    if (occupied[y, x] == -1)
                    {
                        occupied[y, x] = i;
                    }
                    else
                    {
                        //occupied cell, either don't use the
                        //   current vertex, or resize and start again
                        if ((resize) && (size < 1024))
                        {
                            newsize = size * 2;
                            restart = true;
                            break;
                        }
                    }
                }

                //if there is no need to restart, then search through
                // cells for unoccupied cells
                if (!restart)
                {
                    gim = new Vector[size, size];
                    for (y = 0; y < size; y++)
                    {
                        FaceSet linecandidates = new FaceSet();
                        ty = ((y) / (double)(size - 1));
                        if (y == 0) ty = 0.0;
                        if (y == size - 1) ty = 1.0;

                        for (i = 0; i < faces.getCount(); i++)
                        {
                            //compute the bounding box of the triangle
                            //minx=maxx=faces[i].vertices[0].tex.x;
                            miny = maxy = faces[i].vertices[0].tex.y;
                            for (j = 1; j < 3; j++)
                            {
                                //if (faces[i].vertices[j].tex.x>maxx)
                                //	maxx=faces[i].vertices[j].tex.x;
                                //if (faces[i].vertices[j].tex.x<minx)
                                //	minx=faces[i].vertices[j].tex.x;
                                if (faces[i].vertices[j].tex.y > maxy)
                                    maxy = faces[i].vertices[j].tex.y;
                                if (faces[i].vertices[j].tex.y < miny)
                                    miny = faces[i].vertices[j].tex.y;
                            }
                            //check to see if it is in the bounding box
                            if (/*(minx<=tx)&&(tx<=maxx)&&*/(miny <= ty) && (ty <= maxy))
                            {
                                linecandidates.add(faces[i]);
                            }
                        }

                        for (x = 0; x < size; x++)
                        {
                            if (occupied[y, x] != -1)
                            {
                                if (type == NORMAL_MAP)
                                {
                                    gim[y, x] = vertices[occupied[y, x]].normal;
                                }
                                else
                                {
                                    gim[y, x] = vertices[occupied[y, x]].v;
                                }
                            }
                            else
                            if (occupied[y, x] == -1)
                            {
                                //take the unnocupied cell	
                                // and compute in which triangle the
                                // centre of the cell lies
                                FaceSet candidates = new FaceSet();
                                //tx=((x+0.5)/(double)size);
                                //ty=((y+0.5)/(double)size);
                                tx = ((x) / (double)(size - 1));
                                //ty=((y)/(double)(size-1));
                                // handle borders specially
                                if (x == 0) tx = 0.0;
                                if (x == size - 1) tx = 1.0;
                                //if (y==0) ty=0.0;
                                //if (y==size-1) ty=1.0;
                                Vector v = new Vector();
                                v.x = tx; v.y = ty; v.z = 0.0;
                                for (i = 0; i < linecandidates.getCount(); i++)
                                {
                                    //compute the bounding box of the triangle
                                    minx = maxx = linecandidates[i].vertices[0].tex.x;
                                    //miny=maxy=faces[i].vertices[0].tex.y;
                                    for (j = 1; j < 3; j++)
                                    {
                                        if (linecandidates[i].vertices[j].tex.x > maxx)
                                            maxx = linecandidates[i].vertices[j].tex.x;
                                        if (linecandidates[i].vertices[j].tex.x < minx)
                                            minx = linecandidates[i].vertices[j].tex.x;
                                        //if (faces[i].vertices[j].tex.y>maxy)
                                        //	maxy=faces[i].vertices[j].tex.y;
                                        //if (faces[i].vertices[j].tex.y<miny)
                                        //	miny=faces[i].vertices[j].tex.y;
                                    }
                                    //check to see if it is in the bounding box
                                    if ((minx <= tx) && (tx <= maxx)/*&&(miny<=ty)&&(ty<=maxy)*/)
                                    {
                                        //candidates.add(faces[i]);
                                        candidates.add(linecandidates[i]);
                                    }
                                }
                                closest = -1;
                                min = 1e10;
                                Vector v1, v2, v3;
                                for (i = 0; i < candidates.getCount(); i++)
                                {
                                    //compute barycentric coordinates
                                    // and use them to determine the 
                                    // euclidean 3d-space position
                                    //the point is inside the triangle if
                                    // the sum of the barycentric coordinates is 1
                                    // and all coordinates are greater than 0.
                                    f = candidates[i];
                                    v1 = f.vertices[0].tex;
                                    v2 = f.vertices[1].tex;
                                    v3 = f.vertices[2].tex;
                                    basearea = area(v1, v2, v3);
                                    alpha = area(v, v2, v3) / basearea;
                                    beta = area(v1, v, v3) / basearea;
                                    gamma = area(v1, v2, v) / basearea;
                                    //alpha, beta and gamma are all positive by definition	
                                    //test if they add to 1.0
                                    //if (Math.Abs(alpha+beta+gamma-1.0)<1e-7) {
                                    //	v=alpha*f.vertices[0].v+
                                    //		beta*f.vertices[1].v+	
                                    //		gamma*f.vertices[2].v;	
                                    //} else {
                                    if (Math.Abs(alpha + beta + gamma - 1.0) < min)
                                    {
                                        closest = i;
                                        min = Math.Abs(alpha + beta + gamma - 1.0);
                                    }
                                    //}
                                }
                                f = candidates[closest];
                                v1 = f.vertices[0].tex;
                                v2 = f.vertices[1].tex;
                                v3 = f.vertices[2].tex;
                                basearea = area(v1, v2, v3);
                                alpha = area(v, v2, v3) / basearea;
                                beta = area(v1, v, v3) / basearea;
                                gamma = area(v1, v2, v) / basearea;
                                if (type == NORMAL_MAP)
                                {
                                    v = alpha * f.vertices[0].normal +
                                        beta * f.vertices[1].normal +
                                        gamma * f.vertices[2].normal;
                                }
                                else
                                {
                                    v = alpha * f.vertices[0].v +
                                        beta * f.vertices[1].v +
                                        gamma * f.vertices[2].v;
                                }
                                gim[y, x] = v;
                            }
                        }
                    }
                }
            }

            return gim;
        }

        public void WriteGeometryImageRIB(string filename, string tiffname, int size)
        {
            int i;
            StreamWriter f = new StreamWriter(filename);
            Vector v1, v2, v3, v4;

            if (f == null) return;
            f.Write("Format {0} {1} 1\n", size, size);
            f.Write("Clipping 0.1 100.0\n");
            f.Write("PixelSamples 2 2\n");
            f.Write("Sides 2\n");
            f.Write("ShadingRate 0.5\n");
            f.Write("Quantize \"rgba\" 0 0 0 0\n");
            f.Write("Display \"{0}\" \"file\" \"rgba\"\n", tiffname);
            f.Write("Scale {0} {1} 1\n", 2.0 - 1.0 / size, 2.0 - 1.0 / size);
            f.Write("Translate 0 0 1\n");
            f.Write("WorldBegin\n");
            f.Write("\tSurface \"constant\"\n");

            for (i = 0; i < vertices.getCount(); i++)
            {
                vertices[i].tex.x -= 0.5;
                vertices[i].tex.y -= 0.5;
                vertices[i].label = i;
            }
            f.Write("\t PointsPolygons [\n");
            for (i = 0; i < faces.getCount(); i++)
            {
                f.Write("\t\t 3\n");
            }
            f.Write("\t ] [\n");
            for (i = 0; i < faces.getCount(); i++)
            {
                f.Write("\t\t {0} {1} {2}\n", faces[i].vertices[0].label,
                        faces[i].vertices[1].label,
                        faces[i].vertices[2].label);
            }
            f.Write("\t ]\n");
            f.Write("\t \"P\" [\n");
            for (i = 0; i < vertices.getCount(); i++)
            {
                f.Write("\t\t {0} {1} {2} \n",
                    vertices[i].tex.x,
                    vertices[i].tex.y,
                    vertices[i].tex.z);
            }
            f.Write("\t ]\n");
            f.Write("\t \"Cs\" [\n");
            for (i = 0; i < vertices.getCount(); i++)
            {
                f.Write("\t\t {0} {1} {2} \n",
                    vertices[i].v.x,
                    vertices[i].v.y,
                    vertices[i].v.z);
            }
            f.Write("\t ]\n");
            f.Write("\tTranslate 0.0 0.0 0.5\n");
            for (i = 0; i < edges.getCount(); i++)
            {
                if (edges[i].boundary())
                {
                    v1 = edges[i].ends[0].tex * 0.5;
                    v2 = edges[i].ends[0].tex * 1.5;
                    v3 = edges[i].ends[1].tex * 1.5;
                    v4 = edges[i].ends[1].tex * 0.5;
                    f.Write("\t Polygon \"Cs\" [ ");
                    f.Write("{0} {1} {2} ",
                        edges[i].ends[0].v.x,
                        edges[i].ends[0].v.y,
                        edges[i].ends[0].v.z);
                    f.Write("{0} {1} {2} ",
                        edges[i].ends[0].v.x,
                        edges[i].ends[0].v.y,
                        edges[i].ends[0].v.z);
                    f.Write("{0} {1} {2} ",
                        edges[i].ends[1].v.x,
                        edges[i].ends[1].v.y,
                        edges[i].ends[1].v.z);
                    f.Write("{0} {1} {2} ",
                        edges[i].ends[1].v.x,
                        edges[i].ends[1].v.y,
                        edges[i].ends[1].v.z);
                    f.Write("]\n\t\t \"P\" [ ");
                    f.Write("{0} {1} {2} ", v1.x, v1.y, v1.z);
                    f.Write("{0} {1} {2} ", v2.x, v2.y, v2.z);
                    f.Write("{0} {1} {2} ", v3.x, v3.y, v3.z);
                    f.Write("{0} {1} {2} ", v4.x, v4.y, v4.z);
                    f.Write("]\n");
                }
            }
            f.Write("WorldEnd\n");
            f.Close();

            for (i = 0; i < vertices.getCount(); i++)
            {
                vertices[i].tex.x += 0.5;
                vertices[i].tex.y += 0.5;
            }
        }

        public void WriteParamMetaPost(string filename)
        {
            int i;
            StreamWriter f = new StreamWriter(filename);
            if (f == null) return;
            Vector v1, v2;

            f.Write("u:=5cm;\n");
            f.Write("beginfig(1);\n");
            for (i = 0; i < edges.getCount(); i++)
            {
                v1 = edges[i].ends[0].tex;
                v2 = edges[i].ends[1].tex;
                f.Write("\tdraw ({0},{1})--({0},{1});\n", v1.x, v1.y, v2.x, v2.y);
            }
            f.Write("endfig;\n");
            f.Write("end;\n");
            f.Close();
        }

        public void WriteNormalMapRIB(string filename, string tiffname, int size)
        {
            int i;

            StreamWriter f = new StreamWriter(filename);
            Vector v1, v2, v3, v4;
            if (f == null) return;
            f.Write("Format {0} {1} 1\n", size, size);
            f.Write("Clipping 0.1 100.0\n");
            f.Write("PixelSamples 2 2\n");
            f.Write("Sides 2\n");
            f.Write("ShadingRate 0.5\n");
            f.Write("Quantize \"rgba\" 0 0 0 0\n");
            f.Write("Display \"{0}\" \"file\" \"rgba\"\n", tiffname);
            f.Write("Scale {0} {1} 1\n", 2.0 - 1.0 / size, 2.0 - 1.0 / size);
            f.Write("Translate 0 0 1\n");
            f.Write("WorldBegin\n");
            f.Write("\tSurface \"constant\"\n");
            for (i = 0; i < vertices.getCount(); i++)
            {
                vertices[i].tex.x -= 0.5;
                vertices[i].tex.y -= 0.5;
                vertices[i].label = i;
            }

            f.Write("\t PointsPolygons [\n");
            for (i = 0; i < faces.getCount(); i++)
            {
                f.Write("\t\t 3\n");
            }
            f.Write("\t ] [\n");
            for (i = 0; i < faces.getCount(); i++)
            {
                f.Write("\t\t {0} {1} {2}\n", faces[i].vertices[0].label,
                        faces[i].vertices[1].label,
                        faces[i].vertices[2].label);
            }
            f.Write("\t ]\n");
            f.Write("\t \"P\" [\n");
            for (i = 0; i < vertices.getCount(); i++)
            {
                f.Write("\t\t {0} {1} {2} \n",
                    vertices[i].tex.x,
                    vertices[i].tex.y,
                    vertices[i].tex.z);
            }
            f.Write("\t ]\n");
            f.Write("\t \"Cs\" [\n");
            for (i = 0; i < vertices.getCount(); i++)
            {
                f.Write("\t\t {0} {1} {2} \n",
                    vertices[i].normal.x,
                    vertices[i].normal.y,
                    vertices[i].normal.z);
            }
            f.Write("\t ]\n");
            f.Write("\tTranslate 0.0 0.0 0.5\n");
            for (i = 0; i < edges.getCount(); i++)
            {
                if (edges[i].boundary())
                {
                    v1 = edges[i].ends[0].tex * 0.5;
                    v2 = edges[i].ends[0].tex * 1.5;
                    v3 = edges[i].ends[1].tex * 1.5;
                    v4 = edges[i].ends[1].tex * 0.5;
                    f.Write("\t Polygon \"Cs\" [ ");
                    f.Write("{0} {1} {2} ",
                        edges[i].ends[0].normal.x,
                        edges[i].ends[0].normal.y,
                        edges[i].ends[0].normal.z);
                    f.Write("{0} {1} {2} ",
                        edges[i].ends[0].normal.x,
                        edges[i].ends[0].normal.y,
                        edges[i].ends[0].normal.z);
                    f.Write("{0} {1} {2} ",
                        edges[i].ends[1].normal.x,
                        edges[i].ends[1].normal.y,
                        edges[i].ends[1].normal.z);
                    f.Write("{0} {1} {2} ",
                        edges[i].ends[1].normal.x,
                        edges[i].ends[1].normal.y,
                        edges[i].ends[1].normal.z);
                    f.Write("]\n\t\t \"P\" [ ");
                    f.Write("{0} {1} {2} ", v1.x, v1.y, v1.z);
                    f.Write("{0} {1} {2} ", v2.x, v2.y, v2.z);
                    f.Write("{0} {1} {2} ", v3.x, v3.y, v3.z);
                    f.Write("{0} {1} {2} ", v4.x, v4.y, v4.z);
                    f.Write("]\n");
                }
            }
            f.Write("WorldEnd\n");
            f.Close();

            for (i = 0; i < vertices.getCount(); i++)
            {
                vertices[i].tex.x += 0.5;
                vertices[i].tex.y += 0.5;
            }
        }

        public List<int> subset;
        public VertexSet vertices;
        public EdgeSet edges;
        public FaceSet faces;
        public Vector innerRadius = new Vector(1f, 1f, 1f);
        public Vector innerRadius0 = new Vector(1f, 1f, 1f);

        //Compute the genus of the mesh (assuming all vertices and edges are used)
        //And assuming there are no boundaries
        //see Computer Aided Geometric Design (Farin) pg. 395
        public int SimpleGenus()
        {
            //Euler-Poincare formula
            return (vertexCount() - edgeCount() + faceCount() - 2) / 2;
        }

        public int faceCount()
        {
            return faces.getCount();
        }

        public int edgeCount()
        {
            return edges.getCount();
        }

        public int vertexCount()
        {
            return vertices.getCount();
        }

        public void RemoveFace(Face f)
        {
            int i;
            for (i = 0; i < 3; i++)
            {
                f.vertices[i].remove(f);
                f.edges[i].remove(f);
            }
            faces.remove(f);
        }

        public void RemoveFace(int j)
        {
            int i;
            Face f = faces[j];
            for (i = 0; i < 3; i++)
            {
                f.vertices[i].remove(f);
                f.edges[i].remove(f);
            }
            faces.remove(j);
        }

        public int RemoveOrphanedEdges()
        {
            int i, n;
            Edge e;
            i = 0;
            n = 0;
            while (i < edges.getCount())
            {
                if (edges[i].faces.getCount() == 0)
                {
                    edges[i].ends[0].remove(edges[i]);
                    edges[i].ends[1].remove(edges[i]);
                    e = edges[i];
                    edges.remove(i);
                    n++;
                }
                else
                {
                    i++;
                }
            }

            return n;
        }

        public int RemoveOrphanedVertices()
        {
            int i, n;
            Vertex v;
            i = 0;
            n = 0;
            while (i < vertices.getCount())
            {
                if ((vertices[i].edges.getCount() == 0) && (vertices[i].faces.getCount() == 0))
                {
                    v = vertices[i];
                    vertices.remove(i);
                    n++;
                }
                else
                {
                    i++;
                }
            }

            return n;
        }

        private void CheckConnected(Edge e)
        {
            int i;
            Vertex v;

            //if (e==null) cerr << "Broken graph: edge" << endl;
            if (e.label == 1) return;
            e.label = 1;
            v = e.ends[0];
            v.label = 1;
            for (i = 0; i < v.edges.getCount(); i++)
            {
                CheckConnected(v.edges[i]);
            }
            v = e.ends[1];
            v.label = 1;
            for (i = 0; i < v.edges.getCount(); i++)
            {
                CheckConnected(v.edges[i]);
            }
        }

        private void CheckConnected(Face f)
        {
            int i;

            //if (f==null) cerr << "Broken graph: face" << endl;
            if (f == null) return;
            if (f.label == 1) return;
            f.label = 1;

            for (i = 0; i < 3; i++)
            {
                CheckConnected(f.edges[i].faces[0]);
                CheckConnected(f.edges[i].faces[1]);
            }
        }

        public Vector ComputeCentreFromMesh()
        {
            Vector centre = new Vector();
            for(int i=0; i<vertexCount(); i++)
            {
                centre += vertices[i].v;
            }
            return centre/vertexCount();
        }

        public Vector[] GetVertexPositions()
        {
            Vector[] pos = new Vector[vertexCount()];
            for (int i = 0; i < vertexCount(); i++){
                pos[i] = vertices[i].v;
            }
            return pos;
        }

        public string PrintVertices()
        {
            StringBuilder verticesString = new StringBuilder();
            for (int i=0; i<vertexCount(); i++)
            {
                verticesString.AppendLine(vertices[i].v.ToString());
            }
            return verticesString.ToString(); 
        }

        public string PrintTriangles()
        {
            StringBuilder trianglesString = new StringBuilder();
            for (int i = 0; i < faceCount(); i++)
            {
                trianglesString.AppendLine(faces[i].ToString());
            }
            return trianglesString.ToString();
        }

        public bool isPointInConvexMesh(Vector point)
        {
            int i = 0;
            bool value = true;

            Vector centre = ComputeCentreFromMesh();
            
            while (value && i < faceCount())
            {
                value = value && faces[i].ArePointsOnSameSide(point, centre);
                i++;
            }
            return value;
        }

        public void ComputeFaceNormalsOnSubset()
        {
            int i;
            for (i = 0; i < faces.getCount(); i++)
            {
                faces[i].normal = new Vector(0,0,0);
                bool isTriangleInSubset = subset.Contains(faces[i].vertices[0].pos) && subset.Contains(faces[i].vertices[1].pos) && subset.Contains(faces[i].vertices[2].pos);
                if (isTriangleInSubset)
                {
                    faces[i].ComputeNormal();
                }
            }
        }

        public void ComputeVertexNormalsOnSubset()
        {
            int i;
            ComputeFaceNormalsOnSubset();
            for (i = 0; i < vertices.getCount(); i++)
            {
                vertices[i].normal = new Vector(0,0,0);
            }
            for (i = 0; i < subset.Count; i++)
            {
                vertices[subset[i]].ComputeNormal();
            }
        }

        public Vector MeanCurvatureAtPoint(int i)
        {
            Vector Ki = new Vector();
            double Ai = 0;

            bool isIOrigin = true;
            Vertex xi = vertices[i];
            Vertex xj = new Vertex();
            Vertex xk = new Vertex();

            for (int j = 0; j<vertices[i].edges.getCount(); j++) {
                Edge xixj = vertices[i].edges[j];

                isIOrigin = xixj.ends[0] == xi ? true : false;
                xj = isIOrigin ? xixj.ends[1] : xixj.ends[0];
                
                Face T1 = xixj.faces[0];
                Face T2 = xixj.faces[1];

                Vertex A = T1.ThirdVertex(xi, xj);
                Vertex B = T2.ThirdVertex(xi, xj);

                double coTanAxixj = ((xi.v - A.v) * (xj.v - A.v)) / ((xi.v - A.v) ^ (xj.v - A.v)).norm();
                double coTanBxixj = ((xi.v - B.v) * (xj.v - B.v)) / ((xi.v - B.v) ^ (xj.v - B.v)).norm();

                //Console.WriteLine(coTanAxixj + ", " + coTanBxixj);

                Ki += (coTanAxixj + coTanBxixj) * (xi.v - xj.v);    
            }

            Console.WriteLine("Ki: " + Ki);

            // Compute Mixed area 
            for (int j = 0; j<vertices[i].faces.getCount(); j++)
            {
                if(vertices[i].faces[j].vertices[0] == xi)
                {
                    xj = vertices[i].faces[j].vertices[1];
                    xk = vertices[i].faces[j].vertices[2];
                }
                else if (vertices[i].faces[j].vertices[1] == xi)
                {
                    xj = vertices[i].faces[j].vertices[2];
                    xk = vertices[i].faces[j].vertices[0];
                }
                else if (vertices[i].faces[j].vertices[2] == xi)
                {
                    xj = vertices[i].faces[j].vertices[0];
                    xk = vertices[i].faces[j].vertices[1];
                }

                Vector xixj = xj.v - xi.v;
                Vector xjxk = xk.v - xj.v;
                Vector xkxi = xi.v - xk.v;

                if (xixj.sqrNorm() > xjxk.sqrNorm() + xkxi.sqrNorm())
                {
                    Ai += 0.25 * (xixj ^ xjxk).norm()/2;
                }
                else if (xjxk.sqrNorm() > xixj.sqrNorm() + xkxi.sqrNorm())
                {
                    Ai += 0.5 * (xixj ^ xjxk).norm()/2;
                }
                else if (xkxi.sqrNorm() > xixj.sqrNorm() + xjxk.sqrNorm())
                {
                    Ai += 0.25 * (xixj ^ xjxk).norm()/2;
                }
                else
                {
                    // Compute Voronoi Area
                    Vector IJNormalUnit = xixj^(xixj^xkxi);
                    IJNormalUnit.normalize();

                    Vector IKNormalUnit = xkxi^(xixj^xkxi);
                    IKNormalUnit.normalize();

                    //Console.WriteLine("Normal Units: " + IJNormalUnit + ", " + IKNormalUnit);

                    double S = (xixj.norm() + xjxk.norm() + xkxi.norm()) / 2f;
                    double K = Math.Sqrt(S * (S - xkxi.norm()) * (S - xjxk.norm()) * (S - xkxi.norm()));
                    
                    double radius = xkxi.norm() * xjxk.norm() * xkxi.norm() / (4f * K);
                    

                    Vector side1 = Math.Sqrt(radius * radius - xixj.sqrNorm() / 4f) * IJNormalUnit;
                    Vector side2 = Math.Sqrt(radius * radius - xkxi.sqrNorm() / 4f) * IKNormalUnit;

                    Console.WriteLine("Edges: " + xixj + ", " + xkxi);

                    Console.WriteLine("Sides 1&2: " + side1 + ", " + side2);

                    Ai += (xixj.norm() * side1.norm() + xkxi.norm() * side2.norm()) / 4f;

                    
                }
            }
            //Console.WriteLine("Ai: " + Ai);

            Ki /= (2 * Ai);

            return Ki;
        }

        /*
        public Vector MeanCurvatureAtPoint(int i)
        {
            Vector xkxi = new Vector();
            Vector xixj, Axi, Axj, Bxi, Bxj;
            float Ai = 0, coTanAxixj = 0, coTanBxixj = 0;

            for(int j=0; j<vertices[i].edges.getCount(); j++){
                xixj = vertices[i].edges[j].ends[0].v - vertices[i].edges[j].ends[1].v;
                
                //for(int k=0; k<vertices[i].edges[j].faces.getCount(); k++){
                    Axi = vertices[i].edges[j].ends[0].v - vertices[i].edges[j].faces[0].vertices[2].v;
                    Axj = vertices[i].edges[j].ends[1].v - vertices[i].edges[j].faces[0].vertices[2].v;
                    Bxi = vertices[i].edges[j].ends[0].v - vertices[i].edges[j].faces[1].vertices[2].v;
                    Bxj = vertices[i].edges[j].ends[1].v - vertices[i].edges[j].faces[1].vertices[2].v;
                //}

                coTanAxixj = Vector.Dot(Axi, Axj) / Vector.Cross(Axi, Axj).norm();
                coTanBxixj = Vector.Dot(Bxi, Bxj) / Vector.Cross(Bxi, Bxj).norm();

                xkxi += (coTanAxixj + coTanBxixj) * xixj;
                Ai += VoronoiAreaOfTriangle(i, j, thirdPoints[0]);
                Ai += VoronoiAreaOfTriangle(i, j, thirdPoints[1]);
            }
            xkxi /= Ai;

            for(int t=0; t<vertices[i].faces.getCount(); t++){
                
                Vector xixj = vertices[i].faces[t].vertices[1] - vertices[i].faces[t].vertices[0];
                Vector xjxk = vertices[i].faces[t].vertices[2] - vertices[i].faces[t].vertices[1];
                Vector xkxi = vertices[i].faces[t].vertices[0] - vertices[i].faces[t].vertices[2];
            }
        
            return xkxi;
        }
        */

        /*
        public void MeanCurvatureAtPoint(int i)
        {
            int[] thirdPoints;
            Vector xixj, Axi, Axj, Bxi, Bxj;
            int A, B;

            float Ai = 0, coTanAxixj = 0, coTanBxixj = 0;
            Vector xkxi = new Vector();
            for (int k = 0; k < this.particleNeighbourhoods[i].Count; k++)
            {
                int j = this.particleNeighbourhoods[i][k];
                thirdPoints = this.ThirdPointOfTriangle(i, j);
                A = thirdPoints[0];
                B = thirdPoints[1];
                xixj = this.membraneSCEs[i].GetPosition() - this.membraneSCEs[j].GetPosition();
                Axi = this.membraneSCEs[i].GetPosition() - this.membraneSCEs[A].GetPosition();
                Axj = this.membraneSCEs[j].GetPosition() - this.membraneSCEs[A].GetPosition();
                Bxi = this.membraneSCEs[i].GetPosition() - this.membraneSCEs[B].GetPosition();
                Bxj = this.membraneSCEs[j].GetPosition() - this.membraneSCEs[B].GetPosition();

                coTanAxixj = Vector.Dot(Axi, Axj) / Vector.Cross(Axi, Axj).norm();
                coTanBxixj = Vector.Dot(Bxi, Bxj) / Vector.Cross(Bxi, Bxj).norm();

                xkxi += (coTanAxixj + coTanBxixj) * xixj;
                Ai += VoronoiAreaOfTriangle(i, j, thirdPoints[0]);
                Ai += VoronoiAreaOfTriangle(i, j, thirdPoints[1]);
            }
            xkxi /= Ai;


            return xkxi;
        }

        public void GaussianCurvature(int i)
        {

        }
        */

        public void Default2RoundedCube(float r)
        {
            Vector centre = new Vector();
            Vector pointO = vertices[0].v;
            Vector[] targetPos = new Vector[vertexCount()];

            for (int i = 0; i < vertexCount(); i++)
            {
                targetPos[i] = CubizePoint2(vertices[i].v);
                targetPos[i].x = r * targetPos[i].x;
                targetPos[i].z = r * targetPos[i].z;

                vertices[i].v = targetPos[i] + centre;
            }

            innerRadius = new Vector(r, 2 * pointO.norm(), r);
        }

        public Vector CubizePoint2(Vector pos)
        {
            Vector position = pos;

            double x, y, z;
            x = position.x;
            y = position.y;
            z = position.z;

            double fx, fy, fz;
            fx = Math.Abs(x);
            fy = Math.Abs(y);
            fz = Math.Abs(z);

            const float inverseSqrt2 = 0.70710676908493042f;

            if (fy >= fx && fy >= fz)
            {
                double a2 = x * x * 2f;
                double b2 = z * z * 2f;
                double inner = -a2 + b2 - 3;
                double innersqrt = -Math.Sqrt((inner * inner) - 12f * a2);

                if (Math.Abs(x) < Mathf.epsilon)
                {
                    position.x = 0f;
                }
                else
                {
                    position.x = Math.Sqrt(innersqrt + a2 - b2 + 3f) * inverseSqrt2;
                }

                if (z == 0.0 || z == -0.0)
                {
                    position.z = 0f;
                }
                else
                {
                    position.z = Math.Sqrt(innersqrt - a2 + b2 + 3f) * inverseSqrt2;
                }

                if (position.x > 1.0) position.x = 1f;
                if (position.z > 1.0) position.z = 1f;

                if (x < 0) position.x = -position.x;
                if (z < 0) position.z = -position.z;

                if (y > 0)
                {
                    // top face
                    position.y = 1f;
                }
                else
                {
                    // bottom face
                    position.y = -1f;
                }
            }
            else if (fx >= fy && fx >= fz)
            {
                double a2 = y * y * 2f;
                double b2 = z * z * 2f;
                double inner = -a2 + b2 - 3;
                double innersqrt = -Math.Sqrt((inner * inner) - 12f * a2);

                if (Math.Abs(y) < Mathf.epsilon)
                {
                    position.y = 0f;
                }
                else
                {
                    position.y = Math.Sqrt(innersqrt + a2 - b2 + 3f) * inverseSqrt2;
                }

                if (z == 0.0 || z == -0.0)
                {
                    position.z = 0f;
                }
                else
                {
                    position.z = Math.Sqrt(innersqrt - a2 + b2 + 3f) * inverseSqrt2;
                }

                if (position.y > 1.0) position.y = 1f;
                if (position.z > 1.0) position.z = 1f;

                if (y < 0) position.y = -position.y;
                if (z < 0) position.z = -position.z;

                if (x > 0)
                {
                    // right face
                    position.x = 1f;
                }
                else
                {
                    // left face
                    position.x = -1f;
                }
            }
            else
            {
                double a2 = x * x * 2f;
                double b2 = y * y * 2f;
                double inner = -a2 + b2 - 3;
                double innersqrt = -Math.Sqrt((inner * inner) - 12f * a2);

                if (x == 0.0 || x == -0.0)
                {
                    position.x = 0f;
                }
                else
                {
                    position.x = Math.Sqrt(innersqrt + a2 - b2 + 3f) * inverseSqrt2;
                }

                if (y == 0.0 || y == -0.0)
                {
                    position.y = 0f;
                }
                else
                {
                    position.y = Math.Sqrt(innersqrt - a2 + b2 + 3f) * inverseSqrt2;
                }

                if (position.x > 1f) position.x = 1f;
                if (position.y > 1f) position.y = 1f;

                if (x < 0) position.x = -position.x;
                if (y < 0) position.y = -position.y;

                if (z > 0)
                {
                    // front face
                    position.z = 1f;
                }
                else
                {
                    // back face
                    position.z = -1f;
                }
            }

            return position;
        }

        public void ScaleMesh(Vector r)
        {
            Vector centre = ComputeCentreFromMesh();
            for(int i=0; i<vertexCount(); i++)
            {
                vertices[i].v = (vertices[i].v - centre) % r + centre;
            }
        }
    }
}
