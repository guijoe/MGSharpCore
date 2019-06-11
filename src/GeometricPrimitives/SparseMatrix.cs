using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MGSharp.Core.LinearAlgebra
{
    public class SparseMatrix
    {
        public SparseRow[] rows;
        public int maxrows;
        public int maxcols;
        public SparseMatrix()
        {
            maxrows = maxcols = 0;
            rows = null;
        }
        public SparseMatrix(int rows, int cols)
        {
            maxrows = rows;
            maxcols = cols;
            this.rows = new SparseRow[maxrows];
        }

        public SparseMatrix(SparseMatrix M)
        {
            int i;
            Resize(M.maxrows, M.maxcols);
            for (i = 0; i < maxrows; i++)
            {
                if (M.rows[i] != null)
                    rows[i] = new SparseRow(M.rows[i]);
            }
        }

        public double this[int row, int col]
        {
            get { return rows[row].GetValue(col); }
            set
            {
                if (rows[row] == null) rows[row] = new SparseRow(row);
                rows[row].SetValue(col, value);
            }
        }
        public void Resize(int rows, int cols)
        {
            if ((this.maxrows == rows) && (this.maxcols == cols)) return;
            if ((rows < 0) || (cols < 0))
            {
                maxrows = 0;
                maxcols = 0;
                this.rows = null;
            }
            else
            {
                maxrows = rows;
                maxcols = cols;
                this.rows = new SparseRow[maxrows];
            }
        }

        public bool isZero(int row, int col)
        {
            if (rows[row] == null) return true;
            if (rows[row].Exists(col)) return false;
            return true;
        }

        //note: rows should be the same as columns
        //return value <=0 implies failed.
        //     0  = failed to converge
        //     <0 = underdetermined system of equations 	
        public int Jacobi(double[] b, int maxiter, double eps)
        {
            return Jacobi(b, maxiter, eps, null);
        }
        public int Jacobi(double[] b, int maxiter, double eps, double[] guess)
        {
            double[] e, bnew, cb, ob, temp;
            double diff;
            double Mii = 0.0;
            int i, j, col;
            bool converged = false;
            if (guess == null)
                if (!Prepare(0)) return -1;

            e = new double[maxrows];
            bnew = new double[maxrows];
            for (i = 0; i < maxrows; i++) e[i] = b[i];
            if (guess != null)
                for (i = 0; i < maxrows; i++) b[i] = guess[i];

            cb = b;
            ob = bnew;

            while ((!converged) && (maxiter > 0))
            {
                temp = cb; cb = ob; ob = temp;
                diff = 0.0;
                for (i = 0; i < maxrows; i++)
                {
                    cb[i] = e[i];
                    if (rows[i] != null)
                        for (j = 0; j < rows[i].n; j++)
                        {
                            col = rows[i].data[j].col;
                            if (i != col)
                            {
                                cb[i] -= rows[i].data[j].x * ob[col];
                            }
                            else Mii = rows[i].data[j].x;
                        }
                    cb[i] /= Mii;
                    diff += (cb[i] - ob[i]) * (cb[i] - ob[i]);
                }
                maxiter--;
                if (diff < eps) converged = true;
            }

            if (cb != b)
                for (i = 0; i < maxrows; i++) b[i] = cb[i];

            return maxiter;
        }

        //note: rows should be the same as columns
        //return value <=0 implies failed.
        //     0  = failed to converge
        //     <0 = underdetermined system of equations 	
        public int GaussSeidel(double[] b, int maxiter, double eps)
        {
            return GaussSeidel(b, maxiter, eps, null);
        }
        public int GaussSeidel(double[] b, int maxiter, double eps, double[] guess)
        {
            double[] e;
            double diff, bnew;
            double Mii = 0.0;
            int i, j, col;
            bool converged = false;
            if (guess == null)
                if (!Prepare(0)) return -1;

            e = new double[maxrows];
            for (i = 0; i < maxrows; i++) e[i] = b[i];
            if (guess != null)
                for (i = 0; i < maxrows; i++) b[i] = guess[i];

            while ((!converged) && (maxiter > 0))
            {
                diff = 0.0;
                for (i = 0; i < maxrows; i++)
                {
                    bnew = e[i];
                    if (rows[i] != null)
                        for (j = 0; j < rows[i].n; j++)
                        {
                            col = rows[i].data[j].col;
                            if (i != col)
                            {
                                bnew -= rows[i].data[j].x * b[col];
                            }
                            else Mii = rows[i].data[j].x;
                        }
                    bnew /= Mii;
                    diff += (bnew - b[i]) * (bnew - b[i]);
                    b[i] = bnew;
                }
                maxiter--;
                if (diff < eps) converged = true;
            }
            return maxiter;
        }

        //note: rows should be the same as columns
        //return value <=0 implies failed.
        //     0  = failed to converge
        //     <0 = underdetermined system of equations 
        //		A is the areas provided by radiosity
        //	Only suitable for radiosity!
        public int Shooting(double[] e, double[] a, int maxiter, double eps)
        {
            double[] b, db;
            double diff, bnew;
            int i, j, k, col;
            bool converged = false;
            if (!Prepare(0)) return -1;

            b = new double[maxrows];
            db = new double[maxrows];
            for (i = 0; i < maxrows; i++) b[i] = 0;
            for (i = 0; i < maxrows; i++) db[i] = e[i];

            while ((!converged) && (maxiter > 0))
            {
                diff = 0.0;
                j = 0;
                for (i = 1; i < maxrows; i++)
                {
                    if (db[i] * a[i] > db[j] * a[j]) j = i;
                }
                b[j] += db[j];
                for (i = 0; i < maxrows; i++)
                {
                    bnew = db[i];
                    if (rows[i] != null)
                        for (k = 0; k < rows[i].n; k++)
                        {
                            col = rows[i].data[k].col;
                            if (k == col)
                            {
                                bnew += rows[i].data[k].x * b[j];
                            }
                        }
                    diff += (bnew - db[i]) * (bnew - db[i]);
                    db[i] = bnew;
                }
                db[j] = 0;
                maxiter--;
                if (diff < eps) converged = true;
            }
            return maxiter;
        }
        //make sure that the elements on the diagonal are nonzero
        //if this can't be done, then Guass-Siedel and Jacobi
        // iteration will not work
        // very quick check to try ordering, if it fails then it reverts
        // to the backtracking version
        private bool Prepare(int start)
        {
            int i, j, k;
            double max;
            bool found;
            SparseElement col;
            SparseRow trow;

            if (start == 0)
            {
                for (i = 0; i < maxrows; i++)
                {
                    if (rows[i] == null) return false;
                }
            }

            for (i = 0; i < maxrows; i++)
            {
                found = false;
                max = 0.0;
                for (j = i; (!found) && (j < maxrows); j++)
                {
                    if (rows[j] != null)
                        for (k = 0; k < rows[j].n; k++)
                        {
                            col = rows[j].data[k];
                            if (col.col == i)
                            {
                                if ((Math.Abs(col.x) > 1e-8) && (Math.Abs(col.x) > max))
                                {
                                    //found a later row
                                    //swap the rows	
                                    //we speculate that this will work
                                    trow = rows[i];
                                    rows[i] = rows[j];
                                    rows[j] = trow;
                                    found = true;
                                    max = Math.Abs(col.x);
                                    break;
                                }
                            }
                        }
                }
                if (!found) return PrepareBacktrack(start);
            }
            return true;
        }

        //make sure that the elements on the diagonal are nonzero
        //if this can't be done, then Guass-Siedel and Jacobi
        // iteration will not work
        // this version backtracks to test all options
        private bool PrepareBacktrack(int start)
        {
            int i, j, k;
            bool found;
            SparseElement col;
            SparseRow trow;

            if (start == 0)
            {
                for (i = 0; i < maxrows; i++)
                {
                    if (rows[i] == null) return false;
                }
            }

            //found a solution
            if (start >= maxrows) return true;

            //see if we can get a value in the diagonal for this row
            // if so, we are done, otherwise backtrack and try 
            // another option
            i = start;
            found = false;
            //include row i, so that we include the current row
            for (j = i; (!found) && (j < maxrows); j++)
            {
                if (rows[j] != null)
                    for (k = 0; k < rows[j].n; k++)
                    {
                        col = rows[j].data[k];
                        if (col.col == i)
                        {
                            if (Math.Abs(col.x) > 1e-8)
                            {
                                //found a later row
                                //swap the rows	
                                //we speculate that this will work
                                trow = rows[i];
                                rows[i] = rows[j];
                                rows[j] = trow;
                                found = (PrepareBacktrack(start + 1));

                                if (!found)
                                {
                                    //this did not work, undo the row swap
                                    trow = rows[i];
                                    rows[i] = rows[j];
                                    rows[j] = trow;
                                }
                            }
                        }
                    }
            }
            return found;
        }

        public object Clone()
        {
            Object o = new SparseMatrix(this);
            return o;
        }
    }

    public class SparseElement
    {
        public SparseRow row;
        public int col;
        public double x; //the value of the entry
        public SparseElement(SparseRow parent, int col, double x)
        {
            row = parent;
            this.col = col;
            this.x = x;
        }
        public SparseElement(SparseRow parent, SparseElement e)
        {
            row = parent;
            col = e.col;
            x = e.x;
        }
    }

    public class SparseRow
    {
        public int row;
        public SparseElement[] data;
        public int n;
        public SparseRow(int row)
        {
            this.row = row;
            data = null;
            n = 0;
        }

        public SparseRow(SparseRow r)
        {
            int i;
            this.row = r.row;
            this.data = null;
            n = 0;
            for (i = 0; i < r.n; i++)
            {
                SetValue(r.data[i].col, r.data[i].x);
            }
        }

        public void SetValue(int col, double x)
        {
            int i;
            SparseElement[] newdata = null;
            for (i = 0; i < n; i++)
            {
                if (data[i].col == col)
                {
                    data[i].x = x;
                    return;
                }
            }
            //not found, add the element
            newdata = new SparseElement[n + 1];
            for (i = 0; i < n; i++)
            {
                newdata[i] = data[i];
            }
            data = newdata;
            data[n] = new SparseElement(this, col, x);
            n++;
        }

        public double GetValue(int col)
        {
            int i;
            for (i = 0; i < n; i++)
            {
                if (data[i].col == col)
                {
                    return data[i].x;
                }
            }
            //not found
            return 0;
        }

        public bool Exists(int col)
        {
            int i;
            for (i = 0; i < n; i++)
            {
                if (data[i].col == col)
                {
                    return true;
                }
            }
            //not found
            return false;
        }
    }
}
