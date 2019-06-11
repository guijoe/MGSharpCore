using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MGSharp.Core.GeometricPrimitives
{
    public class Set
    {
        public Set()
        {
            n = 0;
            size = 0;
            set = null;
        }

        private void push_back(object x)
        {
            int ns, i;
            object[] newset;
            //Console.WriteLine("n:"+n);
            if (n < size) set[n] = x;
            else
            {
                if (size == 0) ns = 1;
                //grow by ten percent, but at least 1
                else ns = size * 110 / 100 + 1;
                newset = new object[ns];
                for (i = 0; i < size; i++) newset[i] = set[i];
                size = ns;
                set = newset;
                set[n] = x;
            }
        }

        private void pop_back()
        {
            
        }

        public void clear()
        {
            n = 0;
            size = 0;
            set = null;
        }

        public int add(object x)
        {
            int i;
            for (i = 0; i < n; i++)
            {
                if (set[i].Equals(x)) return i;
            }
            push_back(x);
            n++;
            return n - 1;
        }

        public int addWithoutCheck(object x)
        {
            push_back(x);
            n++;
            return n - 1;
        }

        public void remove(object x)
        {
            int i;
            for (i = 0; i < n; i++)
            {
                if (set[i].Equals(x))
                {
                    set[i] = set[n - 1];
                    pop_back();
                    n--;
                    return;
                }
            }
        }

        public void removeAll(int[] all)
        {
            int[] perm = new int[n];
            for(int i=0; i<n; i++)
            {
                perm[i] = i;
            }

            for(int i=0; i<all.Length; i++)
            {
                if (all[i] < 0) { }
                else if (all[i] >= n) { }
                else{
                    set[perm[all[i]]] = set[perm[n - 1]];
                    perm[n - 1] = all[i];
                    pop_back();
                    n--;
                }
            }
        }

        public void remove(int i)
        {
            if (i < 0) return;
            if (i >= n) return;
            set[i] = set[n - 1];
            pop_back();
            n--;
        }

        public int index(object x)
        {
            int i;
            for (i = 0; i < n; i++)
            {
                if (set[i].Equals(x)) return i;
            }
            return -1;
        }

        public bool contains(object x)
        {
            int i = index(x);
            if (i != -1) return true;
            else return false;
        }

        public object match(object x)
        {
            int i = index(x);
            if (i != -1) return set[i];
            else return null;
        }

        public int getCount()
        {
            return n;
        }

        public object this[int i]
        {
            get { return set[i]; }
            set { if (i >= n) n = i + 1; set[i] = value; }
        }

        public void setIndex(int i, object x)
        {
            set[i] = x;
        }

        //union operator
        public static Set operator +(Set s1, Set s2)
        {
            int i;
            Set s = new Set();

            for (i = 0; i < s1.n; i++)
            {
                s.add(s1.set[i]);
            }
            for (i = 0; i < s2.n; i++)
            {
                s.add(s2.set[i]);
            }
            return s;
        }
        //reserve a minimum of so many items
        public void reserve(int s)
        {
            int i;
            object[] newset;
            if (s <= size) return;
            else
            {
                if (s == 0) s = 1;
                newset = new object[s];
                for (i = 0; i < size; i++) newset[i] = set[i];
                size = s;
                set = newset;
            }
        }

        protected object[] set;
        public int n;
        public int size;
    }
}
