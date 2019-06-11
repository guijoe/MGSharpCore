using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using MGSharp.MIConvexHull;

namespace MGSharp.Core.GeometricPrimitives
{
    public class Image2D
    {
        public Bitmap image;
        List<Vector> points;

        public Vector[] positions;
        public int nbOfPositions = 150;
        Vector[] velocities;
        public Bounds2D bounds;
        Color[] oldColors;
        double R;


        public Image2D(string folderName, string imageName)
        {
            string path = folderName + "/" + imageName;
            image = new Bitmap(path);
            positions = new Vector[nbOfPositions];

            for(int i=0; i<positions.Length; i++)
            {
                positions[i] = new Vector(0, 0, 0);
            }
        }
        
        public Bounds2D FindImageBounds()
        {
            points = new List<Vector>();

            for (int x = 0; x < image.Width; x++)
            {
                for (int y = 0; y < image.Height; y++)
                {
                    Color color = image.GetPixel(x, y);
                    if (color.R > 127)
                    {
                        points.Add(new Vector(x, y));
                    }
                }
            }

            Vector min = new Vector(points[0]);
            Vector max = new Vector(points[0]);
            Vector centre = new Vector(points[0]);
            for (int i = 1; i < points.Count; i++)
            {
                if (points[i].x < min.x)
                {
                    min.x = points[i].x;
                }
                if (points[i].y < min.y)
                {
                    min.y = points[i].y;
                }

                if (points[i].x > max.x)
                {
                    max.x = points[i].x;
                }
                if (points[i].y > max.y)
                {
                    max.y = points[i].y;
                }
                centre += points[i];
            }
            centre /= points.Count;

            bounds = new Bounds2D();
            bounds.centre = (max+min)/2;
            bounds.size = max - min + new Vector(1,1);
            bounds.extents = bounds.size / 2;
            bounds.max = max;
            bounds.min = min;

            //Bounds bds = new Bounds(bounds);
            return bounds;
        }

        public void FindCurvature()
        {
            positions = new Vector[nbOfPositions];
            oldColors = new Color[nbOfPositions];

            double theta = 2 * Math.PI / positions.Length;
            //R = bounds.extents.norm()/4;
            //Console.WriteLine("R1: " + R);
            R = bounds.extents.x < bounds.extents.y ? bounds.extents.x : bounds.extents.y;
            R /= 4;
            //Console.WriteLine("R2: " + R);

            for (int i=0; i<positions.Length; i++)
            {
                positions[i] = new Vector();
                positions[i].x = R * Math.Cos(i * theta) + bounds.centre.x;
                positions[i].y = R * Math.Sin(i * theta) + bounds.centre.y;
                //positions[i].y = bounds.extents.y * R * Math.Sin(i * theta) / bounds.extents.x + bounds.centre.y;
                //Console.WriteLine(positions[i]);

            }

            int resolution = 50;
            int NConvolutions = 50;
            double KPM = 100;
            double deltaK = .1;

            resolution = (int) Math.Ceiling(bounds.size.x) +1;
            Console.WriteLine(bounds.size);

            //bounds.DiscretiseSpace(resolution, Kronecker);
            //bounds.DiscretiseSpace(new Vector(bounds.size.y, bounds.size.x), Kronecker);
            bounds.DiscretiseSpace(bounds.size, Kronecker);

            bounds.GaussianConvolution(NConvolutions);
            bounds.ComputeVectorFields(bounds.size, KPM, deltaK);

            //*
            for (int i = 0; i < bounds.size.x; i++)
            {
                for (int j = 0; j < bounds.size.y; j++)
                {
                    //Console.WriteLine(bounds.vPHI[i, j]);
                    int value = (int)(bounds.vPHI[i, j] * 255);
                    Color c = Color.FromArgb(255, value, value, value);
                    if(bounds.X[i,j].x<=image.Width && bounds.X[i, j].x >= 0 && bounds.X[i,j].y <= image.Height && bounds.X[i, j].y >= 0)
                    {
                        //image.SetPixel((int)bounds.X[i, j].x, (int)bounds.X[i, j].y, c);
                        //image.SetPixel((int)bounds.X[i, j].x, (int)bounds.X[i, j].y, Color.White);
                        //image.SetPixel((int)bounds.X[i, j].x + 1, (int)bounds.X[i, j].y + 1, Color.Red);
                    }
                }
            }

            for (int i = 0; i < positions.Length; i++)
            {
                oldColors[i] = image.GetPixel((int)positions[i].x, (int)positions[i].y);
                //image.SetPixel((int)positions[i].x, (int)positions[i].y, Color.Red);
            }
            //*/
        }

        public void WrapImage(Bounds2D bounds)
        {
            for (int i = 0; i < positions.Length; i++)
            {
                Vector IJK = bounds.DichotomicSearch(positions[i]);
                IJK.z = 0;

                //if (IJK.x >= 0 && IJK.x < bounds.resolution && IJK.y >= 0 && IJK.y < bounds.resolution)
                if (IJK.x >= 0 && IJK.x < bounds.resolution.x && IJK.y >= 0 && IJK.y < bounds.resolution.y)
                {
                    Vector normal = ComputeNormal(i);

                    Vector v = -bounds.gradPeronaMalik[(int)IJK.x, (int)IJK.y];
                    if (v.norm() > 0)
                    {
                        v *= .1 / v.norm();
                    }

                    Vector velocity = new Vector();
                    if (bounds.vPHI[(int)IJK.x, (int)IJK.y] > .5)
                    {
                        velocity = (.5) * bounds.vPHI[(int)IJK.x, (int)IJK.y] * normal + .5 * (bounds.gradPeronaMalik[(int)IJK.x, (int)IJK.y] * normal) * normal - .1 * normal;
                    }
                    else
                    {
                        //Console.WriteLine("Hello");
                        //velocity = (.5) * bounds.vPHI[(int)IJK.x, (int)IJK.y] * normal + .5 * (bounds.gradPeronaMalik[(int)IJK.x, (int)IJK.y] * normal) * normal + .1 * normal;
                        //Console.WriteLine(bounds.vPHI[(int)IJK.x, (int)IJK.y]);
                    }

                    
                    //velocity = (.3) * bounds.vPHI[(int)IJK.x, (int)IJK.y] * normal - .7 * v - .25 * normal;

                    //Console.WriteLine(velocity);
                    Vector normalisedVelocity = velocity / velocity.norm();
                    Vector direction = new Vector();

                    
                    positions[i] += velocity;
                    
                    //image.SetPixel((int)positions[i].x, (int)positions[i].y, oldColors[i]);
                    //oldColors[i] = image.GetPixel((int)positions[i].x, (int)positions[i].y);
                    //image.SetPixel((int)positions[i].x, (int)positions[i].y, Color.Red);
                    //*/
                }
            }
        }

        public Vector ComputeNormal(int i)
        {
            Vector normal;
            int j, k;

            if (i == 0)
            {
                j = positions.Length - 1;
                k = i+1;
            }else if (i == positions.Length - 1)
            {
                j = i - 1;
                k = 0;
            }
            else
            {
                j = i - 1;
                k = i + 1;
            }
            Vector tangent = positions[k] - positions[j];
            if (tangent.y == 0) {
                normal = new Vector(0, 1);
            }
            else
            {
                normal = new Vector(1, -tangent.x / tangent.y);
            }
            
            if(normal.x*tangent.y-normal.y*tangent.x < 0)
            {
                normal = -normal;
            }
            //normal.x *= bounds.extents.x;
            //normal.y *= bounds.extents.y;

            normal.normalize();
            //normal *= tangent.norm();
            return normal;
        }

        public void FindConvexHull()
        {
            var vertices = new MIVertex2D[points.Count];

            for (var i = 0; i < vertices.Length; i++) { 
                vertices[i] = new MIVertex2D(points[i].x, points[i].y);
            }

            MIVertex2D[] convexHull = (MIVertex2D[])ConvexHull.Create(vertices).Points;

            for(int i=0; i<convexHull.Length; i++)
            {
                image.SetPixel((int) convexHull[i].Position[0], (int)convexHull[i].Position[1], Color.Red);
            }
        }

        public void Save(string folderName, string imageName, int number)
        {
            imageName = imageName.Insert(imageName.LastIndexOf('.'), "_WithCurvature" + number);
            image.Save(folderName + "/" + imageName);
        }

        public void Save(string folderName, string imageName, string number)
        {
            imageName = imageName.Insert(imageName.LastIndexOf('.'), "_WithCurvature" + number);
            image.Save(folderName + "/" + imageName);
        }

        public void Save(string folderName, string imageName, string number, Rectangle section)
        {
            imageName = imageName.Insert(imageName.LastIndexOf('.'), "_WithCurvature" + number);
            
            //Bitmap croppedImage = CropImage(image, section);
            Bitmap croppedImage = image.Clone(section, image.PixelFormat);
            croppedImage.Save(folderName + "/" + imageName);
        }

        public Bitmap CropImage(Bitmap source, Rectangle section)
        {
            // An empty bitmap which will hold the cropped image
            Bitmap bmp = new Bitmap(section.Width, section.Height);
            
            Graphics g = Graphics.FromImage(bmp);

            // Draw the given area (section) of the source image
            // at location 0,0 on the empty bitmap (bmp)
            g.DrawImage(source, 0, 0, section, GraphicsUnit.Pixel);
            g.DrawImage(source, new Point());

            return bmp;
        }
        
        public double Kronecker(Vector u, Bounds2D bounds)
        {
            double val = 0;
            if(u.x>=0 && u.x < image.Width && u.y >= 0 && u.y < image.Height) { 
                Color color = image.GetPixel((int)u.x, (int)u.y);
                if (color.R > 127)
                {
                    val = 1;
                    //image.SetPixel((int)u.x, (int)u.y, Color.Red);
                }
            }
            return val;
        }

        public void FindSymetryAxisAndRotate(string folderName, string imageName)
        {
            List<Vector> points = new List<Vector>();
            double X_hat = 0;
            double Y_hat = 0;
            double sumXSquare = 0;
            double sumYSquare = 0;
            double sumXY = 0;
            double numPart = 0;

            double a = 0, b = 0, b_num = 0, b_denum = 0;
            int count = 0;


            for (int x = 0; x < image.Width; x++)
            {
                for (int y = 0; y < image.Height; y++)
                {
                    Color color = image.GetPixel(x, y);
                    if (color.R > 127)
                    {

                        X_hat += x;
                        Y_hat += y;
                        sumXSquare += x * x;
                        sumYSquare += y * y;
                        sumXY += x * y;
                        points.Add(new Vector(x, y));
                        count++;
                    }
                }
            }
            X_hat /= count;
            Y_hat /= count;

            for (int i = 0; i < points.Count; i++)
            {
                b_num += (points[i].x - X_hat) * (points[i].y - Y_hat);
                b_denum += (points[i].x - X_hat) * (points[i].x - X_hat);
            }
            b = b_num / b_denum;
            a = Y_hat - b * X_hat;

            Console.WriteLine(a + ", " + b);
            Color red = Color.Red;
            for (int x = 0; x < image.Width; x++)
            {
                int y = (int)Math.Round(b * x + a);
                if (y < image.Height && y > 0)
                {
                    //image.SetPixel(x, y, red);
                }
            }

            b_denum = 2 * (count * sumXY - count * count * X_hat * Y_hat);
            numPart = count * count * (X_hat * X_hat - Y_hat * Y_hat) - count * (sumXSquare - sumYSquare);
            b_num = numPart + Math.Sqrt(numPart * numPart + b_denum * b_denum);
            b = b_num / b_denum;
            a = Y_hat - b * X_hat;

            Console.WriteLine(a + ", " + b);
            bool originFound = false;
            Vector O = new Vector();
            for (int x = 0; x < image.Width; x++)
            {
                int y = (int)Math.Round(b * x + a);
                if (y < image.Height && y > 0)
                {
                    image.SetPixel(x, y, red);
                    O = new Vector(x, y);
                    if (!originFound && points.Contains(O))
                    {
                        originFound = true;
                    }
                }
            }

            double theta = Math.Atan2(b, 1);
            for (int i = 0; i < points.Count; i++)
            {
                Vector OP = points[i] - O;
                int x = (int)(Math.Round(OP.x * Math.Cos(theta) - OP.y * Math.Sin(theta)) + O.x);
                int y = (int)(Math.Round(OP.x * Math.Sin(theta) + OP.y * Math.Cos(theta)) + O.y);

                image.SetPixel(x, y, Color.Yellow);
            }

            Save(folderName, imageName + "_SymetryAxis.png", 0);

        }
    }
}
