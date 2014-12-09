using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Experiment
{
    class Program
    {
        static void Main(string[] args)
        {
            //double[] x = new double[32];
            //int i;

            //// Makes a fancy cubic signal
            //for (i=0;i<32;i++) x[i]=5+i+0.4*i*i-0.02*i*i*i;
  
            //// Prints original sigal x
            //printf("Original signal:\n");
            //for (i=0;i<32;i++) printf("x{0}={1}",i.ToString(),x[i].ToString());
            //printf("\n");

            //// Do the forward 9/7 transform
            //Dwt97.Transform(ref x);
  
            //// Prints the wavelet coefficients
            //printf("Wavelets coefficients:\n");
            //for (i=0;i<32;i++) printf("wc{0}={1}",i.ToString(),x[i].ToString());
            //printf("\n");

            //// Do the inverse 9/7 transform
            //Dwt97.Reverse(ref x); 

            //// Prints the reconstructed signal 
            //printf("Reconstructed signal:\n");
            //for (i=0;i<32;i++) printf("xx{0}={1}",i.ToString(),x[i].ToString());
            //Console.Read();

            byte[,] y = ImageHelper.GetY("pic.jpg");
            Console.Read();

            Bitmap yBmp = ImageHelper.ConvertToY("pic.jpg");
            yBmp.Save("pic.y.png");
            Bitmap yEdge = ImageHelper.Transform(yBmp, ImageHelper.EDGE_LIGHT_KNL);
            yEdge.Save("pic.y.edge_light.png");

            Bitmap yEdgeMed = ImageHelper.Transform(yBmp, ImageHelper.EDGE_MED_KNL);
            yEdgeMed.Save("pic.y.edge_med.png");

            Bitmap yEdgeHev = ImageHelper.Transform(yBmp, ImageHelper.EDGE_HEAVY_KNL);
            yEdgeHev.Save("pic.y.edge_hev.png");

            Bitmap edgeHev = ImageHelper.Transform((Bitmap)Bitmap.FromFile("pic.jpg"), ImageHelper.EDGE_HEAVY_KNL);
            edgeHev.Save("pic.edge_hev.png");
        }

        static void printf(string value)
        {
            Console.WriteLine(value);
        }

        static void printf(string value, params string[] args)
        {
            Console.WriteLine(value, args);
        }
    }
}
