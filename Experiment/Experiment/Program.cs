using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
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

            //byte[,] y = ImageHelper.GetY("pic.jpg");
            //Console.Read();

            //Bitmap raw = new Bitmap("pic.jpg");
            //Bitmap scaled = ImageHelper.ScaleByPowerOf2(raw, 8);

            //Bitmap yBmp = ImageHelper.ConvertToY(scaled);
            //yBmp.Save("pic.y.png");
            //Bitmap yEdge = ImageHelper.Transform(yBmp, ImageHelper.EDGE_LIGHT_KNL);
            //yEdge.Save("pic.y.edge_light.png");

            //Bitmap yEdgeMed = ImageHelper.Transform(yBmp, ImageHelper.EDGE_MED_KNL);
            //yEdgeMed.Save("pic.y.edge_med.png");

            //Bitmap yEdgeHev = ImageHelper.Transform(yBmp, ImageHelper.EDGE_HEAVY_KNL);
            //yEdgeHev.Save("pic.y.edge_hev.png");

            //Bitmap edgeHev = ImageHelper.Transform(scaled, ImageHelper.EDGE_HEAVY_KNL);
            //edgeHev.Save("pic.edge_hev.png");

            //Bitmap reduce = null;
            //double[][] coeffs;
            //ImageHelper.Dwt97Y(yBmp, out reduce, out coeffs); // LL1
            //ImageHelper.Dwt97Y(reduce, out reduce, out coeffs); // LL2
            //ImageHelper.Dwt97Y(reduce, out reduce, out coeffs); // LL3
            //reduce.Save("pic.y.edge_hev.reduce.png");

            //Bitmap bmp_1 = new Bitmap("1.jpg");
            //Bitmap bmp_1_scaled = ImageHelper.ScaleByPowerOf2(bmp_1, 8);
            //Bitmap bmp_1_scaled_y = ImageHelper.ConvertToY(bmp_1_scaled);
            //Bitmap reduce = null;
            //double[][] coeffs;
            //ImageHelper.Dwt97Y(bmp_1_scaled_y, out reduce, out coeffs); // LL1
            //ImageHelper.Dwt97Y(reduce, out reduce, out coeffs); // LL2
            //ImageHelper.Dwt97Y(reduce, out reduce, out coeffs); // LL3
            //reduce.Save("1.LL3.png");
            //for (int r = 0; r < coeffs.Length; r++)
            //{
            //    for (int c = 0; c < coeffs[r].Length; c++)
            //    {
            //        Debug.Write(string.Format("{0}{1}", c>0 ? ",":"" , coeffs[r][c].ToString()));
            //    }
            //    Debug.WriteLine("");
            //}

            //Bitmap bmp_2 = new Bitmap("2.jpg");
            //Bitmap bmp_2_scaled = ImageHelper.ScaleByPowerOf2(bmp_2, 8);
            //Bitmap bmp_2_scaled_y = ImageHelper.ConvertToY(bmp_2_scaled);
            //Bitmap reduce2 = null;
            //double[][] coeffs2;
            //ImageHelper.Dwt97Y(bmp_2_scaled_y, out reduce2, out coeffs2); // LL1
            //ImageHelper.Dwt97Y(reduce2, out reduce2, out coeffs2); // LL2
            //ImageHelper.Dwt97Y(reduce2, out reduce2, out coeffs2); // LL3
            //reduce2.Save("2.LL3.png");
            //for (int r = 0; r < coeffs2.Length; r++)
            //{
            //    for (int c = 0; c < coeffs2[r].Length; c++)
            //    {
            //        Debug.Write(string.Format("{0}{1}", c > 0 ? "," : "", coeffs2[r][c].ToString()));
            //    }
            //    Debug.WriteLine("");
            //}
            //Point p = ImageHelper.FindMotionMaxLL3("1.jpg", "2.jpg", 8);
            //Point p2 = ImageHelper.FindMotionMaxLL3("2.jpg", "3.jpg", 8);

            
            string val = "1:1:0.5 [0.975] 2:1:2.2 [0.765] 3:1";
            Network<bool[]> network = Network<bool[]>.Create(val, 3, (idx, signal) => signal[idx] ? 1 : 0);
            double[] output = network.Update(new bool[] {true, false});
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
