using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Experiment
{
    public static class ImageHelper
    {
        public static readonly float[,] IDENTITY_KNL = new float[3, 3] { { 0, 0, 0 }, { 0, 1, 0 }, { 0, 0, 0 } };
        public static readonly float[,] EDGE_LIGHT_KNL = new float[3, 3] { { 1, 0, -1 }, { 0, 0, 0 }, { -1, 0, 1 } };
        public static readonly float[,] EDGE_MED_KNL = new float[3, 3] { { 0, 1, 0 }, { 1, -4, 1 }, { 0, 1, 0 } };
        public static readonly float[,] EDGE_HEAVY_KNL = new float[3, 3] { { -1, -1, -1 }, { -1, 8, -1 }, { -1, -1, -1 } };
        public static readonly float[,] SHARPEN_KNL = new float[3, 3] { { 0, -1, 0 }, { -1, 5, -1 }, { 0, -1, 0 } };
        public static readonly float[,] BLUR_KNL = new float[3, 3] { { 1f / 9f, 1f / 9f, 1f / 9f }, { 1f / 9f, 1f / 9f, 1f / 9f }, { 1f / 9f, 1f / 9f, 1f / 9f } };
        public static readonly float[,] GAUSSIAN_BLUR_KNL = new float[3, 3] { { 1f / 16f, 2f / 16f, 1f / 16f }, { 2f / 16f, 4f / 16f, 2f / 16f }, { 1f / 16f, 2f / 16f, 1f / 16f } };
        public static readonly float[,] UNSHARP_KNL = new float[5, 5] 
        { 
            { -1f / 256f, -4f / 256f,-6/ 256f , -4f/256f, -1f/256f}, 
            { -4f / 256f, -16f / 256f,-24/ 256f , -16f/256f, -4f/256f}, 
            { -6f / 256f, -24f / 256f,476/ 256f , -24f/256f, -6f/256f},
            { -4f / 256f, -16f / 256f,-24/ 256f , -16f/256f, -4f/256f}, 
            { -1f / 256f, -4f / 256f,-6/ 256f , -4f/256f, -1f/256f}
        };


        public static byte[,] GetY(string bitmapPath)
        {
            return GetY((Bitmap)Bitmap.FromFile(bitmapPath));
        }

        public unsafe static byte[,] GetY(Bitmap rgbImage)
        {
            BitmapData bd = rgbImage.LockBits(new Rectangle(new Point(0, 0), new Size(rgbImage.Width, rgbImage.Height)),
                 System.Drawing.Imaging.ImageLockMode.ReadOnly,
                 rgbImage.PixelFormat);
            byte[,] y = new byte[rgbImage.Height, rgbImage.Width];

            int bpp = 24;

            try
            {
                byte* pixel = (byte*)bd.Scan0;
                byte* pixel0 = pixel;

                for (int h = 0; h < bd.Height; h++)
                {
                    for (int w = 0; w < bd.Width; w++)
                    {

                        pixel = pixel0 + h * bd.Stride + w * bpp / 8;
                        y[h,w] = Clamp(((66 * (byte)pixel[0] + 129 * (byte)(pixel[1]) +  25 * (byte)(pixel[2]) + 128) >> 8) +  16);
                    }
                }
            }
            catch { }
            finally
            {
                rgbImage.UnlockBits(bd);
            }

            return y;
        }

        public unsafe static Bitmap FromY(byte[,] y)
        {
            Bitmap bmp = new Bitmap(y.GetLength(1), y.GetLength(0), PixelFormat.Format24bppRgb);
            BitmapData bd = bmp.LockBits(new Rectangle(new Point(0,0), new Size(bmp.Width, bmp.Height)), 
                System.Drawing.Imaging.ImageLockMode.ReadWrite,
                bmp.PixelFormat);

            int bpp = 24;

            try
            {
                byte* pixel = (byte*)bd.Scan0;
                byte* pixel0 = pixel;
                byte C;

                for (int h = 0; h < bd.Height; h++)
                {
                    for (int w = 0; w < bd.Width; w++)
                    {
                        pixel = pixel0 + h * bd.Stride + w * bpp / 8;
                        C = (byte)(y[h, w] - 16);

                        pixel[0] = C;
                        pixel[1] = C;
                        pixel[2] = C;
                    }
                }
            }
            catch { }
            finally
            {
                bmp.UnlockBits(bd);
            }

            return bmp;
        }

        public static Bitmap ConvertToY(string bitmapPath)
        {
            return ConvertToY((Bitmap)Bitmap.FromFile(bitmapPath));
        }

        public static Bitmap ConvertToY(Bitmap bitmap)
        {
            return FromY(GetY(bitmap));
        }

        public unsafe static Bitmap Transform(Bitmap bitmap, float[,] kernel)
        {
            int kh = kernel.GetLength(0); // kernel height;
            int kw = kernel.GetLength(1); // kernel width;

            if (kh % 2 == 0 || kw % 2 == 0) throw new ArgumentException("Kernel width and height must be odd numbers.");

            Bitmap bmpDest = new Bitmap(bitmap.Width, bitmap.Height, PixelFormat.Format24bppRgb);
            BitmapData bdDest = bmpDest.LockBits(new Rectangle(new Point(0, 0), new Size(bmpDest.Width, bmpDest.Height)),
                System.Drawing.Imaging.ImageLockMode.ReadWrite,
                bmpDest.PixelFormat);
            BitmapData bdSource = bitmap.LockBits(new Rectangle(new Point(0, 0), new Size(bitmap.Width, bitmap.Height)),
                System.Drawing.Imaging.ImageLockMode.ReadOnly,
                bitmap.PixelFormat);

            int bpp = 24;

            try
            {
                byte* pixelDest = (byte*)bdDest.Scan0;
                byte* pixelDest0 = pixelDest;
                byte* pixelSource = (byte*)bdSource.Scan0;
                byte* pixelSource0 = pixelSource;
                

                for (int h = kh / 2 + 1; h < bdDest.Height - (kh /2 + 1); h++)
                {
                    for (int w = kw / 2 + 1; w < bdDest.Width - (kw /2 + 1); w++)
                    {
                        pixelDest = pixelDest0 + h * bdDest.Stride + w * bpp / 8;
                        pixelSource = pixelSource0 + h * bdSource.Stride + w * bpp / 8;
                        Convolve(pixelSource, pixelDest, bdDest.Stride, kernel, bpp);
                    }
                }
            }
            catch { }
            finally
            {
                bmpDest.UnlockBits(bdDest);
                bitmap.UnlockBits(bdSource);
            }

            return bmpDest;
        }

        private unsafe static void Convolve(byte* pixelSource, byte* pixelDest, int stride, float[,] kernel, int bpp)
        {
            int kh = kernel.GetLength(0);
            int kw = kernel.GetLength(1);
            int kh2 = kh / 2;
            int kw2 = kw / 2;
            float sumR = 0f, sumG = 0f, sumB = 0f;
            byte* psk = pixelSource;
            byte* pdk = pixelDest;

            for (int h = -(kh2 + 1); h < (kh2); h++)
            {
                for (int w = -(kw2 + 1); w < (kw2); w++)
                {
                    float k = kernel[h + kh2 + 1,w + kw2 + 1];
                    psk = pixelSource + h * stride + w * bpp / 8;
                    sumR += k * psk[0];
                    sumG += k * psk[1];
                    sumB += k * psk[2];
                }
            }

            pixelDest[0] = Clamp(sumR);
            pixelDest[1] = Clamp(sumG);
            pixelDest[1] = Clamp(sumB);
        }

        public unsafe static void Dwt97Y(Bitmap source, out Bitmap reduction, out double[][] coefficients)
        {
            double logWidth = Math.Log(source.Width, 2);
            if (logWidth % 1 > double.Epsilon) 
                throw new ArgumentException("Image width must be a factor of 2^N.");
            int maxHeight = (source.Height/2)*2; // cuts off off numbers
            reduction = new Bitmap(source.Width/2, source.Height/2, PixelFormat.Format24bppRgb);
            coefficients = new double[source.Height/2][];
            BitmapData bdSource = source.LockBits(new Rectangle(new Point(0, 0), new Size(source.Width, source.Height)),
                System.Drawing.Imaging.ImageLockMode.ReadOnly,
                source.PixelFormat);
            BitmapData bdReduction = reduction.LockBits(new Rectangle(new Point(0, 0), new Size(reduction.Width, reduction.Height)),
                System.Drawing.Imaging.ImageLockMode.ReadWrite,
                reduction.PixelFormat);
            int bpp = 24;

            try
            {
                byte* pixel = (byte*)bdSource.Scan0;
                byte* pixel0 = pixel;

                byte* pixelR = (byte*)bdReduction.Scan0;
                byte* pixelR0 = pixelR;

                for (int h = 0; h < maxHeight; h += 2)
                {
                    double[] row = new double[source.Width];
                    coefficients[h/2] = new double[source.Width/2];
                    for (int w = 0; w < source.Width; w++)
                    {
                        pixel = pixel0 + h * bdSource.Stride + w * bpp / 8;
                        row[w] = (double)(int)Clamp(((66 * (byte)pixel[0] + 129 * (byte)(pixel[1]) + 25 * (byte)(pixel[2]) + 128) >> 8) + 16);
                    }
                    Dwt97.Transform(ref row);
                    int ww = 0;
                    byte C = 0;
                    for (ww = 0; ww < row.Length / 2; ww++)
                    {
                        pixelR = pixelR0 + (h/2) * bdReduction.Stride + ww * bpp / 8;
                        C = (byte)((byte)row[ww] - 16);

                        pixelR[0] = C;
                        pixelR[1] = C;
                        pixelR[2] = C;
                    }
                    int ww0 = ww;
                    for (ww = ww0; ww < row.Length; ww++)
                    {
                        coefficients[h/2][ww - ww0] = row[ww];
                    }
                }
            }
            finally
            {
                source.UnlockBits(bdSource);
                reduction.UnlockBits(bdReduction);
            }

        }

        public static Bitmap ScaleByPowerOf2(Bitmap source, byte factor)
        {
            double newWidth = Math.Pow(2, factor);
            double scaling = newWidth/(double)source.Width;
            double newHeight = source.Height*scaling;
            Size newSize = new Size((int)newWidth, (int)newHeight);
            Bitmap dest = new Bitmap(newSize.Width, newSize.Height, PixelFormat.Format24bppRgb);
            Graphics g = Graphics.FromImage(dest);
            g.DrawImage(source,new Rectangle(new Point(0,0), newSize));
            g.Dispose();
            return dest;
        }

        public static byte Clamp(int value)
        {
            if (value < 0) return 0;
            if (value > 255) return 255;
            return (byte)value;
        }

        public static byte Clamp(float value)
        {
            if (value < 0) return 0;
            if (value > 255) return 255;
            return (byte)value;
        }

        public static Point FindMotionMaxLL3(string bitmapT0Path, string bitmapT1Path, byte scaledBy2Power = 0)
        {
            Bitmap t0 = new Bitmap(bitmapT0Path);
            Bitmap t1 = new Bitmap(bitmapT1Path);
            if (scaledBy2Power > 0)
            {
                t0 = ScaleByPowerOf2(t0, scaledBy2Power);
                t1 = ScaleByPowerOf2(t1, scaledBy2Power);
            }
            return FindMotionMaxLL3(t0, t1);
        }

        public static Point FindMotionMaxLL3(Bitmap imaget0, Bitmap imaget1)
        {
            Bitmap reduce0 = null;
            double[][] coeffs0;
            ImageHelper.Dwt97Y(imaget0, out reduce0, out coeffs0); // LL1
            ImageHelper.Dwt97Y(reduce0, out reduce0, out coeffs0); // LL2
            ImageHelper.Dwt97Y(reduce0, out reduce0, out coeffs0); // LL3

            Bitmap reduce1 = null;
            double[][] coeffs1;
            ImageHelper.Dwt97Y(imaget1, out reduce1, out coeffs1); // LL1
            ImageHelper.Dwt97Y(reduce1, out reduce1, out coeffs1); // LL2
            ImageHelper.Dwt97Y(reduce1, out reduce1, out coeffs1); // LL3

            double[][] absDiffs = new double[coeffs0.Length][];
            for (int r = 0; r < coeffs0.Length; r++)
            {
                absDiffs[r] = new double[coeffs0[r].Length];
                for (int c = 0; c < absDiffs[r].Length; c++)
                {
                    absDiffs[r][c] = Math.Abs(coeffs1[r][c] - coeffs0[r][c]);
                }
            }

            Point maxPoint = new Point(imaget0.Width/2, imaget0.Height/2);
            double[] rowSums = new double[coeffs0.Length];
            for (int r = 0; r < absDiffs.Length; r++)
            {
                for (int c = 0; c < absDiffs[r].Length; c++)
                {
                    rowSums[r] += absDiffs[r][c];
                }
            }

            double[] colSums = new double[coeffs0[0].Length];
            for (int r = 0; r < absDiffs.Length; r++)
            {
                for (int c = 0; c < absDiffs[r].Length; c++)
                {
                    colSums[c] += absDiffs[r][c];
                }
            }

            double max = double.MinValue;
            for (int r = 0; r < rowSums.Length; r++)
            {
                for (int c = 0; c < colSums.Length; c++)
                {
                    double testMax = rowSums[r]*colSums[c];
                    if (testMax > max)
                    {
                        maxPoint.X = c;
                        maxPoint.Y = r;
                        max = testMax;
                    }
                }
            }

            maxPoint.X *= 8;
            maxPoint.Y *= 8;

            return maxPoint;
        }

    }
}
