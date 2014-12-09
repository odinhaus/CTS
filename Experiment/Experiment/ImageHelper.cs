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

    }
}
