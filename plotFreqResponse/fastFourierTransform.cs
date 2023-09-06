using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using FFTWSharp;
using System.Runtime.InteropServices;

namespace fftAndIfft
{
    class fastFourierTransform
    {
        public static double[] ifFT(double[] data)
        {
            // Get the length of the array
            int n = data.Length;
            /* Allocate an unmanaged memory block for the input and output data.
             * (The input and output are of the same length in this case, so we can use just one memory block.) */
            IntPtr ptr = fftw.malloc(n * 8);    // or: n * sizeof(double)
            // Pass the managed input data to the unmanaged memory block
            Marshal.Copy(data, 0, ptr, n);
            // Plan the IFFT and execute it (n/2 because complex numbers are stored as pairs of doubles)
            IntPtr plan = fftw.dft_1d(n / 2, ptr, ptr, fftw_direction.Backward, fftw_flags.Estimate);
            fftw.execute(plan);
            // Create an array to store the output values
            var ifft = new double[n];
            // Pass the unmanaged output data to the managed array
            Marshal.Copy(ptr, ifft, 0, n);
            // Do some cleaning
            fftw.destroy_plan(plan);
            fftw.free(ptr);
            fftw.cleanup();
            // Scale the output values
            for (int i = 0, nh = n / 2; i < n; i++)
                ifft[i] /= nh;
            // Return the IFFT output
            return ifft;
        }
        public static double[] fFT(double[] data, bool real)
        {
            // If the input is real, make it complex
            if (real)
                data = toComplex(data);
            // Get the length of the array
            int n = data.Length;
            /* Allocate an unmanaged memory block for the input and output data.
             * (The input and output are of the same length in this case, so we can use just one memory block.) */
            IntPtr ptr = fftw.malloc(n * 8);    // or: n * sizeof(double)
            // Pass the managed input data to the unmanaged memory block
            Marshal.Copy(data, 0, ptr, n);
            // Plan the FFT and execute it (n/2 because complex numbers are stored as pairs of doubles)
            IntPtr plan = fftw.dft_1d(n / 2, ptr, ptr, fftw_direction.Forward, fftw_flags.Estimate);
            fftw.execute(plan);
            // Create an array to store the output values
            var fft = new double[n];
            // Pass the unmanaged output data to the managed array
            Marshal.Copy(ptr, fft, 0, n);
            // Do some cleaning
            fftw.destroy_plan(plan);
            fftw.free(ptr);
            fftw.cleanup();
            // Return the FFT output
            return fft;
        }
        public static double[] toComplex(double[] real)
        {
            int n = real.Length;
            var comp = new double[n * 2];
            for (int i = 0; i < n; i++)
                comp[2 * i] = real[i];
            return comp;
        }
        public static void DisplayComplex(double[] x)
        {
            if (x.Length % 2 != 0)
                throw new Exception("The number of elements must be even.");
            for (int i = 0, n = x.Length; i < n; i += 2)
                if (x[i + 1] < 0)
                    Console.WriteLine(string.Format("{0} - {1}i", x[i], Math.Abs(x[i + 1])));
                else
                    Console.WriteLine(string.Format("{0} + {1}i", x[i], x[i + 1]));
        }

        /// <summary>
        /// Displays the real parts of complex numbers.
        /// </summary>
        /// <param name="x">An array of complex numbers.</param>
        public static void DisplayReal(double[] x)
        {
            if (x.Length % 2 != 0)
                throw new Exception("The number of elements must be even.");
            for (int i = 0, n = x.Length; i < n; i += 2)
                Console.WriteLine(x[i]);
        }
    }
}
