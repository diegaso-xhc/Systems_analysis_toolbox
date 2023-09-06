using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace System_response_Toolbox
{
    class backupFunctions
    {
        public static void unWrap(double[] input, out double[] unWrVec, double tshold)
        {
            int n = (int)input.Length;
            unWrVec = new double[n];
            unWrVec[0] = input[0];
            double f;
            for (int i = 1; i <= n - 1; i++)
            {
                if (Math.Abs(input[i] - input[i - 1]) > tshold)
                {
                    f = Math.Abs(input[i] - input[i - 1]) / Math.PI;
                    if (Math.Ceiling(f) % 2 != 0) f = Math.Floor(f);
                    else if (Math.Ceiling(f) % 2 == 0) f = Math.Ceiling(f);
                    if (f == 0) f = 2;
                    if ((input[i] - input[i - 1]) > 0)
                        input[i] = input[i] - f * Math.PI;
                    else
                        input[i] = input[i] + f * Math.PI;
                }
                unWrVec[i - 1] = input[i - 1];
            }
            unWrVec[n - 1] = input[n - 1];
        }
        public static double angle(double x, double y)
        {
            double angle;
            if (x < 0)
            {
                if (y > 0) angle = Math.PI - Math.Atan(y / Math.Abs(x));
                else if (y < 0) angle = Math.Atan(Math.Abs(y) / Math.Abs(x)) - Math.PI;
                else angle = Math.PI;
            }
            else if (x > 0)
            {
                if (y > 0) angle = Math.Atan(y / x);
                else if (y < 0) angle = -Math.Atan(Math.Abs(y) / x);
                else angle = 0;
            }
            else
            {
                if (y > 0) angle = Math.PI / 2;
                else if (y < 0) angle = 3 * Math.PI / 2;
                else angle = 0;
            }
            return angle;
        }
        public static void printVector(double[] inputVec)
        {
            Console.WriteLine("Vector:");
            for (int i = 0; i < (int)inputVec.Length; i++)
                Console.WriteLine(inputVec[i]);
        }
        public static void getRoots(double[] polynom, out double[] reig, out double[] imeig)
        {
            int n = polynom.Length - 1;
            double[,] A = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i == 0)
                    {
                        A[i, j] = -1.0 * polynom[j + 1];
                    }
                    else
                    {
                        if (i - 1 == j) A[i, j] = 1;
                        else A[i, j] = 0;
                    }

                }
            }
            reig = new double[n];
            imeig = new double[n];
            double[,] vl = new double[n, n];
            double[,] vr = new double[n, n];
            alglib.rmatrixevd(A, n, 0, out reig, out imeig, out vl, out vr);
        }
        public static double angle0to360(double x, double y)
        {
            double angle;
            if (x < 0)
            {
                if (y > 0) angle = Math.PI - Math.Atan(y / Math.Abs(x));
                else if (y < 0) angle = Math.Atan(Math.Abs(y) / Math.Abs(x)) + Math.PI;
                else angle = Math.PI;
            }
            else if (x > 0)
            {
                if (y > 0) angle = Math.Atan(y / x);
                else if (y < 0) angle = 2 * Math.PI - Math.Atan(Math.Abs(y) / x);
                else angle = 0;
            }
            else
            {
                if (y > 0) angle = Math.PI / 2;
                else if (y < 0) angle = 3 * Math.PI / 2;
                else angle = 0;
            }
            return angle;
        }
    }
}
