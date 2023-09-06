using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace differentialTerms
{
    class DifferentialEquations
    {
        private //Defining the main attributes of the objects within this class
            int yOrder;
            int xOrder;
            double[] yCoefs;
            double[] xCoefs;
            double[,] transferFunctionValues;
            double[,] qMatrixSS;
            double[,] bMatrixSS;
            double[,] cMatrixSS;
            double dUSS;

        public void init(double[] y, double[] x)
        {
            //Remember that the real order is yOrder - 1 and xOrder - 1 
            yOrder = (int)(y.Length);
            xOrder = (int)(x.Length);
            yCoefs = new double[yOrder];
            xCoefs = new double[xOrder];
        }
        public DifferentialEquations()
        {
            yOrder = 0;
            xOrder = 0;
        }
        public DifferentialEquations(double[] y, double[] x)
        {
            this.init(y, x);
            for (int i = 0; i < this.yOrder; i++)
                this.yCoefs[i] = y[i];
            for (int i = 0; i < this.xOrder; i++)
                this.xCoefs[i] = x[i];            
        }
        public void getTransferFunction()
        {
            int tFSize;
            if (this.yOrder >= this.xOrder) { tFSize = this.yOrder; }
            else { tFSize = this.xOrder; }
            this.transferFunctionValues = new double[2, tFSize];
            for(int i = (tFSize - 1); i >= (tFSize - xOrder); i--)
                this.transferFunctionValues[0, i] = this.xCoefs[i - (tFSize - this.xOrder)];
            for (int i = (tFSize - 1); i >= (tFSize - yOrder); i--)
                this.transferFunctionValues[1, i] = this.yCoefs[i - (tFSize - this.yOrder)];
            this.printTransferFunction();
        }
        public void getStateSpace()
        {
            int n = (int)this.yCoefs.Length - 1;
            this.qMatrixSS = new double[n, n];
            this.bMatrixSS = new double[n, 1];
            this.cMatrixSS = new double[1, n];
            this.dUSS = this.transferFunctionValues[0, 0];
            for (int row = 0; row < (n - 1); row++)
                for (int column = 1; column < n; column++)
                {
                    if (column == (row + 1)) { this.qMatrixSS[row, column] = 1; }
                    else { this.qMatrixSS[row, column] = 0; }
                }
            for (int column = 0; column < n; column++)
            {
                this.qMatrixSS[n - 1, column] = (-1) * this.transferFunctionValues[1, n - column];
                if (column != (n - 1)) { this.bMatrixSS[column, 0] = 0; }
                else { this.bMatrixSS[column, 0] = 1; }
                this.cMatrixSS[0, column] = this.transferFunctionValues[0, n - column] - this.transferFunctionValues[1, n - column] * this.transferFunctionValues[0, 0];
            }
        }
        public void printDifferentialEquation()
        {
            int n = (int)this.yCoefs.Length;
            Console.Write("\n\n" + "Differential Equation: ");
            for (int column = 0; column < n; column++)
            {
                if (column != (n - 1) && column != (n - 2)) { Console.Write(yCoefs[column] + " * d" + ((n - 1) - column) + "y/dt" + ((n - 1) - column) + " + "); }
                else if (column == (n - 2)) {Console.Write(yCoefs[column] + " * d" + "y/dt" + " + ");}
                else { Console.Write(yCoefs[column] + " * y"); }
            }
            Console.Write(" = ");
            n = (int)this.xCoefs.Length;
            for (int column = 0; column < n; column++)
            {
                if (column != (n - 1) && column != (n - 2)) { Console.Write(xCoefs[column] + " * d" + ((n - 1) - column) + "x/dt" + ((n - 1) - column) + " + "); }
                else if (column == (n - 2)) { Console.Write(xCoefs[column] + " * d" + "x/dt" + " + "); }
                else { Console.Write(xCoefs[column] + " * x"); }
            }
        
        }
        public void printTransferFunction()
        {
            int tFSize;
            if (this.yOrder >= xOrder) { tFSize = this.yOrder; }
            else { tFSize = this.xOrder; }
            Console.Write("\n\n" + "Transfer Function: ");
            for (int row = 0; row < 3; row++)
            {
                Console.Write("\n\t");
                for (int column = 0; column < tFSize; column++)
                {
                    if (row == 0)
                    {
                        if (column == 0) { Console.Write(this.transferFunctionValues[0, column] + " * s ^ " + (tFSize - column - 1)); }
                        else if (column == (tFSize - 2)) { Console.Write(" + " + this.transferFunctionValues[0, column] + " * s"); }
                        else if (column == (tFSize - 1)) { Console.Write(" + " + this.transferFunctionValues[0, column] + " "); }
                        else { Console.Write(" + " + this.transferFunctionValues[0, column] + " * s ^ " + (tFSize - column - 1)); }
                    }
                    else if (row == 1)
                    {
                        for (int i = 0; i < tFSize; i++)
                            Console.Write("--");
                    }
                    else 
                    {
                        if (column == 0) { Console.Write(this.transferFunctionValues[1, column] + " * s ^ " + (tFSize - column - 1)); }
                        else if (column == (tFSize - 2)) { Console.Write(" + " + this.transferFunctionValues[1, column] + " * s"); }
                        else if (column == (tFSize - 1)) { Console.Write(" + " + this.transferFunctionValues[1, column] + " "); }
                        else { Console.Write(" + " + this.transferFunctionValues[1, column] + " * s ^ " + (tFSize - column - 1)); }                    
                    }
                }
            }        
        }
        public void printSpaceStateMatrices()
        {
            int n = (int)Math.Sqrt(this.qMatrixSS.Length);
            Console.Write("\n\n" + "q Matrix: " + "\n");
            for (int row = 0; row < n; row++)
            {
                for (int column = 0; column < n; column++)
                { 
                    Console.Write(this.qMatrixSS[row, column] + " ");
                }
                Console.WriteLine("");
            }
            n = (int)this.bMatrixSS.Length;
            Console.Write("\n" + "b Matrix: " + "\n");
            for (int row = 0; row < n; row++)
                Console.WriteLine(this.bMatrixSS[row, 0]);
            n = (int)this.cMatrixSS.Length;
            Console.Write("\n" + "c Matrix: " + "\n");
            for (int column = 0; column < n; column++)
                Console.Write(this.cMatrixSS[0, column] + " ");      
        }
        public static void series(DifferentialEquations h1, DifferentialEquations h2)
        {
            DifferentialEquations hSeries = new DifferentialEquations();
            int n1 = (int)(h1.transferFunctionValues.Length / 2);
            int n2 = (int)(h2.transferFunctionValues.Length / 2);
            int numH1 = 0;
            int denH1 = 0;
            int numH2 = 0;
            int denH2 = 0;
            for (int column = 0; column < n1; column++)
            {
                if (h1.transferFunctionValues[0, column] != 0) { numH1++; }
                if (h1.transferFunctionValues[1, column] != 0) { denH1++; }
            }
            
            for (int column = 0; column < n2; column++)
            {
                if (h2.transferFunctionValues[0, column] != 0) { numH2++; }
                if (h2.transferFunctionValues[1, column] != 0) { denH2++; }
            }

            hSeries.xOrder = n1 + n2 - 1;
            hSeries.yOrder = n1 + n2 - 1;

            //hSeries.xOrder = numH1 + numH2 - 1;
            //hSeries.yOrder = denH1 + denH2 - 1;
            int m = 0;
            if (hSeries.xOrder > hSeries.yOrder) { m = hSeries.xOrder; }
            else { m = hSeries.yOrder; }
            hSeries.xOrder = m;
            hSeries.yOrder = m;
            hSeries.transferFunctionValues = new double[2, m];

            for (int i = 0; i < n1; i++)
                for (int j = 0; j < n2; j++)
                {
                    hSeries.transferFunctionValues[0, i + j] += h1.transferFunctionValues[0, i] * h2.transferFunctionValues[0, j];
                    hSeries.transferFunctionValues[1, i + j] += h1.transferFunctionValues[1, i] * h2.transferFunctionValues[1, j];
                }

                    
            hSeries.printTransferFunction();

        }
        public static void feedback(DifferentialEquations g, DifferentialEquations h, int sign)
        {
            DifferentialEquations hFeedback = new DifferentialEquations();
            int n1 = (int)(g.transferFunctionValues.Length / 2);
            int n2 = (int)(h.transferFunctionValues.Length / 2);
            int numG = 0;
            int denG = 0;
            int numH = 0;
            int denH = 0;
            for (int column = 0; column < n1; column++)
            {
                if (g.transferFunctionValues[0, column] != 0) { numG++; }
                if (g.transferFunctionValues[1, column] != 0) { denG++; }
            }

            for (int column = 0; column < n2; column++)
            {
                if (h.transferFunctionValues[0, column] != 0) { numH++; }
                if (h.transferFunctionValues[1, column] != 0) { denH++; }
            }

            hFeedback.xOrder = numG + denH - 1;
            
            int denFd1term;
            int denFd2term;

            if ((denG + denH - 1) >= (numG + numH - 1)) 
            { 
                hFeedback.yOrder = denG + denH - 1;
                denFd1term = denG;
                denFd2term = denH;
            }
            else 
            { 
                hFeedback.yOrder = numG + numH - 1;
                denFd1term = numG;
                denFd2term = numH;
            }
            int m = 0;
            if (hFeedback.xOrder > hFeedback.yOrder) { m = hFeedback.xOrder; }
            else { m = hFeedback.yOrder; }
            hFeedback.xOrder = m;
            hFeedback.yOrder = m;
            hFeedback.transferFunctionValues = new double[2, m];

            for (int i = 0; i < n1; i++)
                for (int j = 0; j < n2; j++)
                    hFeedback.transferFunctionValues[0, i + j] += g.transferFunctionValues[0, i] * h.transferFunctionValues[1, j];

            if (sign == -1)
            {
                for (int i = 0; i < denFd1term; i++)
                    for (int j = 0; j < denFd2term; j++)
                        hFeedback.transferFunctionValues[1, i + j] += g.transferFunctionValues[1, i] * h.transferFunctionValues[1, j] - g.transferFunctionValues[0, i] * h.transferFunctionValues[0, j];
            }
            else
            {
                for (int i = 0; i < n1; i++)
                    for (int j = 0; j < n2; j++)
                        hFeedback.transferFunctionValues[1, i + j] += g.transferFunctionValues[1, i] * h.transferFunctionValues[1, j] + g.transferFunctionValues[0, i] * h.transferFunctionValues[0, j];          
            }
            
            hFeedback.printTransferFunction();

        }
        
        public void timeRespSim(double sT, double[] inputSignal, double x0, double xp0, out double[] yp) 
        //sT: Sampling Time and inputSignal is a vector of 1xN containing information
        //of the time domain input signal, x0 is the initial position and xp0 is the initial speed
        { 
            int qn = (int)Math.Sqrt(this.qMatrixSS.Length);
            double[] qvector = new double[(int)this.qMatrixSS.Length];
            for (int i = 0; i < qn; i++)
                for (int j = 0; j < qn; j++)
                    qvector[i*qn + j] = this.qMatrixSS[i, j];
            Matrix q = new Matrix(qn, qn, qvector);
            int bn = (int)this.bMatrixSS.Length;
            double[] bvector = new double[(int)this.bMatrixSS.Length];
            for (int i = 0; i < bn; i++)
                for (int j = 0; j < 1; j++)
                    bvector[i * 1 + j] = this.bMatrixSS[i, j];
            Matrix b = new Matrix(bn, 1, bvector);
            int cn = (int)this.cMatrixSS.Length;
            double[] cvector = new double[(int)this.bMatrixSS.Length];
            for (int i = 0; i < 1; i++)
                for (int j = 0; j < cn; j++)
                    cvector[i * cn + j] = this.cMatrixSS[i, j];
            Matrix c = new Matrix(1, cn, cvector);
            int n = (int)Math.Sqrt(this.qMatrixSS.Length);
            int m = (int)inputSignal.Length;
            double[] y = new double[n*m];
            yp = new double[1 * m];
            for (int i = 0; i < n; i++)
                for (int j = 1; j < m; j++)
                    y[i * m + j] = 0;
            //for (int j = 0; j < n; j++)
            //    y[j*m] = 0;
            y[0 * m + 0] = x0;
            //y[1 * m + 0] = xp0;            
            Matrix ftdt = new Matrix(n, 1);            
            Matrix DX = new Matrix((int)this.bMatrixSS.Length, (int)inputSignal.Length, y);  
            Matrix Y = new Matrix();
            for (int i = 1; i < m; i++)
            {
                ftdt = (q * DX.getMatrixPart(0, n - 1, i - 1, i - 1) + b * inputSignal[i - 1]) * sT;
                ftdt = DX.getMatrixPart(0, n - 1, i - 1, i - 1) + ftdt;
                DX.setMatrixPart(0, n - 1, i, i, ftdt);
                DX.getMatrixVector();
                //yp[i] = c*DX.outVector[0*m + i];
                Y = c*DX.getMatrixPart(0, n - 1, i, i);                
                Y.getMatrixVector();
                yp[i] = Y.outVector[0];
            }


            //int n = (int)Math.Sqrt(this.qMatrixSS.Length);
            //int m = (int)inputSignal.Length;
            //y = new double[m];            
            //y[0] = x0;
            //double yp0 = xp0;
            //for (int i = 1; i < m; i++)
            //    y[i] = y[i - 1] + sT * (this.qMatrixSS[0, 0] * y[i - 1] + this.bMatrixSS[0, 0] * inputSignal[i - 1]);
        }

        public void rungeKutta4thO(double sT, double[] inputSignal, double x0, double xp0, out double[] yp)
        {
            int qn = (int)Math.Sqrt(this.qMatrixSS.Length);
            double[] qvector = new double[(int)this.qMatrixSS.Length];
            for (int i = 0; i < qn; i++)
                for (int j = 0; j < qn; j++)
                    qvector[i * qn + j] = this.qMatrixSS[i, j];
            Matrix q = new Matrix(qn, qn, qvector);
            int bn = (int)this.bMatrixSS.Length;
            double[] bvector = new double[(int)this.bMatrixSS.Length];
            for (int i = 0; i < bn; i++)
                for (int j = 0; j < 1; j++)
                    bvector[i * 1 + j] = this.bMatrixSS[i, j];
            Matrix b = new Matrix(bn, 1, bvector);
            int cn = (int)this.cMatrixSS.Length;
            double[] cvector = new double[(int)this.bMatrixSS.Length];
            for (int i = 0; i < 1; i++)
                for (int j = 0; j < cn; j++)
                    cvector[i * cn + j] = this.cMatrixSS[i, j];
            Matrix c = new Matrix(1, cn, cvector);
            int n = (int)Math.Sqrt(this.qMatrixSS.Length);
            int m = (int)inputSignal.Length;
            double[] y = new double[n * m];
            yp = new double[1 * m];
            for (int i = 0; i < n; i++)
                for (int j = 1; j < m; j++)
                    y[i * m + j] = 0;
            y[0 * m + 0] = x0;
            //y[1 * m + 0] = xp0;            
            Matrix ftdt = new Matrix(n, 1);
            Matrix DX = new Matrix((int)this.bMatrixSS.Length, (int)inputSignal.Length, y);
            Matrix Y = new Matrix();
            Matrix anR = new Matrix();
            Matrix bnR = new Matrix();
            Matrix cnR = new Matrix();
            Matrix dnR = new Matrix();
            Matrix xNew = new Matrix(n, 1); 
            for (int i = 1; i < m; i++)
            {
                anR = (q * DX.getMatrixPart(0, n - 1, i - 1, i - 1) + b * inputSignal[i - 1]);
                bnR = (q * (DX.getMatrixPart(0, n - 1, i - 1, i - 1) + anR * (sT / 2)) + b * inputSignal[i - 1]);
                cnR = (q * (DX.getMatrixPart(0, n - 1, i - 1, i - 1) + bnR * (sT / 2)) + b * inputSignal[i - 1]);
                dnR = (q * (DX.getMatrixPart(0, n - 1, i - 1, i - 1) + cnR * sT) + b * inputSignal[i - 1]);
                xNew = (anR + bnR * 2 + cnR * 2 + dnR) * (sT / 6);
                xNew = DX.getMatrixPart(0, n - 1, i - 1, i - 1) + xNew;
                DX.setMatrixPart(0, n - 1, i, i, xNew);
                DX.getMatrixVector();
                Y = c * DX.getMatrixPart(0, n - 1, i, i);
                Y.getMatrixVector();
                yp[i] = Y.outVector[0];
                
            }
        }

    }
}
