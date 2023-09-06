using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace differentialTerms
{
    class Matrix
    {
        private //Defining the main attributes of the objects within this class
            double[] vector;
            int rw;
            int cl;
        public
            double[] outVector;

        public void init(int m, int n) //Initializing the attributes of this class by passing two arguments: The number of rows and columns
        {
            vector = new double[m * n]; //This vector will be a one-row vector and will contain the information of a matrix
            outVector = new double[m * n];
            rw = m; //Number of rows
            cl = n; //Number of columns
        }
        public Matrix() //Initializing a Matrix object without the need of passing an argument to the constructor
        {
            rw = 0;
            cl = 0;
        }
        public Matrix(int m, int n)
        {
            init(m, n);
            for (int row = 0; row < m; row++)
                for (int column = 0; column < n; column++)
                    vector[row * n + column] = 0;
        }
        public Matrix(int m, int n, double[] vector) //Initializing a Matrix object by passing not only the number of rows and columns but also its
        //corresponding elements
        {
            init(m, n);
            for (int row = 0; row < m; row++)
                for (int column = 0; column < n; column++)
                    this.vector[row * n + column] = vector[row * n + column];
        }
        public void print() //A function for printing the vectors of each object in a matrix form so that it can be easily understood by the user
        {
            Console.Write("\n");
            for (int row = 0; row < rw; row++)
            {
                for (int column = 0; column < cl; column++)
                    Console.Write("\t" + vector[row * cl + column] + "\t");
                Console.Write("\n");
            }
        }
        public static Matrix operator +(Matrix a, Matrix b) //Overcharge of the "+" Operator
        {
            if (a.rw != b.rw || a.cl != b.cl) return new Matrix();
            else
            {
                Matrix Ans = new Matrix(a.rw, a.cl);
                for (int row = 0; row < a.rw; row++)
                {
                    for (int column = 0; column < a.cl; column++)
                    {
                        Ans.vector[row * a.cl + column] = a.vector[row * a.cl + column] + b.vector[row * a.cl + column];
                    }
                }
                return Ans;
            }
        }
        public static Matrix operator -(Matrix a, Matrix b) //Overcharge of the "-" operator
        {
            if (a.rw != b.rw || a.cl != b.cl) return new Matrix();
            else
            {
                Matrix Ans = new Matrix(a.rw, a.cl);
                for (int row = 0; row < a.rw; row++)
                {
                    for (int column = 0; column < a.cl; column++)
                    {
                        Ans.vector[row * a.cl + column] = a.vector[row * a.cl + column] - b.vector[row * a.cl + column];
                    }
                }
                return Ans;
            }
        }

        public static Matrix operator %(Matrix a, Matrix b)
        {
            if (a.rw != b.rw || a.cl != b.cl) return new Matrix();
            else
            {
                Matrix Ans = new Matrix(a.rw, a.cl);
                for (int row = 0; row < a.rw; row++)
                {
                    for (int column = 0; column < a.cl; column++)
                    {
                        Ans.vector[row * a.cl + column] = a.vector[row * a.cl + column] * b.vector[row * a.cl + column];
                    }
                }
                return Ans;
            }
        }

        public static Matrix operator *(Matrix a, Matrix b) //Overcharge of the "*" operator by multiplying two matrices
        {
            if (a.cl != b.rw) return new Matrix();
            else
            {
                Matrix Ans = new Matrix(a.rw, b.cl);
                for (int row = 0; row < a.rw; row++)
                {
                    for (int column = 0; column < b.cl; column++)
                    {
                        for (int k = 0; k < a.cl; k++)
                        {
                            Ans.vector[row * b.cl + column] += a.vector[row * a.cl + k] * b.vector[k * b.cl + column];
                        }
                    }
                }
                return Ans;
            }
        }
        public static Matrix operator *(Matrix a, double b) //Overcharge of the "*" by multiplying a matrix by a scalar value
        {
            Matrix Ans = new Matrix(a.rw, a.cl);
            for (int row = 0; row < a.rw; row++)
            {
                for (int column = 0; column < a.cl; column++)
                {
                    Ans.vector[row * a.cl + column] = a.vector[row * a.cl + column] * b;
                }
            }
            return Ans;
        }
        public static Matrix operator /(Matrix a, Matrix b) //Overcharge of the "/" operator
        {
            Matrix Ans = new Matrix(a.rw, a.cl);
            Ans = a * b.inverseMatrix();
            return Ans;
        }
        public Matrix transpose() //Function for getting the transpose matrix of a given matrix
        {
            Matrix Ans = new Matrix(this.cl, this.rw);
            for (int row = 0; row < this.cl; row++)
            {
                for (int column = 0; column < this.rw; column++)
                {
                    Ans.vector[row * this.rw + column] = this.vector[column * this.cl + row];
                }
            }
            return Ans;
        }
        public double determinant() //Function for getting the determinant value of a given matrix
        {
            double detVal = 0;
            int col2;
            if (this.cl == this.rw)
            {
                if (this.cl == 1)
                    detVal = this.vector[0];
                else if (this.cl == 2)
                    detVal = determinantOperation(this.vector);
                else
                {
                    for (int column = 0; column < this.cl; column++)
                    {
                        Matrix newVector;
                        double[] v = new double[(this.rw - 1) * (this.cl - 1)];
                        for (int row = 0; row < this.rw - 1; row++)
                        {
                            for (int col = 0; col < this.cl - 1; col++)
                            {
                                if (col >= column) { col2 = col + 1; }
                                else { col2 = col; }
                                v[row * (this.cl - 1) + col] = this.vector[(row + 1) * (this.cl) + col2];
                            }
                        }
                        newVector = new Matrix(this.rw - 1, this.cl - 1, v);
                        if (column % 2 == 0)
                            detVal += newVector.determinant() * this.vector[column];
                        else
                            detVal -= newVector.determinant() * this.vector[column];
                    }
                }
            }
            return detVal;
        }

        public double determinantOperation(double[] twoByTwoVector) //Function for obtaining the determinant value of a two by two matrix
        {
            double detVal;
            detVal = twoByTwoVector[0] * twoByTwoVector[3] - twoByTwoVector[1] * twoByTwoVector[2];
            return detVal;
        }
        public Matrix cofactorMatrix() //Method for getting the cofactor Matrix of an object within this class
        {
            Matrix Ans = new Matrix(this.rw, this.cl);
            int row2, col2;
            if (this.cl == this.rw)
            {
                for (int row = 0; row < this.rw; row++)
                {
                    for (int column = 0; column < this.cl; column++)
                    {
                        Matrix newVector;
                        double[] v = new double[this.rw * this.cl];
                        for (int rowd = 0; rowd < this.rw - 1; rowd++)
                        {
                            if (rowd >= row) { row2 = rowd + 1; }
                            else { row2 = rowd; }

                            for (int columnd = 0; columnd < this.cl - 1; columnd++)
                            {
                                if (columnd >= column) { col2 = columnd + 1; }
                                else { col2 = columnd; }

                                v[rowd * (this.cl - 1) + columnd] = this.vector[row2 * (this.cl) + col2];
                            }
                        }
                        newVector = new Matrix(this.rw - 1, this.cl - 1, v);
                        Ans.vector[row * this.cl + column] = (long)Math.Pow(-1, (row + column)) * newVector.determinant();
                    }
                }
            }
            return Ans;
        }
        public Matrix inverseMatrix() //Method for obtaining the inverse matrix of an object
        {
            Matrix Ans = new Matrix(this.rw, this.cl);
            if (this.cl == this.rw)
            {
                Ans = this.cofactorMatrix().transpose() * (1 / this.determinant());
            }
            return Ans;
        }
        public void getMatrixVector()
        { 
             for (int i = 0; i < this.rw; i++)
                for (int j = 0; j < this.cl; j++)
                    this.outVector[i * this.cl + j] = this.vector[i * this.cl + j];        
        }
        public Matrix setMatrixVector(int row, int col, double nVal)
        {
            Matrix newM = new Matrix(this.rw, this.cl, this.vector);
            newM.vector[row * this.cl + col] = nVal;
            return newM;
        }
        public Matrix getMatrixPart(int row1, int row2, int col1, int col2)
        {
            int m = row2 - row1 + 1;
            int n = col2 - col1 + 1;
            double[] tempOut = new double[m*n];
            for (int i = row1; i <= row2; i++)
                for (int j = col1; j <= col2; j++)
                    tempOut[(i - row1) * n + (j - col1)] = this.vector[i * this.cl + j];
            Matrix newM = new Matrix(m, n, tempOut);            
            return newM;
        }
        public void setMatrixPart(int row1, int row2, int col1, int col2, Matrix a)
        {
            int m = row2 - row1 + 1;
            int n = col2 - col1 + 1;
            for (int i = row1; i <= row2; i++)
                for (int j = col1; j <= col2; j++)
                    this.vector[i * this.cl + j] = a.vector[(i - row1) * n + (j - col1)];                    
        }

        
    }
}
