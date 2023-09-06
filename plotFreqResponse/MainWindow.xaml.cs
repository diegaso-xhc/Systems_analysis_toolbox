using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using Gigasoft.ProEssentials.Enums;
using differentialTerms;
using fftAndIfft;

namespace System_response_Toolbox
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public
            double dt;
            double tf;
            double x0;
            double[] input;
            double[] output;
            double t0 = 0;
            double[] x;
            double[] y;
            double fSInput; //Frequency for the sine input if chosen
            DifferentialEquations diff1;

        public MainWindow()
        {
            InitializeComponent();            
            //setupGraphPego();
        }

        private void setupGraphPego()
        {
            //double[] x = new double[2];
            //double[] y = new double[3];
            //x[0] = 1; x[1] = 1;
            //y[0] = 1; y[1] = 1; y[2] = 1;            
            this.diff1 = null;
            this.diff1 = new DifferentialEquations(this.y, this.x);            
            this.diff1.getTransferFunction();
            this.diff1.getStateSpace();            
            //double tf = 10;
            //double dt = 0.001;
           
            double xp0 = 0.0;                        
            
            
            //for (int i = 0; i < (int)input.Length; i++)
            //    input[i] = 1 * Math.Sin(2 * 3.14159 * 1 * (t0 + i * this.dt));
            //for (int i = 0; i < (int)input1.Length; i++)
            //{
            //    if (i == 0) input1[i] = 0;
            //    else input1[i] = 1;
            //}
            this.output = new double[(int)this.input.Length];
            //diff1.timeRespSim(dt, input1, x0, xp0, out output);
            this.diff1.rungeKutta4thO(this.dt, this.input, this.x0, xp0, out this.output);
            //diff1.rungeKutta4thO(dt, input1, input2, x0, xp0, out output);


            /*graph1.PeFunction.Reset();
            graph1.PeData.Subsets = 2;
            graph1.PeString.MainTitle = "Plotting Time Domain Response";
            graph1.PeString.SubTitle = "First Test";
            graph1.PeData.Y.Clear();
            graph1.PeData.X.Clear();
            graph1.PeData.NullDataValue = -99999;
            graph1.PeData.NullDataValueX = -99999;
            graph1.PePlot.Allow.Bar = false;
            graph1.PeFont.FontSize = Gigasoft.ProEssentials.Enums.FontSize.Medium;
            graph1.PeFont.SizeGlobalCntl = 1.2F;
            graph1.PeColor.SubsetColors[0] = Colors.Green;
            graph1.PeColor.SubsetColors[1] = Colors.Red;
            graph1.PePlot.SubsetLineTypes[0] = LineType.MediumThinSolid;
            graph1.PePlot.Allow.Line = true;
            graph1.PePlot.Method = SGraphPlottingMethod.Line;
            graph1.PeGrid.Configure.AutoPadBeyondZeroTX = true;
            graph1.PeGrid.Configure.AutoPadBeyondZeroX = true;
            graph1.PeGrid.Configure.ManualScaleControlX = ManualScaleControl.MinMax;
            graph1.PeGrid.Configure.ManualMinX = 0;
            graph1.PeGrid.Configure.ManualMaxX = this.tf;
            graph1.PeGrid.Configure.ManualScaleControlY = ManualScaleControl.MinMax;
            //graph1.PeGrid.Configure.ManualMinY = -2;
            //graph1.PeGrid.Configure.ManualMaxY = 2;
            if (this.input.Max() > this.output.Max()) graph1.PeGrid.Configure.ManualMaxY = (float)this.input.Max() + 1.0f;
            else graph1.PeGrid.Configure.ManualMaxY = (float)this.output.Max() + 1.0f;
            if (this.input.Min() < this.output.Min()) graph1.PeGrid.Configure.ManualMinY = (float)this.input.Min() - 1.0f;
            else graph1.PeGrid.Configure.ManualMinY = (float)this.output.Min() - 1.0f;
            int n = (int)(this.tf / this.dt);
            graph1.PeData.Points = n;
            for (int i = 0; i < n; i++)
            {
                graph1.PeData.X[0, i] = (float)(i * this.dt);
                graph1.PeData.Y[0, i] = (float)this.input[i];
                graph1.PeData.X[1, i] = (float)(i * this.dt);
                graph1.PeData.Y[1, i] = (float)this.output[i];
            }
            
            

            ////////////////////////This section is for plotting sine waves while testing the library////////////////////
            //int minVal = -4;
            //int maxVal = 4;
            //float gap = 0.001f;
            //int n = (int)((maxVal - minVal) / gap);
            //graph1.PeData.Points = n;
            //for (int i = 0; i < n; i++)
            //{
            //    graph1.PeData.X[0, i] = minVal + i * gap;
            //    graph1.PeData.Y[0, i] = (float)Math.Sin(2 * 3.1416 * (minVal + i * gap));
            //    graph1.PeData.X[1, i] = minVal + i * gap;
            //    graph1.PeData.Y[1, i] = (float)Math.Cos(2 * 3.1416 * 2*(minVal + i * gap));
            //}
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
            graph1.PeFunction.ReinitializeResetImage();            
            graph1.Invalidate();      */      

        }

        
        private void plotTimeResponse_Click(object sender, RoutedEventArgs e)
        {
            setupGraphPego();
        }

        private void sampleTime_TextChanged(object sender, TextChangedEventArgs e)
        {

        }

        private void finalTime_TextChanged(object sender, TextChangedEventArgs e)
        {

        }

        private void Initial_Condition_TextChanged(object sender, TextChangedEventArgs e)
        {

        }

        private void stepResponse_Click(object sender, RoutedEventArgs e)
        {
            this.dt = Convert.ToDouble(sampleTime.Text);
            this.tf = Convert.ToDouble(finalTime.Text);
            this.x0 = Convert.ToDouble(Initial_Condition.Text);
            this.input = null;
            this.input = new double[(int)(this.tf / this.dt)];
            for (int i = 0; i < (int)input.Length; i++)
            {
                if (i == 0) this.input[i] = 0;
                else this.input[i] = 1;
            }
        }

        private void sineResponse_Click(object sender, RoutedEventArgs e)
        {
            this.dt = Convert.ToDouble(sampleTime.Text);
            this.tf = Convert.ToDouble(finalTime.Text);
            this.x0 = Convert.ToDouble(Initial_Condition.Text);
            this.fSInput = Convert.ToDouble(fSineInp.Text);
            this.input = null;
            this.input = new double[(int)(this.tf / this.dt)];
            for (int i = 0; i < (int)this.input.Length; i++)
                this.input[i] = 1 * Math.Sin(2 * Math.PI * fSInput * (this.t0 + i * this.dt));
        }

        private void plotFrequencySpectrumFcn()
        {
            int npoints = (int)(this.tf / this.dt);            
            double[] amplitude = new double[npoints / 2];
            double[] freq = new double[npoints / 2];            
            //for (int i = 0; i < npoints; i++)
            //    this.input[i] = Math.Sin(2 * Math.PI * 100 * i * this.dt) + Math.Sin(2 * Math.PI * 340 * i * this.dt);
            var dft = fastFourierTransform.fFT(this.input, true);
            for (int i = 0; i < npoints; i += 2)
                amplitude[i / 2] = Math.Sqrt(dft[i] * dft[i] + dft[i + 1] * dft[i + 1]);
            for (int i = 0; i < (int)npoints / 2; i++)
                freq[i] = i / this.tf;
            /*
            freqPlot.PeFunction.Reset();
            freqPlot.PeData.Subsets = 1;
            freqPlot.PeString.MainTitle = "Plotting Frequency Domain Response";
            freqPlot.PeString.SubTitle = "First Test";
            freqPlot.PeData.Y.Clear();
            freqPlot.PeData.X.Clear();
            freqPlot.PeData.NullDataValue = -99999;
            freqPlot.PeData.NullDataValueX = -99999;
            freqPlot.PePlot.Allow.Bar = true;
            freqPlot.PeFont.FontSize = Gigasoft.ProEssentials.Enums.FontSize.Medium;
            freqPlot.PeFont.SizeGlobalCntl = 1.2F;
            freqPlot.PeColor.SubsetColors[0] = Colors.Green;
            freqPlot.PePlot.Method = SGraphPlottingMethod.Bar;
            freqPlot.PePlot.Method = SGraphPlottingMethod.Line;
            freqPlot.PePlot.Option.BarWidth = 0.001;
            freqPlot.PeGrid.Configure.AutoPadBeyondZeroTX = true;
            freqPlot.PeGrid.Configure.AutoPadBeyondZeroX = true;
            freqPlot.PeGrid.Configure.ManualScaleControlX = ManualScaleControl.MinMax;
            freqPlot.PeGrid.Configure.ManualMinX = 0;
            freqPlot.PeGrid.Configure.ManualMaxX = npoints / (2 * this.tf);
            freqPlot.PeGrid.Configure.ManualScaleControlY = ManualScaleControl.MinMax;
            freqPlot.PeGrid.Configure.ManualMaxY = (float)amplitude.Max() + 1.0f;
            freqPlot.PeGrid.Configure.ManualMinY = 0;
            int n = (int)amplitude.Length;
            freqPlot.PeData.Points = n;
            for (int i = 0; i < n; i++)
            {
                freqPlot.PeData.X[0, i] = (float)(freq[i]);
                freqPlot.PeData.Y[0, i] = (float)(amplitude[i]);
            }
            freqPlot.PeFunction.ReinitializeResetImage();
            freqPlot.Invalidate();   */
        }

        private void plotFrequencySpectrum_Click(object sender, RoutedEventArgs e)
        {
            plotFrequencySpectrumFcn();
        }

        private void plotBodeDiagram(double fini, double fend, double shiftT)
        {
            this.dt = Convert.ToDouble(sampleTime.Text);
            this.tf = Convert.ToDouble(finalTime.Text);
            this.x0 = Convert.ToDouble(Initial_Condition.Text);
            this.input = new double[(int)(this.tf / this.dt)];
            double k = (fend - fini) / shiftT;
            k = k / 2;
            for (int i = 0; i < (int)input.Length; i++)
                input[i] = 1 * Math.Sin(2 * Math.PI * (fini * i * this.dt + k * (i * this.dt) * (i * this.dt)));

            this.diff1 = null;
            this.diff1 = new DifferentialEquations(this.y, this.x);            
            this.diff1.getTransferFunction();
            this.diff1.getStateSpace();            
            double xp0 = 0.0;            
            this.output = new double[(int)this.input.Length];
            this.diff1.rungeKutta4thO(this.dt, this.input, this.x0, xp0, out this.output);
            int npoints = (int)(this.tf / this.dt);
            double[] amplitudeIn = new double[npoints / 2];
            double[] amplitudeOut = new double[npoints / 2];
            double[] ampBode = new double[npoints / 2];
            double[] freq = new double[npoints / 2];

            double[] phaseIn = new double[npoints / 2];
            double[] phaseOut = new double[npoints / 2];

            var dft = fastFourierTransform.fFT(this.input, true);
            for (int i = 0; i < npoints; i += 2)
            {
                amplitudeIn[i / 2] = Math.Sqrt(dft[i] * dft[i] + dft[i + 1] * dft[i + 1]);
                phaseIn[i / 2] = backupFunctions.angle(dft[i], dft[i + 1]); 
            }
            for (int i = 0; i < (int)npoints / 2; i++)
                freq[i] = i / (this.tf);
            dft = null;
            dft = fastFourierTransform.fFT(this.output, true);
            for (int i = 0; i < npoints; i += 2)
            {
                amplitudeOut[i / 2] = Math.Sqrt(dft[i] * dft[i] + dft[i + 1] * dft[i + 1]);
                phaseOut[i / 2] = backupFunctions.angle(dft[i], dft[i + 1]); 
            }
            
            double[] phaseInUnw = new double[npoints / 2];
            double[] phaseOutUnw = new double[npoints / 2];
            double[] phaseDiff = new double[npoints / 2];
            backupFunctions.unWrap(phaseIn, out phaseInUnw, Math.PI);
            backupFunctions.unWrap(phaseOut, out phaseOutUnw, Math.PI);

            for (int i = 0; i < (int)ampBode.Length; i++)
            {
                ampBode[i] = 20 * Math.Log10(amplitudeOut[i] / amplitudeIn[i]);
                phaseDiff[i] = (phaseOutUnw[i] - phaseInUnw[i]) * 180 / Math.PI;
            }

            /*
            freqPlot.PeFunction.Reset();
            freqPlot.PeString.MainTitle = "Plotting Frequency Domain Response";
            freqPlot.PeString.SubTitle = "Second Test";
            freqPlot.PeGrid.Configure.YAxisScaleControl = Gigasoft.ProEssentials.Enums.ScaleControl.Normal;
            freqPlot.PeGrid.Configure.XAxisScaleControl = Gigasoft.ProEssentials.Enums.ScaleControl.Log;
            freqPlot.PeColor.BitmapGradientMode = true;
            // logRight.PeColor.QuickStyle = QuickStyle.DarkNoBorder;
            freqPlot.PeGrid.LineControl = GridLineControl.Both;
            freqPlot.PeFont.FontSize = Gigasoft.ProEssentials.Enums.FontSize.Medium;
            freqPlot.PePlot.Method = SGraphPlottingMethod.Area;


            

            freqPlotPhase.PeFunction.Reset();
            freqPlotPhase.PeString.MainTitle = "Phase Diagram";
            freqPlotPhase.PeString.SubTitle = "Second Test";
            freqPlotPhase.PeGrid.Configure.YAxisScaleControl = Gigasoft.ProEssentials.Enums.ScaleControl.Normal;
            freqPlotPhase.PeGrid.Configure.XAxisScaleControl = Gigasoft.ProEssentials.Enums.ScaleControl.Log;
            freqPlotPhase.PeColor.BitmapGradientMode = true;
            // logRight.PeColor.QuickStyle = QuickStyle.DarkNoBorder;
            int n = (int)ampBode.Length/3;
            freqPlotPhase.PeGrid.LineControl = GridLineControl.Both;
            freqPlotPhase.PeFont.FontSize = Gigasoft.ProEssentials.Enums.FontSize.Medium;
            freqPlotPhase.PePlot.Method = SGraphPlottingMethod.Area;
            freqPlotPhase.PeGrid.Configure.ManualMinX = (float)freq.Min();
            freqPlotPhase.PeGrid.Configure.ManualMaxX = (float)freq[n - 1];
            freqPlotPhase.PeGrid.Configure.ManualScaleControlY = ManualScaleControl.MinMax;
            freqPlotPhase.PeGrid.Configure.ManualMaxY = (float)phaseDiff.Max() + 5.0f;
            freqPlotPhase.PeGrid.Configure.ManualMinY = (float)phaseDiff.Min() - 5.0f;
            // Change data so it varies over wider range // 
            
            freqPlot.PeData.Points = n;
            freqPlotPhase.PeData.Points = n;
            for (int i = 0; i < n/2; i++)
            {
                freqPlot.PeData.X[0, i] = (float)(freq[i]*2*Math.PI);
                freqPlot.PeData.Y[0, i] = (float)(ampBode[i]);
                freqPlotPhase.PeData.X[0, i] = (float)(freq[i] * 2 * Math.PI);
                freqPlotPhase.PeData.Y[0, i] = (float)(phaseDiff[i]);
            }

            freqPlot.PePlot.Method = SGraphPlottingMethod.Line;
            freqPlot.PePlot.MarkDataPoints = false;
            freqPlot.PeLegend.SubsetPointTypes[0] = PointType.DotSolid;
            freqPlot.PeLegend.SubsetPointTypes[1] = PointType.DotSolid;
            freqPlot.PeString.XAxisLabel = "Frequency in [rad / s]";
            
            freqPlotPhase.PePlot.Method = SGraphPlottingMethod.Line;
            freqPlotPhase.PePlot.MarkDataPoints = false;
            freqPlotPhase.PeLegend.SubsetPointTypes[0] = PointType.DotSolid;
            freqPlotPhase.PeLegend.SubsetPointTypes[1] = PointType.DotSolid;
            freqPlotPhase.PeString.XAxisLabel = "Frequency in [rad / s]";

            freqPlot.PeFunction.ReinitializeResetImage();
            freqPlotPhase.PeFunction.ReinitializeResetImage();
            freqPlot.Invalidate();
            freqPlotPhase.Invalidate();*/
        }

        private void initial_frequency_TextChanged(object sender, TextChangedEventArgs e)
        {

        }

        private void final_frequency_TextChanged(object sender, TextChangedEventArgs e)
        {

        }

        private void shift_time_TextChanged(object sender, TextChangedEventArgs e)
        {

        }

        private void plotBode_Click(object sender, RoutedEventArgs e)
        {
            double f0 = Convert.ToDouble(initial_frequency.Text);
            double f1 = Convert.ToDouble(final_frequency.Text);
            double shiftT = Convert.ToDouble(shift_time.Text);
            plotBodeDiagram(f0, f1, shiftT);
        }

        private void createModel_Click(object sender, RoutedEventArgs e)
        {
            this.x = null;
            this.y = null;
            // Create a FlowDocument to contain content for the RichTextBox.            
            char[] delimiterChars = {'[',',',' ',']','\r','\n'};
            TextRange textRange = new TextRange(numTF.Document.ContentStart, numTF.Document.ContentEnd);
            string text = textRange.Text;
            string[] data = text.Split(delimiterChars);
            data = data.Where(xt => !string.IsNullOrEmpty(xt)).ToArray();
            int nData = 0;            
            for (int i = 0; i < data.Length; i++)
            {
                if (data[i] != ";") nData++;
                else break;
            }
            int dData = (data.Length - 1) - nData; 
            this.x = new double[nData];
            this.y = new double[dData];
            for (int i = 0; i < data.Length; i++)
            {
                if (i < nData) x[i] = Convert.ToDouble(data[i]);
                else if (i > nData ) y[i - nData - 1] = Convert.ToDouble(data[i]);
            }           
        }

        private void numTF_TextChanged(object sender, TextChangedEventArgs e)
        {

        }

        private void fSineInp_TextChanged(object sender, TextChangedEventArgs e)
        {

        }

        private void thoericBode(double wmax, double prn) //max frequency in rad/s and the precision
        { 

            int nL = (int)this.x.Length;
            int dL = (int)this.y.Length;
            double[] nRe;
            double[] nIm;
            if (nL >= 2)
            {
                nRe = new double[nL - 1]; //Getting the number of real roots we will get for the numerator of the Transfer Function
                nIm = new double[nL - 1]; //Getting the number of imaginary roots we will get for the numerator of the Transfer Function
            }
            else
            {
                nRe = new double[1];
                nIm = new double[1];
            }
            double[] dRe;
            double[] dIm;
            if (dL >= 2)
            {
                dRe = new double[dL - 1]; //Getting the number of real roots we will get for the denominator of the Transfer Function
                dIm = new double[dL - 1]; //Getting the number of imaginary roots we will get for the denominator of the Transfer Function
            }
            else 
            {
                dRe = new double[1];
                dIm = new double[1];
            }
            
            if (nL >= 2)
            {
                double px = this.x[0];
                for (int i = 0; i < (int)this.x.Length; i++)
                    this.x[i] = this.x[i] / px;
                backupFunctions.getRoots(this.x, out nRe, out nIm);
                for (int i = 0; i < (int)this.x.Length; i++)
                    this.x[i] = this.x[i] * px;                
            }
            else nRe[0] = this.x[0];
            if (dL >= 2)
            {
                double py = this.y[0];
                for (int i = 0; i < (int)this.y.Length; i++)
                    this.y[i] = this.y[i] / py;
                backupFunctions.getRoots(this.y, out dRe, out dIm);
                for (int i = 0; i < (int)this.y.Length; i++)
                    this.y[i] = this.y[i] * py;
            }
            else dRe[0] = this.y[0];
            int wL = (int)(wmax / prn);
            double[] nMag = new double[wL];
            double[] nPh = new double[wL];
            double[] mag = new double[wL];
            double[] ph = new double[wL];
            for (int i = 0; i < wL; i++)
            {
                nMag[i] = this.x[0]; 
                mag[i] = this.y[0];
                if (this.x[0] >= 0) nPh[i] = 0;
                else nPh[i] = Math.PI;
                if (this.y[0] >= 0) ph[i] = 0;
                else ph[i] = Math.PI;
            }
            if (nL >= 2)
            {
                int i = 0;
                for (double w = 0; w < wmax - prn; w += prn)
                {
                    for (int j = 0; j < nL - 1; j++)
                    {
                        nMag[i] *= Math.Sqrt(nRe[j] * nRe[j] + (w - nIm[j]) * (w - nIm[j]));
                        nPh[i] += backupFunctions.angle(-1.0*nRe[j], (w - nIm[j]));
                    }
                    i++;
                }
            }
            else
            {
                int i = 0;
                for (double w = 0; w < wmax - prn; w += prn)
                {
                    nMag[i] = this.x[0]; 
                    i++;
                }
            }
            if (dL >= 2)
            {
                int i = 0;
                for (double w = 0; w < wmax - prn; w += prn)
                {
                    for (int j = 0; j < dL - 1; j++)
                    {
                        mag[i] *= Math.Sqrt(dRe[j] * dRe[j] + (w - dIm[j]) * (w - dIm[j]));
                        ph[i] += backupFunctions.angle(-1.0 * dRe[j], (w - dIm[j]));
                    }
                    mag[i] = nMag[i] / mag[i];
                    ph[i] = nPh[i] - ph[i];
                    ph[i] = ph[i] * 180 / Math.PI;
                    if (Math.Abs(ph[i]) > 360)
                    {
                        if (ph[i] > 0) ph[i] -= 360;
                        else ph[i] += 360;
                    }
                    i++;                    
                }
            }
            else
            {
                int i = 0;
                for (double w = 0; w < wmax - prn; w += prn)
                {
                    mag[i] = nMag[i] / this.y[0]; 
                    i++;                    
                }
            }
                
            double[] freq = new double[wL];
            for (int i = 0; i < wL; i++)
                freq[i] = prn * i;

            /*
            freqPlot.PeFunction.Reset();
            freqPlot.PeString.MainTitle = "Plotting Frequency Domain Response";
            freqPlot.PeString.SubTitle = "Second Test";
            freqPlot.PeGrid.Configure.YAxisScaleControl = Gigasoft.ProEssentials.Enums.ScaleControl.Normal;
            freqPlot.PeGrid.Configure.XAxisScaleControl = Gigasoft.ProEssentials.Enums.ScaleControl.Log;
            freqPlot.PeColor.BitmapGradientMode = true;
            // logRight.PeColor.QuickStyle = QuickStyle.DarkNoBorder;
            freqPlot.PeGrid.LineControl = GridLineControl.Both;
            freqPlot.PeFont.FontSize = Gigasoft.ProEssentials.Enums.FontSize.Medium;
            freqPlot.PePlot.Method = SGraphPlottingMethod.Area;

            freqPlotPhase.PeFunction.Reset();
            freqPlotPhase.PeString.MainTitle = "Phase Diagram";
            freqPlotPhase.PeString.SubTitle = "Second Test";
            freqPlotPhase.PeGrid.Configure.YAxisScaleControl = Gigasoft.ProEssentials.Enums.ScaleControl.Normal;
            freqPlotPhase.PeGrid.Configure.XAxisScaleControl = Gigasoft.ProEssentials.Enums.ScaleControl.Log;
            freqPlotPhase.PeColor.BitmapGradientMode = true;
            // logRight.PeColor.QuickStyle = QuickStyle.DarkNoBorder;
            freqPlotPhase.PeGrid.LineControl = GridLineControl.Both;
            freqPlotPhase.PeFont.FontSize = Gigasoft.ProEssentials.Enums.FontSize.Medium;
            freqPlotPhase.PePlot.Method = SGraphPlottingMethod.Area;
            //freqPlotPhase.PeGrid.Configure.ManualMinX = (float)ph.Min();
            //freqPlotPhase.PeGrid.Configure.ManualMaxX = (float)ph.Max();
            //freqPlotPhase.PeGrid.Configure.ManualScaleControlY = ManualScaleControl.MinMax;
            //freqPlotPhase.PeGrid.Configure.ManualMaxY = (float)wmax;
            //freqPlotPhase.PeGrid.Configure.ManualMinY = (float)0;
            // Change data so it varies over wider range // 

            freqPlot.PeData.Points = wL;            
            freqPlotPhase.PeData.Points = wL;

            for (int i = 0; i < wL; i++)
            {
                freqPlot.PeData.X[0, i] = (float)(freq[i]);
                freqPlot.PeData.Y[0, i] = (float)(20*Math.Log10(mag[i]));
                freqPlotPhase.PeData.X[0, i] = (float)(freq[i]);
                freqPlotPhase.PeData.Y[0, i] = (float)(ph[i]);
            }

            freqPlot.PePlot.Method = SGraphPlottingMethod.Line;
            freqPlot.PePlot.MarkDataPoints = false;
            freqPlot.PeLegend.SubsetPointTypes[0] = PointType.DotSolid;
            freqPlot.PeLegend.SubsetPointTypes[1] = PointType.DotSolid;
            freqPlot.PeString.XAxisLabel = "Frequency in [rad / s]";

            freqPlotPhase.PePlot.Method = SGraphPlottingMethod.Line;
            freqPlotPhase.PePlot.MarkDataPoints = false;
            freqPlotPhase.PeLegend.SubsetPointTypes[0] = PointType.DotSolid;
            freqPlotPhase.PeLegend.SubsetPointTypes[1] = PointType.DotSolid;
            freqPlotPhase.PeString.XAxisLabel = "Frequency in [rad / s]";

            freqPlotPhase.PeFunction.ReinitializeResetImage();
            freqPlot.PeFunction.ReinitializeResetImage();            
            freqPlotPhase.Invalidate();
            freqPlot.Invalidate();*/
            
        }

        private void getNichols(double wmax, double prn)
        {
            int nL = (int)this.x.Length;
            int dL = (int)this.y.Length;
            double[] nRe;
            double[] nIm;
            if (nL >= 2)
            {
                nRe = new double[nL - 1]; //Getting the number of real roots we will get for the numerator of the Transfer Function
                nIm = new double[nL - 1]; //Getting the number of imaginary roots we will get for the numerator of the Transfer Function
            }
            else
            {
                nRe = new double[1];
                nIm = new double[1];
            }
            double[] dRe;
            double[] dIm;
            if (dL >= 2)
            {
                dRe = new double[dL - 1]; //Getting the number of real roots we will get for the denominator of the Transfer Function
                dIm = new double[dL - 1]; //Getting the number of imaginary roots we will get for the denominator of the Transfer Function
            }
            else
            {
                dRe = new double[1];
                dIm = new double[1];
            }

            if (nL >= 2)
            {
                double px = this.x[0];
                for (int i = 0; i < (int)this.x.Length; i++)
                    this.x[i] = this.x[i] / px;
                backupFunctions.getRoots(this.x, out nRe, out nIm);
                for (int i = 0; i < (int)this.x.Length; i++)
                    this.x[i] = this.x[i] * px;
            }
            else nRe[0] = this.x[0];
            if (dL >= 2)
            {
                double py = this.y[0];
                for (int i = 0; i < (int)this.y.Length; i++)
                    this.y[i] = this.y[i] / py;
                backupFunctions.getRoots(this.y, out dRe, out dIm);
                for (int i = 0; i < (int)this.y.Length; i++)
                    this.y[i] = this.y[i] * py;
            }
            else dRe[0] = this.y[0];
            int wL = (int)(wmax / prn);
            double[] nMag = new double[wL];
            double[] nPh = new double[wL];
            double[] mag = new double[wL];
            double[] ph = new double[wL];
            for (int i = 0; i < wL; i++)
            {
                nMag[i] = this.x[0];
                mag[i] = this.y[0];
                if (this.x[0] >= 0) nPh[i] = 0;
                else nPh[i] = Math.PI;
                if (this.y[0] >= 0) ph[i] = 0;
                else ph[i] = Math.PI;
            }
            if (nL >= 2)
            {
                int i = 0;
                for (double w = 0; w < wmax - prn; w += prn)
                {
                    for (int j = 0; j < nL - 1; j++)
                    {
                        nMag[i] *= Math.Sqrt(nRe[j] * nRe[j] + (w - nIm[j]) * (w - nIm[j]));
                        nPh[i] += backupFunctions.angle(-1.0 * nRe[j], (w - nIm[j]));
                    }
                    i++;
                }
            }
            else
            {
                int i = 0;
                for (double w = 0; w < wmax - prn; w += prn)
                {
                    nMag[i] = this.x[0];
                    i++;
                }
            }
            if (dL >= 2)
            {
                int i = 0;
                for (double w = 0; w < wmax - prn; w += prn)
                {
                    for (int j = 0; j < dL - 1; j++)
                    {
                        mag[i] *= Math.Sqrt(dRe[j] * dRe[j] + (w - dIm[j]) * (w - dIm[j]));
                        ph[i] += backupFunctions.angle(-1.0 * dRe[j], (w - dIm[j]));
                    }
                    mag[i] = nMag[i] / mag[i];
                    ph[i] = nPh[i] - ph[i];
                    ph[i] = ph[i] * 180 / Math.PI;
                    if (Math.Abs(ph[i]) > 360)
                    {
                        if (ph[i] > 0) ph[i] -= 360;
                        else ph[i] += 360;
                    }
                    i++;
                }
            }
            else
            {
                int i = 0;
                for (double w = 0; w < wmax - prn; w += prn)
                {
                    mag[i] = nMag[i] / this.y[0];
                    i++;
                }
            }

            double[] freq = new double[wL];
            for (int i = 0; i < wL; i++)
                freq[i] = prn * i;
            /*
            graph1.PeFunction.Reset();
            graph1.PeString.MainTitle = "Plotting Frequency Domain Response";
            graph1.PeString.SubTitle = "Second Test";
            graph1.PeGrid.Configure.YAxisScaleControl = Gigasoft.ProEssentials.Enums.ScaleControl.Normal;
            graph1.PeGrid.Configure.XAxisScaleControl = Gigasoft.ProEssentials.Enums.ScaleControl.Normal;
            graph1.PeColor.BitmapGradientMode = true;
            // logRight.PeColor.QuickStyle = QuickStyle.DarkNoBorder;
            graph1.PeGrid.LineControl = GridLineControl.Both;
            graph1.PeFont.FontSize = Gigasoft.ProEssentials.Enums.FontSize.Medium;
            graph1.PePlot.Method = SGraphPlottingMethod.Area;

            graph1.PeData.Points = wL;
            
            for (int i = 0; i < wL; i++)
            {
                graph1.PeData.X[0, i] = (float)(ph[i]);
                graph1.PeData.Y[0, i] = (float)(20 * Math.Log10(mag[i]));            
            }

            graph1.PePlot.Method = SGraphPlottingMethod.Line;
            graph1.PePlot.MarkDataPoints = false;
            graph1.PeLegend.SubsetPointTypes[0] = PointType.DotSolid;
            graph1.PeLegend.SubsetPointTypes[1] = PointType.DotSolid;
            graph1.PeString.XAxisLabel = "Phase in [degree]";

            graph1.PeFunction.ReinitializeResetImage();
            graph1.Invalidate();*/
        }
        
        private void plotTheoricBode_Click(object sender, RoutedEventArgs e)
        {
            double wM = Convert.ToDouble(wMax.Text);
            double p = Convert.ToDouble(precision.Text);
            thoericBode(wM, p);            
        }

        private void TextBox_TextChanged_1(object sender, TextChangedEventArgs e)
        {

        }

        private void precision_TextChanged(object sender, TextChangedEventArgs e)
        {

        }

        private void nicholsChart_Click(object sender, RoutedEventArgs e)
        {
            double wM = Convert.ToDouble(wMax.Text);
            double p = Convert.ToDouble(precision.Text);
            getNichols(wM, p);            
        }
    }
}
