using GAZDIN_CALC_APP.SplineAPI;
using System.IO;
using System.Numerics;
using System.Windows;

namespace GAZDIN_CALC_APP
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private double w_www;

        public MainWindow()
        {
            InitializeComponent();

            PlotBuilder.OnGraphDraw += ParseStrings;
        }

        //private Func<double, double> _calculationFunc;

		private void Button_Click(object sender, RoutedEventArgs e)
		{
			PlotBuilder.Clear();

            if ((bool)IsAnalytic.IsChecked)
            {
				ParseStrings(out Vector2 start, out Vector2 end, out Vector2 crit);

				double[] x = { start.X, crit.X, end.X };
				double[] y = { start.Y, crit.Y, end.Y };
				double[] derrivatives = { 0, 0, 0 };

				Func<double, double> func = HermiteSpline.Interpolate(x, y, derrivatives);
				PlotBuilder.BuildCurve(func, start.X, end.X, (end.X - start.X) / 100);
				PlotBuilder.BuildLine(start, new Vector2(start.X, -start.Y));
				PlotBuilder.BuildLine(end, new Vector2(end.X, -end.Y));

				y = new double[] { -start.Y, -crit.Y, -end.Y };
				func = HermiteSpline.Interpolate(x, y, derrivatives);
				PlotBuilder.BuildCurve(func, start.X, end.X, (end.X - start.X) / 100);
			}
            else
            {
                string path = Path.Text;
                var lines = File.ReadAllLines(path);

				double[] x = new double[lines.Length];
				double[] y = new double[lines.Length];

				for (int i = 0; i < lines.Length; i++)
                {
                    x[i] = double.Parse(lines[i].Split(' ')[0]) * 100;
					y[i] = double.Parse(lines[i].Split(' ')[2]) * 100;
				}

                double lDerrivative = 0;
                double rDerrivative = 0;

				Func<double, double> func = CubicSpline.Interpolate(x, y, lDerrivative, rDerrivative);

                PlotBuilder.BuildCurve(func, x[0], x[lines.Length - 1], (x[lines.Length - 1] - x[0]) / 180);
				PlotBuilder.BuildLine(new Vector2((float)x[0], (float)y[0]), new Vector2((float)x[0], -(float)y[0]));
				PlotBuilder.BuildLine(new Vector2((float)x[x.Length - 1], (float)y[y.Length - 1]), new Vector2((float)x[x.Length - 1], -(float)y[y.Length - 1]));

				for (int i = 0; i < lines.Length; i++)
				{
					y[i] = -double.Parse(lines[i].Split(' ')[2]) * 100;
				}

                func = CubicSpline.Interpolate(x, y, rDerrivative, lDerrivative);
				PlotBuilder.BuildCurve(func, x[0], x[lines.Length - 1], (x[lines.Length - 1] - x[0]) / 180);
			}

            PlotBuilder.Refresh();
		}

		private void ParseStrings(out Vector2 start, out Vector2 end, out Vector2 crit)
		{
            string line = startVec.Text;
            start = new Vector2(float.Parse(line.Split(", ")[0]), float.Parse(line.Split(", ")[1]));

            line = endVec.Text;
            end = new Vector2(float.Parse(line.Split(", ")[0]), float.Parse(line.Split(", ")[1]));

            line = critPoint.Text;
            crit = new Vector2(float.Parse(line.Split(", ")[0]), float.Parse(line.Split(", ")[1]));
		}
        private void ParseStrings(out float k, out float c, out float T_0, out float P_0, out float RO_0, out float epsilon, out float numOfParts)
        {
            k = float.Parse(this.k.Text);
            c = float.Parse(this.c.Text);
            T_0 = float.Parse(this.T_0.Text);
            P_0 = float.Parse(this.P_0.Text);
            RO_0 = float.Parse(this.RO_0.Text);
            epsilon = float.Parse(this.epsilon.Text);
            numOfParts = float.Parse(this.numOfParts.Text);
        }

        private void ParseStrings(out double T, out double P, out double RO, out bool isBase)
        {
			T = double.Parse(this.T_0.Text);
			P = double.Parse(this.P_0.Text);
			RO = double.Parse(this.RO_0.Text);

            isBase = (bool)IsBase.IsChecked;
        }

        private void Button_Click_1(object sender, RoutedEventArgs e)
		{
            PlotBuilder.Clear();
            
            ParseStrings(out float k, out float c, out float T_0, out float P_0, out float RO_0, out float epsilon, out float numOfParts);

            float d_crit = 0;
            double xLength = 0;
            double crit_x = 0;

            Func<double, double> func;

			if (!(bool)IsAnalytic.IsChecked)
			{
				var lines = File.ReadAllLines(Path.Text);
				double[] xss = new double[lines.Length];
                double[] ys = new double[lines.Length];

				for (int i = 0; i < lines.Length; i++)
				{
					xss[i] = double.Parse(lines[i].Split(' ')[0]);
                    ys[i] = double.Parse(lines[i].Split(' ')[2]);
				}

				d_crit = (float)ys.Min();

                float[] th = new float[ys.Length];
                for (int i = 0; i < ys.Length; i++) th[i] = (float)ys[i];
 
                int index = MathHelper.IndexOf<float>(th, d_crit);

				crit_x = xss[index];                

                xLength = xss.Last();

                float lDerrivative = 0, rDerrivative = 0;

                func = CubicSpline.Interpolate(xss, ys, lDerrivative, rDerrivative);
			}
            else
            {
                d_crit = float.Parse(critPoint.Text.Split(", ")[1]) / 100;
                xLength = double.Parse(endVec.Text.Split(", ")[0]) / 100;
                crit_x = double.Parse(critPoint.Text.Split(", ")[0]) / 100;

                ParseStrings(out Vector2 start, out Vector2 end, out Vector2 crit);
                start /= 100;
                end /= 100;
                crit /= 100;

                double[] xss = { start.X, crit.X, end.X };
                double[] ys = { start.Y, crit.Y, end.Y };
                double[] derrivatives = { 0, 0, 0 };

                func = HermiteSpline.Interpolate(xss, ys, derrivatives);
            }
            
			List<double> T = new List<double>();
			List<double> P = new List<double>();
			List<double> RO = new List<double>();
			List<double> M = new List<double>();
			List<double> w = new List<double>();
            List<double> xs = new List<double>();

			for (double x = 0; x < xLength; x += (xLength / numOfParts))
            {
                float d = (float)func(x);

                bool isOverSound = x < crit_x ? false : true;
                
                Dictionary<string, double> results = MathAPI.CalculatePart(k, epsilon, d_crit, T_0, P_0, RO_0, d, c, isOverSound);

                if(results == null)
                {
                    continue;
                }

                T.Add(results["T"]);
                P.Add(results["P"]);
                RO.Add(results["RO"]);
                M.Add(results["M"]);
                w.Add(results["w"]);
                xs.Add(x);
            }



            PlotBuilder.BuildResults(T, P, RO, M, w, xs);

            PlotBuilder.Refresh();
		}
	}
}