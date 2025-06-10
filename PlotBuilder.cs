using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using System.Numerics;

namespace GAZDIN_CALC_APP
{
	/// <summary>
	/// Use this static class for building graph
	/// </summary>
	public class PlotBuilder
	{
		private enum FuncType
		{
			T,
			P,
			RO,
			M,
			w
		}

		public PlotBuilder()
		{
			Plot = new PlotModel() { Title = "test" };
			ResultPlot = new PlotModel() { Title = "Результат расчетов" };
			
			var linearAxis1 = new LinearAxis();
			linearAxis1.Title = "Давление";
			linearAxis1.TitleColor = OxyColor.FromRgb(0, 0, 255);
			linearAxis1.PositionTier = 1;
			linearAxis1.AxisDistance = 25;
			linearAxis1.Key = "P";
			linearAxis1.AxislineColor = OxyColor.FromRgb(0, 0, 255);
			var linearAxis2 = new LinearAxis();
			linearAxis2.Title = "Температура";
			linearAxis2.TitleColor = OxyColor.FromRgb(255, 0, 0);
			linearAxis2.PositionTier = 2;
			linearAxis2.AxisDistance = 25;
			linearAxis2.Key = "T";
			linearAxis2.AxislineColor = OxyColor.FromRgb(255, 0, 0);
			var linearAxis3 = new LinearAxis();
			linearAxis3.Title = "Скорость потока";
			linearAxis3.TitleColor = OxyColor.FromRgb(0, 255, 0);
			linearAxis3.PositionTier = 3;
			linearAxis3.AxisDistance = 25;
			linearAxis3.Key = "w";
			linearAxis3.AxislineColor = OxyColor.FromRgb(0, 255, 0);
			var linearAxis4 = new LinearAxis();
			linearAxis4.Title = "Плотность";
			linearAxis4.TitleColor = OxyColor.FromRgb(100, 100, 100);
			linearAxis4.PositionTier = 4;
			linearAxis4.AxisDistance = 25;
			linearAxis4.Key = "RO";
			linearAxis4.AxislineColor = OxyColor.FromRgb(100, 100, 100);
			var linearAxis5 = new LinearAxis();
			linearAxis5.Title = "Число Маха";
			linearAxis5.TitleColor = OxyColor.FromRgb(0, 125, 255);
			linearAxis5.PositionTier = 5;
			linearAxis5.AxisDistance = 25;
			linearAxis5.Key = "M";
			linearAxis5.AxislineColor = OxyColor.FromRgb(0, 125, 255);

			ResultPlot.Axes.Add(linearAxis1);
			ResultPlot.Axes.Add(linearAxis2);
			ResultPlot.Axes.Add(linearAxis3);
			ResultPlot.Axes.Add(linearAxis4);
			ResultPlot.Axes.Add(linearAxis5);
		}

		/// <summary>
		/// Func for building line
		/// </summary>
		/// <param name="start">The start point of line</param>
		/// <param name="end">The end point of line</param>
		public static void BuildLine(Vector2 start, Vector2 end)
		{
			LineSeries ls = new LineSeries();
			ls.Points.Add(new DataPoint(start.X, start.Y));
			ls.Points.Add(new DataPoint(end.X, end.Y));
			ls.Color = OxyColor.FromRgb(135, 145, 100);

			Plot.Series.Add(ls);
		}

		public static void BuildResults(List<double> T, List<double> P, List<double> RO, List<double> M, List<double> w, List<double> xs)
		{
			T.RemoveAll(double.IsNaN);
			P.RemoveAll(double.IsNaN);
			RO.RemoveAll(double.IsNaN);
			M.RemoveAll(double.IsNaN);
			w.RemoveAll(double.IsNaN);

			for(int i = 0; i < T.Count; i++)
			{
				if (i != T.Count - 1)
				{
					BuildResultLine(new Vector2((float)xs[i], (float)T[i]), new Vector2((float)xs[i + 1], (float)T[i + 1]), FuncType.T);
					BuildResultLine(new Vector2((float)xs[i], (float)P[i]), new Vector2((float)xs[i + 1], (float)P[i + 1]), FuncType.P);
					BuildResultLine(new Vector2((float)xs[i], (float)RO[i]), new Vector2((float)xs[i + 1], (float)RO[i + 1]), FuncType.RO);
					BuildResultLine(new Vector2((float)xs[i], (float)M[i]), new Vector2((float)xs[i + 1], (float)M[i + 1]), FuncType.M);
					BuildResultLine(new Vector2((float)xs[i], (float)w[i]), new Vector2((float)xs[i + 1], (float)w[i + 1]), FuncType.w);
				}
				else
					break;
			}
		}

		private static void BuildResultLine(Vector2 start, Vector2 end, FuncType type)
		{
			LineSeries ls = new LineSeries();
			switch (type)
			{
				case FuncType.T:
					ls.Points.Add(new DataPoint(start.X, start.Y));
					ls.Points.Add(new DataPoint(end.X, end.Y));
					ls.YAxisKey = "T";
					ls.Color = OxyColor.FromRgb(255, 0, 0);
					ResultPlot.Series.Add(ls);

					break;
				case FuncType.P:
					ls = new LineSeries();
					ls.Points.Add(new DataPoint(start.X, start.Y));
					ls.Points.Add(new DataPoint(end.X, end.Y));
					ls.YAxisKey = "P";
					ls.Color = OxyColor.FromRgb(0, 0, 255);
					ResultPlot.Series.Add(ls);

					break;
				case FuncType.RO:
					ls.Points.Add(new DataPoint(start.X, start.Y));
					ls.Points.Add(new DataPoint(end.X, end.Y));
					ls.YAxisKey = "RO";
					ls.Color = OxyColor.FromRgb(100, 100, 100);
					ResultPlot.Series.Add(ls);

					break;
				case FuncType.M:
					ls.Points.Add(new DataPoint(start.X, start.Y));
					ls.Points.Add(new DataPoint(end.X, end.Y));
					ls.YAxisKey = "M";
					ls.Color = OxyColor.FromRgb(0, 125, 255);
					ResultPlot.Series.Add(ls);

					break;
				case FuncType.w:
					ls.Points.Add(new DataPoint(start.X, start.Y));
					ls.Points.Add(new DataPoint(end.X, end.Y));
					ls.YAxisKey = "w";
					ls.Color = OxyColor.FromRgb(0, 255, 0);
					ResultPlot.Series.Add(ls);

					break;
			}
		}

		public static void BuildCurve(Func<double, double> func, double start, double end, double shag)
		{
			var fs = new FunctionSeries(func, start, end, shag);
			Plot.Series.Add(fs);
		}

		/// <summary>
		/// Refreshing the graph
		/// </summary>
		public static void Refresh()
		{
			Plot.InvalidatePlot(true);
			ResultPlot.InvalidatePlot(true);
		}

		/// <summary>
		/// Clears the graph
		/// </summary>
		public static void Clear()
		{
			Plot.Series.Clear();
			Refresh();
		}

		public static PlotModel Plot { get; private set; }
		public static PlotModel ResultPlot { get; private set; }
	}
}
