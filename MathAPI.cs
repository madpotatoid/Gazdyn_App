using System.Numerics;

namespace GAZDIN_CALC_APP
{
	public static class MathAPI
	{
		/// <summary>
		/// Считает газдин функции для конкретного сечения
		/// </summary>
		/// <param name="k">Аддиабата</param>
		/// <param name="epsilon">Точность</param>
		/// <param name="d_crit">Диаметр критического сечения</param>
		/// <param name="T_0">Температура торможения</param>
		/// <param name="P_0">Давление торможения</param>
		/// <param name="RO_0">Плотность торможения</param>
		/// <param name="d">Текущий диаметр</param>
		/// <param name="c">Скорость звука в среде</param>
		/// <param name="isOverSound">Число за критическим диаметром или нет</param>
		/// <returns>Словарь, где хранятся функции торможения, плотности давления, числа маха и скорости газа для конкретного сечения</returns>

		private static float _ideal_oversound_lambda = 0;

		public static Dictionary<string, double> CalculatePart(float k, float epsilon, float d_crit, float T_0, float P_0, float RO_0, float d, double c, bool isOverSound)
		{
			Dictionary<string, double> results = new Dictionary<string, double>();

			float q_labda = MathF.Pow(d_crit / d, 2);

			float T = 0;
			float P = 0;
			float RO = 0;

			float lambda = 0;
			
			if (!isOverSound)
			{
				float lambda_last = 0.02f;
				lambda = q_labda / (MathF.Pow((k + 1) / 2, 1 / (k - 1)) * MathF.Pow(1 - (k - 1) / (k + 1) * MathF.Pow(lambda_last, 2), 1 / (k + 1)));
				while (MathF.Abs(lambda - lambda_last) > epsilon)
				{
					lambda_last = lambda;
					float smth = MathF.Pow(lambda_last, 2);
					float arg1 = 1 - ((k - 1) / (k + 1)) * smth;
					float arg2 = 1 / (k + 1);
					float smth2 = MathF.Pow(arg1, arg2);
					float smth3 = MathF.Pow((k + 1) / 2, 1 / (k + 1));
					lambda = q_labda / (smth3 * smth2);
				}
			}
			else
			{
				float lambda_last = 1.01f;
				float arg1 = (k + 1) / (k - 1);
				float arg2 = MathF.Pow(lambda_last, k - 1) * (k + 1) / 2;
				float arg3 = MathF.Pow(q_labda, k - 1);
				float arg4 = arg3 / arg2;
				float arg5 = 1 - arg4;
				lambda = MathF.Sqrt(arg1 * arg5);

				//float lambda_last = 1.01f;
				//float arg1 = (k + 1) / (k - 1);
				//float arg2 = MathF.Pow(lambda_last, k + 1) * MathF.Pow((k + 1) / 2, (k + 1) / (k - 1));
				//float arg3 = q_labda;
				//float arg4 = arg3 / arg2;
				//float arg5 = 1 - arg4;
				//lambda = MathF.Sqrt(arg1 * arg5);

				while (MathF.Abs(lambda - lambda_last) > epsilon)
				{
					lambda_last = lambda;
					float arg1_1 = (k + 1) / (k - 1);
					float arg2_1 = MathF.Pow(lambda_last, k - 1) * (k + 1) / 2;
					float arg3_1 = MathF.Pow(q_labda, k - 1);
					float arg4_1 = arg3_1 / arg2_1;
					float arg5_1 = 1 - arg4_1;
					lambda = MathF.Sqrt(arg1_1 * arg5_1);


					//lambda_last = lambda;
					//float arg1_1 = (k + 1) / (k - 1);
					//float arg2_1 = MathF.Pow(lambda_last, k + 1) * MathF.Pow((k + 1) / 2, (k + 1) / (k - 1));
					//float arg3_1 = q_labda;
					//float arg4_1 = arg3_1 / arg2_1;
					//float arg5_1 = 1 - arg4_1;
					//lambda = MathF.Sqrt(arg1_1 * arg5_1);
				}
			}

			if (isOverSound)
			{
				Console.WriteLine(" ");
			}

			//if (_ideal_oversound_lambda > lambda && isOverSound) 
			//	return null;

			//_ideal_oversound_lambda = lambda;

			T = (1 - ((k - 1) / (k + 1) * MathF.Pow(lambda, 2))) * T_0;
			P = MathF.Pow((1 - ((k - 1) / (k + 1) * MathF.Pow(lambda, 2))), k / (k - 1)) * P_0;
			RO = MathF.Pow((1 - ((k - 1) / (k + 1) * MathF.Pow(lambda, 2))), k / (k - 1)) * RO_0;

			double M = MathF.Sqrt(((T_0 - T) / T) * (2 / (k - 1)));
			double w = M * c;

			Console.WriteLine($"{T}; {P}; {RO}; {M}; {w}; {isOverSound}");

			results.Add("T", T);
			results.Add("P", P);
			results.Add("RO", RO);
			results.Add("M", M);
			results.Add("w", w);

			return results;
		}
	}

	namespace SplineAPI
	{
		public class HermiteSpline
		{
			private readonly double[] _x;
			private readonly double[] _y;
			private readonly double[] _derivatives;
			private readonly HermiteSegment[] _segments;

			private struct HermiteSegment
			{
				public double x0;
				public double x1;
				public double y0;
				public double y1;
				public double m0;
				public double m1;
			}

			public HermiteSpline(double[] x, double[] y, double[] derivatives)
			{
				if (x == null || y == null || derivatives == null)
					throw new ArgumentNullException("Input arrays cannot be null");

				if (x.Length != y.Length || x.Length != derivatives.Length)
					throw new ArgumentException("All input arrays must have the same length");

				if (x.Length < 2)
					throw new ArgumentException("At least two points are required");

				for (int i = 0; i < x.Length - 1; i++)
				{
					if (x[i] >= x[i + 1])
						throw new ArgumentException("x-coordinates must be strictly increasing");
				}

				_x = (double[])x.Clone();
				_y = (double[])y.Clone();
				_derivatives = (double[])derivatives.Clone();
				Array.Sort(_x, _y);
				Array.Sort(_x, _derivatives);

				_segments = new HermiteSegment[_x.Length - 1];
				for (int i = 0; i < _segments.Length; i++)
				{
					_segments[i] = new HermiteSegment
					{
						x0 = _x[i],
						x1 = _x[i + 1],
						y0 = _y[i],
						y1 = _y[i + 1],
						m0 = _derivatives[i],
						m1 = _derivatives[i + 1]
					};
				}
			}

			private double Evaluate(double t)
			{
				if (t < _x[0] || t > _x[_x.Length - 1])
					return double.NaN;

				int idx = Array.BinarySearch(_x, t);
				if (idx < 0) idx = ~idx - 1;
				if (idx >= _segments.Length) idx = _segments.Length - 1;

				var seg = _segments[idx];
				return HermiteInterpolation(t, seg);
			}

			private double HermiteInterpolation(double x, HermiteSegment seg)
			{
				double t = (x - seg.x0) / (seg.x1 - seg.x0);
				double t2 = t * t;
				double t3 = t2 * t;

				double h00 = 2 * t3 - 3 * t2 + 1;
				double h10 = t3 - 2 * t2 + t;
				double h01 = -2 * t3 + 3 * t2;
				double h11 = t3 - t2;

				return h00 * seg.y0 +
					   h10 * (seg.x1 - seg.x0) * seg.m0 +
					   h01 * seg.y1 +
					   h11 * (seg.x1 - seg.x0) * seg.m1;
			}

			public static Func<double, double> Interpolate(double[] x, double[] y, double[] derivatives)
			{
				var spline = new HermiteSpline(x, y, derivatives);
				return spline.Evaluate;
			}
		}

		public class CubicSpline
		{
			private readonly double[] _x;
			private readonly double[] _y;
			private readonly SplineSegment[] _segments;

			private struct SplineSegment
			{
				public double x0;
				public double x1;
				public double a;
				public double b;
				public double c;
				public double d;
			}

			private CubicSpline(double[] x, double[] y, double leftDerivative, double rightDerivative)
			{
				_x = x;
				_y = y;
				int n = x.Length;
				double[] h = new double[n - 1];

				for (int i = 0; i < n - 1; i++)
				{
					h[i] = x[i + 1] - x[i];
					if (h[i] <= 0)
						throw new ArgumentException("x must be strictly increasing");
				}

				double[] diag = new double[n]; // Главная диагональ
				double[] sub = new double[n];  // Нижняя диагональ
				double[] sup = new double[n];  // Верхняя диагональ
				double[] rhs = new double[n];  // Правая часть

				// Левое граничное условие (производная)
				diag[0] = 2 * h[0];
				sup[0] = h[0];
				rhs[0] = 6 * ((y[1] - y[0]) / h[0] - leftDerivative);

				// Правое граничное условие (производная)
				diag[n - 1] = 2 * h[n - 2];
				sub[n - 1] = h[n - 2];
				rhs[n - 1] = 6 * (rightDerivative - (y[n - 1] - y[n - 2]) / h[n - 2]);

				// Внутренние точки
				for (int i = 1; i < n - 1; i++)
				{
					double sum = h[i - 1] + h[i];
					sub[i] = h[i - 1];
					diag[i] = 2 * sum;
					sup[i] = h[i];
					rhs[i] = 6 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
				}

				double[] M = SolveTridiagonal(sub, diag, sup, rhs);

				_segments = new SplineSegment[n - 1];
				for (int i = 0; i < n - 1; i++)
				{
					double dx = h[i];
					_segments[i] = new SplineSegment
					{
						x0 = x[i],
						x1 = x[i + 1],
						a = y[i],
						b = (y[i + 1] - y[i]) / dx - (M[i + 1] + 2 * M[i]) * dx / 6,
						c = M[i] / 2,
						d = (M[i + 1] - M[i]) / (6 * dx)
					};
				}
			}

			private double[] SolveTridiagonal(double[] a, double[] d, double[] c, double[] b)
			{
				int n = d.Length;
				double[] x = new double[n];

				// Прямой ход
				for (int i = 1; i < n; i++)
				{
					double m = a[i] / d[i - 1];
					d[i] -= m * c[i - 1];
					b[i] -= m * b[i - 1];
				}

				// Обратный ход
				x[n - 1] = b[n - 1] / d[n - 1];
				for (int i = n - 2; i >= 0; i--)
					x[i] = (b[i] - c[i] * x[i + 1]) / d[i];

				return x;
			}

			private double Evaluate(double t)
			{
				if (t < _x[0] || t > _x[_x.Length - 1])
					return double.NaN;

				int idx = Array.BinarySearch(_x, t);
				if (idx < 0)
					idx = ~idx - 1;
				else if (idx == _x.Length - 1)
					idx--;

				var s = _segments[idx];
				double dx = t - s.x0;
				return s.a + s.b * dx + s.c * dx * dx + s.d * dx * dx * dx;
			}

			public static Func<double, double> Interpolate(
				double[] x, double[] y,
				double leftDerivative, double rightDerivative)
			{
				if (x == null || y == null)
					throw new ArgumentNullException();
				if (x.Length != y.Length)
					throw new ArgumentException("Arrays must have the same length");
				if (x.Length < 2)
					throw new ArgumentException("At least two points required");

				double[] xSorted = (double[])x.Clone();
				double[] ySorted = (double[])y.Clone();
				Array.Sort(xSorted, ySorted);

				var spline = new CubicSpline(xSorted, ySorted, leftDerivative, rightDerivative);
				return spline.Evaluate;
			}
		}
	}
}
