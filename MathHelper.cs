using System.Numerics;

namespace GAZDIN_CALC_APP
{
	public static class MathHelper
	{
		public static Vector2 CreateVector2(Vector2 start, Vector2 end)
		{
			return end - start;
		}

		public static int IndexOf<T>(T[] array, T element)
		{
			for(int i = 0; i < array.Length; i++)
			{
				if (array[i].Equals(element))
					return i;
			}

			return -1;
		}
	}
}
