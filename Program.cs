using static System.Collections.Specialized.BitVector32;

namespace interpolation
{
    internal class Program
    {
        private static double a = 0;
        private static double b = 4;
        static void Main(string[] args)
        {
            Console.WriteLine("введите количество узлов");
            int.TryParse(Console.ReadLine(), out int n);
            var equalNodes = getEqualSectionNodes(n, out double section);
            var nonEqualNodes = getNonEqualSectionNodes(n);
            var equalSectionfunctionValues = getFunctionValues(equalNodes);
            var nonEqualSectionfunctionValues = getFunctionValues(nonEqualNodes);
            var newtonInterpolar = Newton(equalNodes, equalSectionfunctionValues, section);
            var lagrangeInterpolar = Lagrange(nonEqualNodes, nonEqualSectionfunctionValues);
            for (int i = 0; i< equalNodes.Length; i++)
            {
                Console.WriteLine($"{equalNodes[i]};  {equalSectionfunctionValues[i]};  {newtonInterpolar[i]}");
            }
        }

        private static double[] getEqualSectionNodes(int n, out double section)
        {
            section = (b - a) / ((double)n - 1);
            double[] nodes = new double[n];
            nodes[0] = a;
            nodes[n - 1] = b;
            for (int i = 1; i < n - 1; i++)
            {
                nodes[i] = nodes[i-1] + section;
            }
            return nodes;
        }

        private static double[] getNonEqualSectionNodes(int n)
        {
            Random rnd = new Random();
            double[] nodes = new double[n];
            nodes[0] = a;
            nodes[n - 1] = b;

            double totalDistance = b - a;
            double remainingDistance = totalDistance;

            for (int i = 1; i < n - 1; i++)
            {
                double minDistance = remainingDistance / (n - i);
                double maxDistance = remainingDistance * (1 - Math.Pow(1 - 1.0 / n, n - i));
                nodes[i] = nodes[i - 1] + minDistance + rnd.NextDouble() * (maxDistance - minDistance);
                remainingDistance -= (nodes[i] - nodes[i - 1]);
            }

            return nodes;
        }

        private static double[] getFunctionValues(double[] nodes)
        {
            int n = nodes.Length;
            double[] functionValues = new double[n];
            for (int i = 0; i < n; i++)
            {
                functionValues[i] = 1 / 1 + Math.Exp(-nodes[i]);
            }
            return functionValues;
        }

        private static double[] Lagrange(double[] nodes, double[] values)
        {
            int n = nodes.Length;
            double[] interpolarValues = new double[n];
            for (int i = 0; i < n; i++)
            {
                for (int q = 0; q < n; q++)
                {
                    var l = getL(q, nodes[i], nodes);
                    interpolarValues[i] += values[q] * l;
                }
            }
            return interpolarValues;
        }

        private static double getL(int index, double x, double[] nodes)
        {
            double l = 1;
            double q = 1;
            for (int i = 0; i < nodes.Length; i++)
            {
                if (nodes[i] != nodes[index])
                {
                    l *= x - nodes[i];
                    q *= nodes[index] - nodes[i];
                }
            }
            return l/q;
        }

        private static double[] Newton(double[] nodes, double[] values, double h)
        {
            int n = nodes.Length;
            double[] interpolarValues = new double[n];
            for (int i = 0; i < n/2; i++)
            {
                double t = (nodes[i] - nodes[0])/h;
                interpolarValues[i] += values[0] + getDelta(values, 1) * t;
                for (int j = 1; j < n - 3; j++)
                {
                    double factorial = t;
                    double del = 1;
                    for (int q = 1; q < j; q++) {
                        factorial *= t - q;
                        del *= q;
                    }
                    factorial /= del;
                    interpolarValues[i] += getDelta(values, j) * factorial;
                }
            }

            for (int i = n - 1; i >= n/2; i--)
            {
                double t = (nodes[i] - nodes[0]) / h;
                interpolarValues[i] += values[n - 1] + getDelta(values, 1, n - 1) * t;
                for (int j = 0; j < n - 3; j++)
                {
                    double factorial = t;
                    double del = 1;
                    for (int q = 1; q < j; q++)
                    {
                        factorial *= t + q;
                        del *= q;
                    }
                    factorial /= del;
                    interpolarValues[i] += getDelta(values, j, i) * factorial;
                }
            }
            return interpolarValues;
        }

        private static double getDelta(double[] values, int power, int index = 0)
        {
            double delta = values[index + power];
            for (int i = 0; i < power; i++)
            {
                double factorial = power;
                double del = 1;
                for (int q = 1; q < i; q++)
                {
                    factorial *= power - q;
                    del += q;
                }
                factorial /= del;
                if (i%2 == 0)
                    delta -= factorial * values[index + power - i];
                else
                    delta += factorial * values[index + power - i];
            }
            return delta;
        }
    }
}
