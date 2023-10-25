using System;

namespace ev2pt2
{
    class Program
    {
        // Programa que incluye los tres métodos numéricos complementarios para analizar dos series de datos.
        // Los datos corresponden a la estatura de dos recién nacidos (Aníbal y María) durante el primer año de vida.
        
        // Variables de Aníbal: y(t) = x1 * ln(t) + x2
        static double[] anibalx = { 50, 55, 60, 61, 65, 67, 69, 70, 72, 73, 74, 76 };
        static double[] anibalt = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

        // Variables de María: y(t) = x1 * cos(t/8) + x2 * e^(t/10)
        static double[] mariax = { 49, 57, 59, 61, 63, 65, 67, 69, 70, 71, 72, 74 };
        static double[] mariat = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };


        static void Main(string[] args)
        {
            Console.WriteLine("Método de mínimos cuadrados");
            Console.WriteLine("funcion de Aníbal");
            var anibalx1 = minimoscuadradosx(anibalx, anibalt);
            var anibalx2 = minimoscuadradosx2(anibalx, anibalt);
            Console.WriteLine("y(t) = {0} * ln(t) + {1}", anibalx1, anibalx2 + "\n");

            Console.WriteLine("--------------------------------------");

            Console.WriteLine("funcion de María");
            var mariax1 = minimoscuadradosx(mariax, mariat);
            var mariax2 = minimoscuadradosx2(mariax, mariat);
            Console.WriteLine("y(t) = {0} * cos(t/8) + {1} * e^(t/10)", mariax1, mariax2 + " \n");

        }

        // Método de mínimos cuadrados para x1
        public static double minimoscuadradosx(double[] x, double[] t)
        {
            //se cambio el numero de incognitas a 2
            int datos = t.Length, incognitas = 2;
            double[,] jacobiana = new double[datos, incognitas];
            double[,] matriz = new double[incognitas, incognitas + 1];

            if (x == anibalx)
            {
                for (int i = 0; i < datos; i++)
                {
                    // PRIMERA COLUMNA: DERIVADA PARCIAL DE X1 EVALUADA EN TIEMPO I
                    jacobiana[i, 0] = Math.Log(t[i]); // Adaptar
                    // SEGUNDA COLUMNA: DERIVADA PARCIAL DE X2 EVALUADA EN TIEMPO I
                    jacobiana[i, 1] = 1; // Adaptar
                }
            }
            else if (x == mariax)
            {
                for (int i = 0; i < datos; i++)
                {
                    // PRIMERA COLUMNA: DERIVADA PARCIAL DE X1 EVALUADA EN TIEMPO I
                    jacobiana[i, 0] = Math.Cos(t[i] / 8); // Adaptar
                    // SEGUNDA COLUMNA: DERIVADA PARCIAL DE X2 EVALUADA EN TIEMPO I
                    jacobiana[i, 1] = Math.Pow(Math.E, t[i] / 10); // Adaptar
                }
            }

            // Encontrar la parte cuadrada de la matriz de Gauss
            // Multiplicar matriz jacobiana * jacobiana T

            for (int i = 0; i < incognitas; i++)
            {
                for (int j = 0; j < incognitas; j++)
                {
                    for (int k = 0; k < datos; k++)
                    {
                        matriz[i, j] += jacobiana[k, i] * jacobiana[k, j];
                    }
                }
            }

            for (int i = 0; i < incognitas; i++)
            {
                for (int j = 0; j < datos; j++)
                {
                    matriz[i, incognitas] += jacobiana[j, i] * x[j];
                }
            }

            // Eliminación gaussiana
            for (int i = 0; i < incognitas; i++)
            {
                double pivote = matriz[i, i];
                if (double.IsNaN(pivote))
                {
                    Console.WriteLine("La matriz no tiene solución.");
                    return 0; // Salir del programa si es NaN
                }
                for (int j = i; j < incognitas + 1; j++)
                {
                    matriz[i, j] /= pivote;
                }
                for (int r = 0; r < incognitas; r++)
                {
                    if (r != i)
                    {
                        double factor = matriz[r, i];
                        for (int k = i; k < incognitas + 1; k++)
                        {
                            matriz[r, k] -= factor * matriz[i, k];
                        }
                    }
                }
            }
            return matriz[0, incognitas];
        }

        // Método de mínimos cuadrados para x2
        public static double minimoscuadradosx2(double[] x, double[] t)
        {
            int datos = t.Length, incognitas = 2;
            double[,] jacobiana = new double[datos, incognitas];
            double[,] matriz = new double[incognitas, incognitas + 1];

            if (x == anibalx)
            {
                for (int i = 0; i < datos; i++)
                {
                    // PRIMERA COLUMNA: DERIVADA PARCIAL DE X1 EVALUADA EN TIEMPO I
                    jacobiana[i, 0] = Math.Log(t[i]); // Adaptar
                    // SEGUNDA COLUMNA: DERIVADA PARCIAL DE X2 EVALUADA EN TIEMPO I
                    jacobiana[i, 1] = 1; // Adaptar
                }
            }
            else if (x == mariax)
            {
                for (int i = 0; i < datos; i++)
                {
                    // PRIMERA COLUMNA: DERIVADA PARCIAL DE X1 EVALUADA EN TIEMPO I
                    jacobiana[i, 0] = Math.Cos(t[i] / 8); // Adaptar
                    // SEGUNDA COLUMNA: DERIVADA PARCIAL DE X2 EVALUADA EN TIEMPO I
                    jacobiana[i, 1] = Math.Pow(Math.E, t[i] / 10); // Adaptar
                }
            }

            // Encontrar la parte cuadrada de la matriz de Gauss
            // Multiplicar matriz jacobiana * jacobiana T

            for (int i = 0; i < incognitas; i++)
            {
                for (int j = 0; j < incognitas; j++)
                {
                    for (int k = 0; k < datos; k++)
                    {
                        matriz[i, j] += jacobiana[k, i] * jacobiana[k, j];
                    }
                }
            }

            for (int i = 0; i < incognitas; i++)
            {
                for (int j = 0; j < datos; j++)
                {
                    matriz[i, incognitas] += jacobiana[j, i] * x[j];
                }
            }

            // Eliminación gaussiana
            for (int i = 0; i < incognitas; i++)
            {
                double pivote = matriz[i, i];
                if (double.IsNaN(pivote))
                {
                    Console.WriteLine("La matriz no tiene solución.");
                    return 0; // Salir del programa si es NaN
                }
                for (int j = i; j < incognitas + 1; j++)
                {
                    matriz[i, j] /= pivote;
                }
                for (int r = 0; r < incognitas; r++)
                {
                    if (r != i)
                    {
                        double factor = matriz[r, i];
                        for (int k = i; k < incognitas + 1; k++)
                        {
                            matriz[r, k] -= factor * matriz[i, k];
                        }
                    }
                }
            }
            return matriz[1, incognitas];
        }


    } 
}
