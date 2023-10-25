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
            Console.WriteLine("y(t) = {0} * cos(t/8) + {1} * e^(t/10)", mariax1, mariax2);
            Console.WriteLine("");

            Console.WriteLine("--------------------------------------");
            Console.WriteLine("Ecuaciones no lineales");
            ecuacionesnolineales(anibalx1, anibalx2, mariax1, mariax2);

            Console.WriteLine("--------------------------------------");
            Console.WriteLine("Bisección");
            Console.WriteLine(biseccionanibal(anibalx1, anibalx2));
            Console.WriteLine(biseccionmaria(mariax1, mariax2));

            Console.WriteLine("--------------------------------------");
            Console.WriteLine("Secante");
            Console.WriteLine(secanteAnibal(anibalx1, anibalx2));
            Console.WriteLine(secanteMaria(mariax1, mariax2));
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


        //método para el sistema de ecuaciones no lineales que indique en qué momento los dos bebés tendrán la misma estatura, por primera vez
        static void ecuacionesnolineales(double ax1, double ax2, double mx1, double mx2){
        int ren = 2, col = ren + 1;
        double[,] matriz = new double[ren, col];
        double x = 4, y = 57, criterio = 0.00001;
        double pivote, factor;
        double funcionanibal = 1, funcionmaria = 1;
        while( Math.Abs(funcionanibal) > criterio  || Math.Abs(funcionmaria)>criterio){
            funcionanibal = ax1 * Math.Log(x) + ax2 -y;
            funcionmaria =  mx1 * Math.Cos(x/8) + mx2 * Math.Exp(x/10) -y;
            //PRIMERA FUNCION
            matriz[0,0] = ax1 / x; //ADAPTAR
            matriz[0,1] = -1; //ADAPTAR
            matriz[0,2] = -funcionanibal; //ADAPTAR

            //SEGUNDA FUNCION
            matriz[1,0] = mx1 * (-1/8) * Math.Sin(x/8) + mx2 * (1/10) * Math.Exp(x/10); //ADAPTAR 
            matriz[1,1] = -1;//ADAPTAR
            matriz[1,2] = -funcionmaria;//ADAPTAR
            
            // INICIA GAUSS
            //ESTRUCTURA PARA RECORRER MATRIZ
            
            for (int i = 0; i < ren; i++) //PARA RECORRER RENGLONES
            {
            
                //SELECCIONAR PIVOTE
                pivote = matriz[ i ,  i ];
            
                for (int j = 0; j < col; j++)//PARA RECORRER COLUMNAS
                {
            
                    //OPERACION PARA DIVIDIR TODAS LAS CASILLAS DE ESE RENGLON ENTRE EL PIVOTE
                    matriz[i, j] = matriz[i, j] / pivote; //matriz[i,j] /= pivote
                }
                
                //VOLVER A RECORRER LA MATRIZ PARA HACER CERO (RESTA) EL RESTO DE LA COLUMNA
                for (int r = 0; r < ren; r++) //RECORRE RENGLONES
                {
                if(r!=i) //SI NO ES EL RENGLON DEL PIVOTE ENTONCES SE HACE LA RESTA
                {
                    //SELECCIONAR FACTOR
                    factor = matriz[r , i ];
                    for (int k = 0; k < col; k++) //RECORRER COLUMNAS
                    {
                        matriz[r, k] = matriz[r, k] - factor * matriz [ i , k ]; 
                    }
                }
            }
            }
            //TERMINA GAUSS
            // RECORRER COORDENADAS
            x += matriz[0,2];
            y += matriz[1,2];
        }

        //IMPRIMIR EL VALOR X y Y
        Console.WriteLine("Los bebes tendran la misma estatura en el mes " + x + " con una estatura de " + y);
        }

        public static String biseccionanibal(double x1, double x2)
        {
            Console.WriteLine("Metodo de biseccion para anibal");
            double y_inicial, x_inicial = 0, y_final , x_final = 0, x_intermedia=1, y_intermedia=1, paso = 1, criterio = 0.0001;
            //y (t) = x1 * ln(t) + x2
            do    
            {
            x_inicial = x_final; //AVANZA AL SIGUIENTE INTERVALO
            y_inicial = x1 * Math.Log(x_inicial) + x2 - 60; //ADAPTAR
            x_final = x_inicial + paso;
            y_final = x1 * Math.Log(x_final) + x2-60; //ADAPTAR
            }
            while(y_inicial * y_final >0);
            //BISECCION
            while( Math.Abs(y_intermedia) > criterio)
            {
            //CALCULAR x_intermedia
            x_intermedia = (x_inicial + x_final) / 2;
            //CALUCLAR y_internedia y (t) = x1 * ln(t) + x2
            y_intermedia = x1 * Math.Log(x_intermedia) + x2 - 60; //ADAPTAR

            //COMPARAR CON Y INICIAL O Y FINAL LA Y INTERMEDIA, PARA VALIDAR
            if(y_final * y_intermedia > 0) //intermedia y final tienen el mismo signo
            {
                x_final = x_intermedia;
                y_final = y_intermedia;
            }
            else
            {
                x_inicial = x_intermedia;
                y_inicial = y_intermedia;
            }
            }

            return "Raíz: "+ x_intermedia;
        }

                public static String biseccionmaria(double x1, double x2)
        {
            Console.WriteLine("Metodo de biseccion para maria");
            double y_inicial, x_inicial = 0, y_final , x_final = 0, x_intermedia=1, y_intermedia=1, paso = 1, criterio = 0.0001;
            //y(t) = x1 * cos(t/8) + x2 * e^(t/10)
            do    
            {
            x_inicial = x_final; //AVANZA AL SIGUIENTE INTERVALO
            y_inicial = x1 * Math.Cos(x_inicial/8) + x2 * Math.Exp(x_inicial/10) - 60; //ADAPTAR
            x_final = x_inicial + paso;
            y_final = x1 * Math.Cos(x_final/8) + x2 * Math.Exp(x_final/10) - 60; //ADAPTAR
            }
            while(y_inicial * y_final >0);
            //BISECCION
            while( Math.Abs(y_intermedia) > criterio)
            {
            //CALCULAR x_intermedia
            x_intermedia = (x_inicial + x_final) / 2;
            //CALUCLAR y_internedia y (t) = x1 * ln(t) + x2
            y_intermedia = x1 * Math.Cos(x_intermedia/8) + x2 * Math.Exp(x_intermedia/10) - 60; //ADAPTAR

            //COMPARAR CON Y INICIAL O Y FINAL LA Y INTERMEDIA, PARA VALIDAR
            if(y_final * y_intermedia > 0) //intermedia y final tienen el mismo signo
            {
                x_final = x_intermedia;
                y_final = y_intermedia;
            }
            else
            {
                x_inicial = x_intermedia;
                y_inicial = y_intermedia;
            }
            }

            return "Raíz: "+ x_intermedia;
        }

        //Metodo de la secante
        
                public static String secanteAnibal(double x1, double x2)
            {
            Console.WriteLine("Metodo de secante para anibal");      
            double y_inicial, x_inicial = 0, y_final , x_final = 0 , x_intermedia=1, y_intermedia=1, paso = 1, criterio = 0.0001;
            do    
            {
            x_inicial = x_final; //AVANZA AL SIGUIENTE INTERVALO
            y_inicial = x1 * Math.Log(x_inicial) + x2 - 60; //ADAPTAR
            x_final = x_inicial + paso;
            y_final = x1 * Math.Log(x_final) + x2 - 60; //ADAPTAR
            }
            while(y_inicial * y_final >0);

            //BISECCION
            while( Math.Abs(y_intermedia) > criterio)
            {
            //CALCULAR x_intermedia
            x_intermedia = x_final - ( x_inicial - x_final) * y_final  / (y_inicial - y_final);
            //CALUCLAR y_internedia
            y_intermedia = x1 * Math.Log(x_intermedia) + x2 - 60; //ADAPTAR

            //COMPARAR CON Y INICIAL O Y FINAL LA Y INTERMEDIA, PARA VALIDAR
            if(y_final * y_intermedia > 0) //intermedia y final tienen el mismo signo
            {
                x_final = x_intermedia;
                y_final = y_intermedia;
            }
            else
            {
                x_inicial = x_intermedia;
                y_inicial = y_intermedia;
            }
            }

            return "Raíz: "+ x_intermedia;
        } 

                public static String secanteMaria(double x1, double x2)
            {
            Console.WriteLine("Metodo de secante para maria");      
            double y_inicial, x_inicial = 0, y_final , x_final = 0 , x_intermedia=1, y_intermedia=1, paso = 1, criterio = 0.0001;
            do    
            {
            x_inicial = x_final; //AVANZA AL SIGUIENTE INTERVALO
            y_inicial = x1 * Math.Cos(x_inicial/8) + x2 * Math.Exp(x_inicial/10) - 60; //ADAPTAR
            x_final = x_inicial + paso;
            y_final = x1 * Math.Cos(x_final/8) + x2 * Math.Exp(x_final/10) - 60; //ADAPTAR
            }
            while(y_inicial * y_final >0);

            //BISECCION
            while( Math.Abs(y_intermedia) > criterio)
            {
            //CALCULAR x_intermedia
            x_intermedia = x_final - ( x_inicial - x_final) * y_final  / (y_inicial - y_final);
            //CALUCLAR y_internedia
            y_intermedia = x1 * Math.Cos(x_intermedia/8) + x2 * Math.Exp(x_intermedia/10) - 60; //ADAPTAR

            //COMPARAR CON Y INICIAL O Y FINAL LA Y INTERMEDIA, PARA VALIDAR
            if(y_final * y_intermedia > 0) //intermedia y final tienen el mismo signo
            {
                x_final = x_intermedia;
                y_final = y_intermedia;
            }
            else
            {
                x_inicial = x_intermedia;
                y_inicial = y_intermedia;
            }
            }

            return "Raíz: "+ x_intermedia;
        } 
    }
}
