package Task1DirichletProblem;

import org.mariuszgromada.math.mxparser.Function;
import org.mariuszgromada.math.mxparser.License;
import dnl.utils.text.table.TextTable;

import java.util.ArrayList;
import java.util.Arrays;

import static java.lang.Math.*;
import static org.mariuszgromada.math.mxparser.mathcollection.MathFunctions.ln;



public class Main6 {
    static int kMax = (int)1e10;


    static int N = 10, M = 10;
    static double xMin = 0, xMax = 1, yMin = 0, yMax = 1;
    static double eps = 1e-3;
    static double[][] U;
    static double[][] UPrev;
    static double[][] UTrue;
    static double hx;
    static double hy;
    static double delta;
    static double Delta;
    static double tau;
    static double[][] F;
    static double startNorm;
    static Function nu;
    static Function f;
    static double[][] LhUTrue;
    static double[][] Lh;
    static double[][] Phi_curr;
    static double[][] w_cap;
    static double[][] w;

    static int pk = 15; // переменная, после которой будет логироваться инфа через каждые 10 итераций
    static double ksi;
    static double rho;
    static double startIncoherence;
    static double norm_Uk_Uprev = 1;
    static double spRPrev = 1;
    static double gamma_1;
    static double gamma_2;
    static double omega;
    static double eta;
    static double cappa_1;
    static double cappa_2;


    static int n = 16;
    static ArrayList<Double> tau_k;

    static ArrayList<String> iterations = new ArrayList<>();
    static ArrayList<String> normIncoherence = new ArrayList<>();
    static ArrayList<String> incoherence =  new ArrayList<>();
    static ArrayList<String> normInaccuracy =  new ArrayList<>();
    static ArrayList<String> inaccuracy =  new ArrayList<>();
    static ArrayList<String> normNeighbour =  new ArrayList<>();
    static ArrayList<String> inaccuracy_k =  new ArrayList<>();
    static ArrayList<String> rho_k = new ArrayList<>();




    public static void main(String[] args) {
        // строчка для библиотеки
        License.iConfirmNonCommercialUse("xSICHx");

        System.out.println("Simple iteration method. Variant 6");
        System.out.println("N = " + N + ", M = " + M);
        System.out.println("eps = " + eps);
        System.out.println("xMin = "+ xMin + ", xMax = "+xMax+", yMin = "+yMin+", yMax = " + yMax);
        System.out.println("p = " + 1 + ", q = " + 1);
        System.out.println("f(x,y) = -12*(x*(y^3) + x^3*y)");
        System.out.println("nu(x,y) = 2*(x^3)*(y^3)");
        System.out.println("u*(x,y) = 2*x^3*y^3");

        f = new Function("f(x,y) = -12*(x*(y^3) + x^3*y)");
        nu = new Function(("f(x,y) = 2*(x^3)*(y^3)"));
        U = new double[N+1][M+1];
        hx = (xMax - xMin) / N;
        hy = (yMax- yMin) / M;
        delta = calc_delta(xMax - xMin, yMax- yMin, N, M);
        Delta = calc_Delta(xMax - xMin, yMax- yMin, N, M);
        omega = 2/sqrt(delta*Delta);
        eta = delta/Delta;
        gamma_1 = delta/(2+2*sqrt(eta));
        gamma_2 = delta/(4*sqrt(eta));
        tau = 2/(gamma_1 + gamma_2);
        ksi = delta/Delta;
        rho = (1-ksi)/(1+ksi);
        cappa_1 = omega/(hx*hx);
        cappa_2 = omega/(hy*hy);
        tau_k = new ArrayList<>();
        constructTau_k();


        UPrev = new double[N+1][M+1];

        // Заполняем границы
        fillBorder();

        // Посчитаем F
        F = new double[N+1][M+1];
        calcF();

        // Посчитаем настояещее U
        UTrue = new double[N+1][M+1];
        for (int i = 0; i <= N; i++) {
            for (int j = 0; j <= M ; j++) {
                UTrue[i][j] = nu.calculate(hx*i + xMin, hy*j+yMin);
            }
        }

        startNorm = calcNorm(U, UTrue);

        // ||F - Au*||
        LhUTrue = new double[N+1][M+1];
        calcLu(LhUTrue, UTrue);
        System.out.println("1. Measure of approximation ||F - Lu*||: " + calcNorm(LhUTrue, F));

        // ||F - Au_0||
        Lh = new double[N+1][M+1];
//        calcLu(Lh, U);
        startIncoherence = calcNorm(Lh, F);
        System.out.println("2. Discrepancy norm for u_0 ||F - Lu_0||: " + startIncoherence);
//        System.out.println("3. m >= " + ln(2/eps)/(2*sqrt(2) * sqrt(sqrt(eta))));
        System.out.println("3. m >= " + ln(1/eps)/ln(1/rho)/2);
        // Начальное Phi
        Phi_curr = new double[N+1][M+1];
//        calc_Phi();
        //TODO Lu

        //Начальное w_cap
        w_cap = new double[N+1][M+1];
        init_w_cap();
//        calc_w_cap();

        //Начальное w
        w = new double[N+1][M+1];
        init_w();
//        calc_w();

        // Выполним итерации пока не будет выполнено условие
        int i = 0;
        while ( (stopCondition() && (i < kMax)) || i % n != 0 ){
            // swap
            double[][] temp = UPrev;
            UPrev = U;
            U = temp;

            doIteration(i);
            double currSpRad = calcSpRad();

            if (i < pk || i % 10 == 0){
                logInfo(i, currSpRad);
            }
            i++;

        }
        logInfo(i, 1);
        // Напечатать все таблички
        printLog();

    }


    static void calc_w(){
        for (int i = 1; i < N; i++) {
            for (int j = 1; j < M; j++) {
                w[i][j] = (cappa_1*w[i+1][j] + cappa_2*w[i][j+1] + w_cap[i][j])/
                        (1 + cappa_1 + cappa_2);
            }
        }

    }


    static void init_w(){
        for (int i = 0; i <= N; i++) {
            for (int j = 0; j <= M; j++) {
                w[i][j] = 0;
            }
        }


    }

    static void calc_w_cap(){
        for (int i = 1; i < N; i++) {
            for (int j = 1; j < M; j++) {
                w_cap[i][j] = (cappa_1*w_cap[i-1][j] + cappa_2*w_cap[i][j-1] + Phi_curr[i][j])/
                        (1 + cappa_1 + cappa_2);
            }
        }

    }


    static void init_w_cap(){
        for (int i = 0; i <= N; i++) {
            for (int j = 0; j <= M; j++) {
                w_cap[i][j] = 0;
            }
        }


    }



    static void calc_Phi(){
        for (int i = 0; i <= N; i++) {
            for (int j = 0; j <= M; j++) {
                Phi_curr[i][j] = -Lh[i][j] + F[i][j];
            }
        }
    }

    static double calc_delta(double lx, double ly, int N, int M){
        double hx = lx/N;
        double hy = ly/M;
        double snx = sin(Math.PI*hx/(2*lx));
        double sny = sin(Math.PI*hy/(2*ly));

        return 4*snx*snx/(hx*hx) + 4*sny*sny/(hy*hy);
    }
    static double calc_Delta(double lx, double ly, int N, int M){
        double hx = lx/N;
        double hy = ly/M;
        double snx = cos(Math.PI*hx/(2*lx));
        double sny = cos(Math.PI*hy/(2*ly));

//        return 4*snx*snx/(hx*hx) + 4*sny*sny/(hy*hy);
        return 4/(hx*hx) + 4/(hy*hy);
    }

    static void doIteration(int index){
        calcLu(Lh, U);
        calc_Phi();
        calc_w_cap();
        calc_w();
        for (int i = 1; i < N; i++) {
            for (int j = 1; j < M; j++) {
//                tau = tau_k.get((index) % n);
                tau = 2/(Delta + delta);
                U[i][j] = UPrev[i][j] + tau*(w[i][j]);
            }
        }
    }

    // Norm
    static double calcNorm(double[][] U1,  double[][] U2){
        double mx = abs( U1[2][2] - U2[2][2]);

        for (int i = 1; i < N; i++) {
            for (int j = 1; j < M ; j++) {
//                double asdsa = abs(U1[i][j] - U2[i][j]);
                mx = Math.max(abs(U1[i][j] - U2[i][j]), mx);
            }
        }
        return mx;
    }
    static boolean stopCondition(){
        double sm = calcNorm(U, UTrue);
        return !(sm/startNorm < eps);
    }

    static void fillBorder(){
        for (int i = 0; i <= N; ++i) {
            UPrev[i][0] = nu.calculate(hx*i + xMin, yMin);
            UPrev[i][M] = nu.calculate(hx*i + xMin, yMax);
            U[i][0] = UPrev[i][0];
            U[i][M] = UPrev[i][M];
        }
        for (int i = 0; i <= M; ++i) {
            UPrev[0][i] = nu.calculate(xMin, hy*i+yMin);
            UPrev[N][i] = nu.calculate(xMax, hy*i+yMin);
            U[0][i] = UPrev[0][i];
            U[N][i] = UPrev[N][i];
        }
    }

    static void calcF(){
        for (int i = 0; i <= N; i++) {
            for (int j = 0; j <= M; j++) {
                F[i][j] = f.calculate(hx*i + xMin, hy*j+yMin);
            }
        }
    }

    static void calcLu(double[][] L, double[][]u){
        for (int i = 1; i < N; i++) {
            for (int j = 1; j < M; j++) {
                L[i][j] =
                        -(
                                (u[i + 1][j] - u[i][j]) / (hx * hx) -
                                        (u[i][j] - u[i-1][j]) / (hx * hx) +
                                        (u[i][j + 1] - u[i][j]) / (hy * hy) -
                                        (u[i][j] - u[i][j - 1]) / (hy * hy)
                        );
            }
        }
    }

    static void logInfo(int i, double spRad){
        iterations.add(String.valueOf(i));
//        System.out.print(i+" ");

//        calcLu(Lh, U);
        double norm = calcNorm(Lh, F);
        normIncoherence.add(String.format("%.2e", norm));
//        System.out.print(norm+" ");

        incoherence.add(String.format("%.2e", norm/startIncoherence));
//        System.out.print(norm/startIncoherence+" ");

        norm = calcNorm(U, UTrue);
        normInaccuracy.add(String.format("%.2e", norm));
//        System.out.print(norm+" ");

        inaccuracy.add(String.format("%.2e", norm/startNorm));
//        System.out.print(norm/startNorm+" ");

        norm = calcNorm(U, UPrev);
        normNeighbour.add(String.format("%.2e", norm));
//        System.out.print(norm+" ");

//        inaccuracy_k.add(String.format("%.2e", norm*rho/(1-rho)));
////        System.out.println(norm*rho/(1-rho));
//
//        if (i == 0)
//            rho_k.add("-");
//        else
//            rho_k.add(String.format("%.4f", spRad));
    }

    public static String[][] to2DArray(ArrayList<ArrayList<String>> list) {
//        // Определяем размеры двумерного массива
//        int rows = list.size();
//        int cols = list.get(0).size();
//
//        // Создаем массив
//        String[][] array = new String[rows][cols];
//
//        // Копируем данные из ArrayList в массив
//        for (int i = 0; i < rows; i++) {
//            for (int j = 0; j < cols; j++) {
//                array[i][j] = list.get(i).get(j);
//            }
//        }
//
//        return array;

        // Определяем размеры исходного списка
        int rows = list.size();
        int cols = list.get(0).size();

        // Создаем транспонированный массив (размеры меняются местами)
        String[][] transposedArray = new String[cols][rows];

        // Копируем данные, меняя местами индексы строк и столбцов
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                transposedArray[j][i] = list.get(i).get(j);
            }
        }

        return transposedArray;
    }


    public static void printUTable(String name, double[][] u){
        System.out.println(name);
        String[] columns = new String[N+2];
        columns[0] = "x/y";
        for (int i = 1; i <= N+1; i++) {
            columns[i] = String.format("%.2f", hx*(i-1) + xMin);
        }

        String[][] tab = new String[N+1][M+2];
        for (int i = 0; i <= M; i++) {
            tab[i][0] = String.format("%.2f", hy*i + yMin);
        }

        for (int i = 0; i <= N; i++) {
            for (int j = 0; j <= M; j++) {
                tab[i][j+1] = String.format("%.2e", u[i][j]);
            }
        }

        TextTable table = new TextTable(columns, tab);
        table.printTable();

    }
    public static void printLog(){
        // Заполнение основной таблицы
        String[][] dataTable = to2DArray(new ArrayList<>(
                Arrays.asList(iterations,
                        normIncoherence,
                        incoherence,
                        normInaccuracy,
                        inaccuracy,
                        normNeighbour
//                        ,
//                        inaccuracy_k,
//                        rho_k
                )
        ));


        // Таблица для вывода результатов
        String[] tableColumns = {
                "k",
                "||F - AU_k",
                "rel.d.",
                "||Uk - u*||",
                "rel.error",
                "||U_k - U_(k-1)||"
//                ,
//                "apost.est.",
//                "sp.rad._k"
        };
        TextTable resultsTable;
        resultsTable = new TextTable(tableColumns, dataTable);
        resultsTable.printTable();
        System.out.println();


        printUTable("Exact solution:", UTrue);
        System.out.println();
        printUTable("Approximate solution:", U);
    }


    public static double calcSpRad(){
        double norm = calcNorm(U, UPrev);
        double spR = norm / norm_Uk_Uprev;
        double res = sqrt(spR*spRPrev);
//        double res = spR;
        norm_Uk_Uprev = norm;
        spRPrev = spR;
        return res;
    }

    public static void constructTau_k(){
        int n_curr = 2;
        ArrayList<Integer> oPrev = new ArrayList<>(){};
        ArrayList<Integer> o = new ArrayList<>();
        ArrayList<Integer> tmp = new ArrayList<>();
        oPrev.add(1);
        o.add(1);

        while(n_curr <= n) {
            tmp = o;
            o = oPrev;
            oPrev = tmp;
            o.clear();

            for (int i = 0; i < n_curr / 2; i++) {
                o.add(oPrev.get(i));
                o.add(4 * (n_curr / 2) - o.get(2 * i));
            }
//            System.out.println(o);
            n_curr *= 2;
        }

        // конструируем tau_k
        for (int i = 0; i < n; i++) {
            double cs = cos((double) o.get(i) /(2*n)*PI);
            double val = 2/(Delta + delta + (Delta - delta) * cs);
            tau_k.add(val);
        }
    }
}
