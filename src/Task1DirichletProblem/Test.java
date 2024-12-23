package Task1DirichletProblem;

import java.util.ArrayList;

import static java.lang.Math.PI;
import static java.lang.Math.cos;

public class Test {
    public static void main(String[] args) {
//        int[] a2 = {1, 15, 7, 9, 3, 13, 5, 11};
//        for (int i = 0; i < a2.length; i++) {
//            System.out.print(cos( (double) a2[i]/(2* a2.length) * PI ) + " ");
//            System.out.println(cos( (double) (2 * i - 1) /(2*a2.length) *PI));
//        }
        constructTau_k();
    }

    public static void constructTau_k(){
        int n_curr = 2;
        ArrayList<Integer> oPrev = new ArrayList<>(){};
        ArrayList<Integer> o = new ArrayList<>();
        ArrayList<Integer> tmp = new ArrayList<>();
        oPrev.add(1);
        o.add(1);

        while(n_curr <= 16) {
            tmp = o;
            o = oPrev;
            oPrev = tmp;
            o.clear();

            for (int i = 0; i < n_curr / 2; i++) {
                o.add(oPrev.get(i));
                o.add(4 * (n_curr / 2) - o.get(2 * i));
            }
            System.out.println(o);
            n_curr *= 2;
        }
    }
}
