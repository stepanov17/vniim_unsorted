
import java.util.ArrayList;
import java.util.TreeMap;

class Result {
    public boolean found = false;
    public double value = Double.NEGATIVE_INFINITY;
    public int indices[] = new int[]{};
}


// Largest consistent subset (Maurice G Cox)
public class LCS {

    // Chi2 critical values, n = 1 to 40
    private static final double CHI2CRIT[] = {
           3.841 // 1
        ,  5.991 // 2
        ,  7.815 // 3
        ,  9.448 // 4
        , 11.070 // 5
        , 12.592 // 6
        , 14.067 // 7
        , 15.507 // 8
        , 16.919 // 9
        , 18.307 // 10
        , 19.675 // 11
        , 21.026 // 12
        , 22.362 // 13
        , 23.685 // 14
        , 24.996 // 15
        , 26.296 // 16
        , 27.587 // 17
        , 28.869 // 18
        , 30.144 // 19
        , 31.410 // 20
        , 32.671 // 21
        , 33.924 // 22
        , 35.172 // 23
        , 36.415 // 24
        , 37.652 // 25
        , 38.885 // 26
        , 40.113 // 27
        , 41.337 // 28
        , 42.557 // 29
        , 43.773 // 30
        , 44.985 // 31
        , 46.194 // 32
        , 47.400 // 33
        , 48.602 // 34
        , 49.802 // 35
        , 50.998 // 36
        , 52.192 // 37
        , 53.384 // 38
        , 54.572 // 39
        , 55.758 // 40
    };

    // return a pair of xRef, uRef
    private double[] getRef(double x[], double u[]) {

        double res[] = new double[2];

        double s = 0.;
        for (double uv: u) { s += 1. / (uv * uv); }
        double uRef2 = 1. / s;
        res[1] = Math.sqrt(uRef2);

        s = 0.;
        for (int i = 0; i < x.length; ++i) {
            s += x[i] * uRef2 / (u[i] * u[i]);
        }
        res[0] = s;

        return res;
    }

    private double getChi2(double x[], double u[]) {

        double ref[] = getRef(x, u);
        double xRef = ref[0];

        int n = x.length;
        double s = 0.;
        for (int i = 0; i < n; ++i) {
            s += ((x[i] - xRef) * (x[i] - xRef)) / (u[i] * u[i]);
        }

        return s;
    }

    private boolean testChi2(double x[], double u[]) {

        double s = getChi2(x, u);
        return (s <= CHI2CRIT[x.length - 1]);
    }

    private void checkData(double x[], double u[]) {

        if (x.length != u.length) {
            throw new RuntimeException("inconsistent x and u dimensions");
        }

        if (x.length < 2) {
            throw new RuntimeException("at least two x values are expected");
        }

        if (x.length > CHI2CRIT.length + 1) {
            throw new RuntimeException("too much data");
        }

        for (double v: u) { if (v <= 0.) { throw new RuntimeException("nonpositive u value: " + v); } }
    }

    public Result getLCS(double x[], double u[], int ind[]) {

        Result res = new Result();

        int l = x.length;
        if (l < 1) { return res; }

        double s = getChi2(x, u);
        if (testChi2(x, u)) {
            res.found = true;
            res.value = s;
            res.indices = ind;
            return res;
        }

        //map length of subres indices length -> res array list
        TreeMap<Integer, ArrayList<Result>> resMap = new TreeMap<>();

        // initialise
        for (int i = 1; i < l; ++i) { resMap.put(i, new ArrayList<>()); }

        for (int i = 0; i < l; ++i) {
            double subx[] = new double[l - 1];
            double subu[] = new double[l - 1];
            int  subind[] = new int   [l - 1];
            int ii = 0;
            for (int j = 0; j < l; ++j) {
                if (j != i) {
                    subx[ii] = x[j];
                    subu[ii] = u[j];
                    subind[ii] = ind[j];
                    ++ii;
                }
            }
            Result r = getLCS(subx, subu, subind);
            if (r.found) {
                int nConsistent = r.indices.length;
                resMap.get(nConsistent).add(r);
            }
        }

        for (int i = l - 1; i > 0; --i) {
            ArrayList<Result> results = resMap.get(i);
            int nRes = results.size();
            if (nRes > 0) {
                int i0 = 0;
                double s0 = results.get(i0).value;
                for (int j = 1; j < nRes; ++j) {
                    double v = results.get(j).value;
                    if (v < s0) {
                        s0 = v;
                        i0 = j;
                    }
                }
                return results.get(i0);
            }
        }

        return res;
    }

    static double[] getSubArr(double arr[], int indices[]) {

        double a[] = new double[indices.length];
        for (int i = 0; i < indices.length; ++i) { a[i] = arr[indices[i]]; }
        return a;
    }

    public static void main(String[] args) {

        double x[] = new double[]{0.0, 1.5, 0.1, 0.1, 2.0, 0.5, 0.1};
        double u[] = new double[]{0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2};

        LCS lcs = new LCS();
        lcs.checkData(x, u);

        int ind[] = new int[x.length];
        for (int i = 0; i < x.length; ++i) { ind[i] = i; }

        Result res = lcs.getLCS(x, u, ind);
        if (res.found) {

            int indices[] = res.indices;
            System.out.print("LCS: " );
            for (int i: indices) { System.out.print(i + " "); }
            System.out.println("");

            double ref[] = lcs.getRef(getSubArr(x, indices), getSubArr(u, indices));
            System.out.println("xRef = " + ref[0] + ", uRef = " + ref[1]);

        } else {
            System.out.println("LCS not found");
        }
    }
}
