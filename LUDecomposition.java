import java.util.Arrays;

public class LUDecomposition {

    // Method to compute the LU Decomposition
    public static void luDecomposition(double[][] A, double[][] L, double[][] U) {
        int n = A.length;

        // Initialize L and U
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    L[i][j] = 1.0;
                } else {
                    L[i][j] = 0.0;
                }
                U[i][j] = 0.0;
            }
        }

        // Perform the LU Decomposition
        for (int j = 0; j < n; j++) {
            for (int i = 0; i <= j; i++) {
                double sum = 0.0;
                for (int k = 0; k < i; k++) {
                    sum += L[i][k] * U[k][j];
                }
                U[i][j] = A[i][j] - sum;
            }
            for (int i = j + 1; i < n; i++) {
                double sum = 0.0;
                for (int k = 0; k < j; k++) {
                    sum += L[i][k] * U[k][j];
                }
                L[i][j] = (A[i][j] - sum) / U[j][j];
            }
        }
    }

    // Main method for testing
    public static void main(String[] args) {
        double[][] A = {
            {4, 3},
            {6, 3}
        };

        int n = A.length;

        double[][] L = new double[n][n];
        double[][] U = new double[n][n];

        luDecomposition(A, L, U);

        System.out.println("Matrix L:");
        for (double[] row : L) {
            System.out.println(Arrays.toString(row));
        }

        System.out.println("Matrix U:");
        for (double[] row : U) {
            System.out.println(Arrays.toString(row));
        }
    }
}

/*
Expected output:
Matrix L:
[1.0, 0.0]
[1.5, 1.0]
Matrix U:
[4.0, 3.0]
[0.0, -1.5]
*/
