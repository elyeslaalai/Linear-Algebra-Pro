import java.util.Arrays;

public class QRDecomposition {

    // Method to compute the QR Decomposition using the Gram-Schmidt process
    public static void qrDecomposition(double[][] A, double[][] Q, double[][] R) {
        int n = A.length;
        int m = A[0].length;

        double[][] u = new double[n][m];
        double[][] e = new double[n][m];

        // Copy A to u
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                u[i][j] = A[i][j];
            }
        }

        for (int i = 0; i < m; i++) {
            // Calculate the norm of u_i
            double norm = 0;
            for (int j = 0; j < n; j++) {
                norm += u[j][i] * u[j][i];
            }
            norm = Math.sqrt(norm);

            // Normalize u_i to get e_i
            for (int j = 0; j < n; j++) {
                e[j][i] = u[j][i] / norm;
            }

            // Calculate the i-th row of R
            for (int j = i; j < m; j++) {
                R[i][j] = 0;
                for (int k = 0; k < n; k++) {
                    R[i][j] += e[k][i] * A[k][j];
                }
            }

            // Update the remaining columns of u
            for (int j = i + 1; j < m; j++) {
                for (int k = 0; k < n; k++) {
                    u[k][j] -= R[i][j] * e[k][i];
                }
            }
        }

        // Copy e to Q
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                Q[i][j] = e[i][j];
            }
        }
    }

    // Main method for testing
    public static void main(String[] args) {
        double[][] A = {
            {1, 1, 0},
            {1, 0, 1},
            {0, 1, 1}
        };

        int n = A.length;
        int m = A[0].length;

        double[][] Q = new double[n][m];
        double[][] R = new double[m][m];

        qrDecomposition(A, Q, R);

        System.out.println("Matrix Q:");
        for (double[] row : Q) {
            System.out.println(Arrays.toString(row));
        }

        System.out.println("Matrix R:");
        for (double[] row : R) {
            System.out.println(Arrays.toString(row));
        }
    }
}
