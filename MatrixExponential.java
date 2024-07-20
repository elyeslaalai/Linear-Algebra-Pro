package LinearAlgebra;

import java.util.Arrays;

public class MatrixExponential {

    // Method to compute the matrix exponential e^(A * t) using the series expansion
    public static double[][] matrixExponential(double[][] matrix, double t) {
        int n = matrix.length;
        double[][] result = new double[n][n];
        double[][] term = new double[n][n];
        double[][] identity = new double[n][n];

        // Initialize result as identity matrix
        for (int i = 0; i < n; i++) {
            result[i][i] = 1.0;
            identity[i][i] = 1.0;
        }

        // Initialize term as identity matrix
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                term[i][j] = identity[i][j];
            }
        }

        int k = 1;
        double factorial = 1.0;
        double[][] scaledMatrix = multiplyMatrixByScalar(matrix, t);

        while (true) {
            term = multiplyMatrices(term, scaledMatrix);
            factorial *= k;
            double[][] newTerm = divideMatrixByScalar(term, factorial);

            if (isMatrixZero(newTerm)) {
                break;
            }

            result = addMatrices(result, newTerm);
            k++;
        }

        return result;
    }

    // Method to multiply two matrices
    private static double[][] multiplyMatrices(double[][] matrixA, double[][] matrixB) {
        int n = matrixA.length;
        double[][] result = new double[n][n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    result[i][j] += matrixA[i][k] * matrixB[k][j];
                }
            }
        }

        return result;
    }

    // Method to add two matrices
    private static double[][] addMatrices(double[][] matrixA, double[][] matrixB) {
        int n = matrixA.length;
        double[][] result = new double[n][n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result[i][j] = matrixA[i][j] + matrixB[i][j];
            }
        }

        return result;
    }

    // Method to divide a matrix by a scalar
    private static double[][] divideMatrixByScalar(double[][] matrix, double scalar) {
        int n = matrix.length;
        double[][] result = new double[n][n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result[i][j] = matrix[i][j] / scalar;
            }
        }

        return result;
    }

    // Method to multiply a matrix by a scalar
    private static double[][] multiplyMatrixByScalar(double[][] matrix, double scalar) {
        int n = matrix.length;
        double[][] result = new double[n][n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result[i][j] = matrix[i][j] * scalar;
            }
        }

        return result;
    }

    // Method to check if a matrix is effectively zero
    private static boolean isMatrixZero(double[][] matrix) {
        int n = matrix.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (Math.abs(matrix[i][j]) > 1e-10) {
                    return false;
                }
            }
        }
        return true;
    }

    // Main method for testing
    public static void main(String[] args) {
        double[][] matrix = {
            {0, 1},
            {-1, 0}
        };

        double t = 1.0; // Exponent value
        double[][] result = matrixExponential(matrix, t);

        System.out.println("Matrix Exponential e^(A * t):");
        for (double[] row : result) {
            System.out.println(Arrays.toString(row));
        }
    }
}
