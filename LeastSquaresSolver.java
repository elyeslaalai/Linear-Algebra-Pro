package LinearAlgebra;

import java.util.Arrays;

public class LeastSquaresSolver {

    // Method to transpose a matrix
    public static double[][] transpose(double[][] matrix) {
        int rows = matrix.length;
        int cols = matrix[0].length;
        double[][] transposedMatrix = new double[cols][rows];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                transposedMatrix[j][i] = matrix[i][j];
            }
        }
        return transposedMatrix;
    }

    // Method to multiply two matrices
    public static double[][] multiplyMatrices(double[][] matrixA, double[][] matrixB) {
        int rowsA = matrixA.length;
        int colsA = matrixA[0].length;
        int colsB = matrixB[0].length;
        double[][] result = new double[rowsA][colsB];

        for (int i = 0; i < rowsA; i++) {
            for (int j = 0; j < colsB; j++) {
                for (int k = 0; k < colsA; k++) {
                    result[i][j] += matrixA[i][k] * matrixB[k][j];
                }
            }
        }
        return result;
    }

    // Method to multiply a matrix by a vector
    public static double[] multiplyMatrixVector(double[][] matrix, double[] vector) {
        int rows = matrix.length;
        int cols = matrix[0].length;
        double[] result = new double[rows];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i] += matrix[i][j] * vector[j];
            }
        }
        return result;
    }

    // Method to solve a system of linear equations using Gaussian elimination
    public static double[] gaussianElimination(double[][] matrix, double[] vector) {
        int n = matrix.length;

        // Augment the matrix with the vector
        double[][] augmentedMatrix = new double[n][n + 1];
        for (int i = 0; i < n; i++) {
            System.arraycopy(matrix[i], 0, augmentedMatrix[i], 0, n);
            augmentedMatrix[i][n] = vector[i];
        }

        // Apply Gaussian elimination
        for (int i = 0; i < n; i++) {
            // Find the pivot row and swap
            int max = i;
            for (int j = i + 1; j < n; j++) {
                if (Math.abs(augmentedMatrix[j][i]) > Math.abs(augmentedMatrix[max][i])) {
                    max = j;
                }
            }
            double[] temp = augmentedMatrix[i];
            augmentedMatrix[i] = augmentedMatrix[max];
            augmentedMatrix[max] = temp;

            // Make the diagonal element 1
            for (int j = i + 1; j < n + 1; j++) {
                augmentedMatrix[i][j] /= augmentedMatrix[i][i];
            }

            // Eliminate the column entries below the pivot
            for (int j = i + 1; j < n; j++) {
                for (int k = i + 1; k < n + 1; k++) {
                    augmentedMatrix[j][k] -= augmentedMatrix[j][i] * augmentedMatrix[i][k];
                }
            }
        }

        // Back substitution
        double[] solution = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            solution[i] = augmentedMatrix[i][n];
            for (int j = i + 1; j < n; j++) {
                solution[i] -= augmentedMatrix[i][j] * solution[j];
            }
        }
        return solution;
    }

    // Method to solve the least squares problem Ax = b
    public static double[] leastSquares(double[][] A, double[] b) {
        double[][] ATranspose = transpose(A);
        double[][] ATA = multiplyMatrices(ATranspose, A);
        double[] ATb = multiplyMatrixVector(ATranspose, b);
        return gaussianElimination(ATA, ATb);
    }

    // Main method for testing
    public static void main(String[] args) {
        double[][] A = {
            {1, 1, 1},
            {2, 3, 5},
            {4, 0, 5},
            {6, 7, 8}
        };

        double[] b = {6, -4, 27, 20};

        double[] solution = leastSquares(A, b);

        System.out.println("Least Squares Solution:");
        System.out.println(Arrays.toString(solution));
    }
}
