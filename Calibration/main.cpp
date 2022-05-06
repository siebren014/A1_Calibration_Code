/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "calibration.h"
#include "matrix_algo.h"
#include <easy3d/fileio/resources.h>
#include <easy3d/util/logging.h>

using namespace easy3d;

void VectorTest() {

    Vector vec = Vector(10, 88);
    Vector fleurtje = Vector(vec);
    for (int i = 0; i < 10; i++) {
        fleurtje[i] = (i+1)*10.0;
    }

    Vector panda = Vector(2);
    panda[0] = 5;
    panda[1] = 2;

    Vector a = fleurtje;

    for (int i = 0; i < 10; i++) {
        printf("%e %e %e\n",vec[i],fleurtje[i],a[i]);
    }
    printf("dimension %d\n",vec.dimension());
    printf("size %d\n",vec.size());

    fleurtje.resize(2);
    for (int i = 0; i < fleurtje.size(); i++) {
        printf("%e \n", fleurtje[i]);
    }


    printf("squared length %e\n",fleurtje.length2());
    printf(" length %e (%e)\n",fleurtje.length(), sqrt(500));
    printf(" norm %e \n",fleurtje.norm());
    printf(" squared distance %e \n",fleurtje.distance2(panda));

    printf(" normalized fleurtje\n");
    fleurtje.normalize();
    for (int i = 0; i < fleurtje.size(); i++) {
        printf("%e \n", fleurtje[i]);
    }

    Vector vec1 = Vector(2);
    Vector vec2 = Vector(2);
    Vector vec3 = Vector(2);

    vec1[0] = 1.0;
    vec1[1] = 2.0;
    vec2[0] = 3.0;
    vec2[1] = 4.0;

    //vec3 = vec1 + vec2;
    //vec3 = vec1 - vec2;
    //vec3 = vec1[0] * vec2[0];
    //vec3 = vec1[0] / vec2[0];
    vec3 = 3.3 * vec1;
    for (int i = 0; i < vec3.size(); i++) {
        printf("%e %e %e\n", vec1[i], vec2[i], vec3[i]);
    }

    Vector2D v2da = Vector2D(3,8);
    Vector2D v2db = Vector2D(1,8);
    printf("squared length %e\n",v2da.length2());
    printf(" length %e (%e)\n",v2da.length(), sqrt(73));
    printf(" norm %e \n",v2da.norm());
    printf(" squared distance %e\n",v2da.distance2(v2db));

    Vector3D v3da = Vector3D(3,8, 1);
    Vector3D v3db = Vector3D(1,8, 2);
    printf("squared length %e\n",v3da.length2());
    printf(" length %e (%e)\n",v3da.length(), sqrt(74));
    printf(" norm %e \n",v3da.norm());
    printf(" squared distance %e (%e)\n",v3da.distance2(v3db), (v3db-v3da).length2());

    printf("v2db\n");
    for (int i = 0; i < v2db.size(); i++) {
        printf("%e \n", v2db[i]);
    }
    printf("homogeneous v3db\n");
    v3db = v2db.homogeneous();
    for (int i = 0; i < v3db.size(); i++) {
        printf("%e \n", v3db[i]);
    }

    double x = v2db.x();

    printf("v3db\n");
    for (int i = 0; i < v3db.size(); i++) {
        printf("%e \n", v3db[i]);
    }
    printf("homogeneous v3db\n");
    Vector4D v4db = v3db.homogeneous();
    for (int i = 0; i < v4db.size(); i++) {
        printf("%e \n", v4db[i]);
    }

    // print to screen
    vec = Vector(4,88);
    std::cout << vec << std::endl;

    // print to file
    std::ofstream stream_out;
    std::string fileOut = "c:\\tmp\\vector.dat";
    stream_out.open(fileOut);
    if (stream_out.is_open()) {
        stream_out << vec << std::endl;
        stream_out.close();
    }

    // read from screen
    printf("input vector4 \n");
    std::cin >> vec;
    //std::cin >> vec[0] >> vec[1] >> vec[2];
    std::cout << vec << std::endl;

    // read from file
    std::ifstream stream_in;
    std::string fileIn = "c:\\tmp\\vector.dat";
    stream_in.open(fileIn);
    if (stream_in.is_open()) {
        std::string line;
        // each vector on its own row
        while (getline(stream_in, line)) {
            std::istringstream iss(line);
            iss >> vec;
        }
        stream_in.close();
    }
    std::cout << vec << std::endl;
}

void MatrixTest() {
    Matrix A = Matrix(2,2,0.0);

    // print to screen
    Matrix mat = Matrix(2,2, 0.0);
    std::cout << mat << std::endl;

    // print to file
    std::ofstream stream_out;
    std::string fileOut = "c:\\tmp\\matrix.dat";
    stream_out.open(fileOut);
    if (stream_out.is_open()) {
        stream_out << mat << std::endl;
        stream_out.close();
    }

    // read from screen
    printf("input nrow, ncol, matrix \n");
    std::cin >> mat;
    std::cout << mat << std::endl;

    // read from file (single matrix in file
    std::ifstream stream_in;
    std::string fileIn = "c:\\tmp\\matrix.dat";
    stream_in.open(fileIn);
    if (stream_in.is_open()) {
        stream_in >> mat;
        stream_in.close();
    }
    std::cout << mat << std::endl;

}
int main(int argc, char** argv) {
    // the model file.
    const std::string model_file = resource::directory() + "/data/corner.obj";

    try {

        printf("hallo\n");

        //VectorTest();
        //MatrixTest();

        // read matrix from file
        std::ifstream stream_in;
        std::string fileIn = "C:\\Users\\marie\\Documents\\q4geomatics\\geo1016\\A1_Calibration\\A1_Calibration_Code\\resources\\data\\test_data_1(6_points)-test1.txt";
        Matrix mat = Matrix(1,1,0.0);
        stream_in.open(fileIn);
        if (stream_in.is_open()) {
            stream_in >> mat;
            stream_in.close();
        }
        std::cout << mat << std::endl;

        ////Matrix P = mat.resize(mat.rows(),3); did not work
        //Matrix P = Matrix(mat.rows(), mat.cols()-2, 0.0);
        //for (int i = 0; i < P.rows(); i++) {
        //    for (int j =0; j < P.cols(); j++) {
        //        P[i][j] = mat[i][j];
        //    }
        //}
        //std::cout << P << std::endl;

        int m = mat.rows()*2;
        int n = 12;

        Matrix A = Matrix(m, n, 0.0);
        int ii = 0;
        for (int i = 0; i < mat.rows(); i++) {
            A[ii][0] = mat[i][0];
            A[ii][1] = mat[i][1];
            A[ii][2] = mat[i][2];
            A[ii][3] = 1.0;
            A[ii][4] = 0.0;
            A[ii][5] = 0.0;
            A[ii][6] = 0.0;
            A[ii][7] = 0.0;
            A[ii][8] = -1.0*mat[i][0]*mat[i][3];
            A[ii][9] = -1.0*mat[i][1]*mat[i][3];
            A[ii][10] = -1.0*mat[i][2]*mat[i][3];
            A[ii][11] = -1.0*mat[i][3];
            ii += 1;
            A[ii][0] = 0.0;
            A[ii][1] = 0.0;
            A[ii][2] = 0.0;
            A[ii][3] = 0.0;
            A[ii][4] = mat[i][0];
            A[ii][5] = mat[i][1];
            A[ii][6] = mat[i][2];
            A[ii][7] = 1.0;
            A[ii][8] = -1.0*mat[i][0]*mat[i][4];
            A[ii][9] = -1.0*mat[i][1]*mat[i][4];
            A[ii][10] = -1.0*mat[i][2]*mat[i][4];
            A[ii][11] = -1.0*mat[i][4];
            ii += 1;
        }
        // print to file
        std::ofstream stream_out;
        std::string fileOut = "c:\\tmp\\fullmatrix.dat";
        stream_out.open(fileOut);
        if (stream_out.is_open()) {
            stream_out << A << std::endl;
            stream_out.close();
        }

        Vector M = Vector(n, 0.0);
        Matrix U = Matrix(m,m,0.0);
        Matrix S = Matrix(m,n, 0.0);
        Matrix V = Matrix(n,n,0.0);

        // A input matrix to be decomposed M by N matrix
        // U left side M by M orthogonal matrix.
        // S middle M by N diagonal matrix, with zero elements outside of its main diagonal.
        // V The right side N by N orthogonal matrix V.

        svd_decompose(A, U, S, V);
        for (int i = 0; i < n; i++) {
            M[i] = V[i][n-1];
        }

        std::ofstream stream_out1;
        std::string fileOut1 = "c:\\tmp\\mmatrix.dat";
        stream_out1.open(fileOut1);
        if (stream_out1.is_open()) {
            stream_out1 << M << std::endl;
            stream_out1.close();
        }

        Vector N = Vector(m,5.0);
        N = mult(A,M);
        std::cout << N << "\n" << std::endl;


        //Calibration viewer("Calibration", model_file);

        // Run the viewer
        //viewer.run();
    } catch (const std::runtime_error &e) {
        LOG(ERROR) << "caught a fatal error: " + std::string(e.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
