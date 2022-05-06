#include "calibration.h"
#include "matrix_algo.h"

using namespace easy3d;


/**
 * TODO: Finish this function for calibrating a camera from the corresponding 3D-2D point pairs.
 *       You may define a few functions for some sub-tasks.
 * @return True on success, otherwise false. On success, the camera parameters are returned by
 */

bool isvalid(const std::vector<Vector3D>& points_3d, const std::vector<Vector2D>& points_2d){

    if (points_3d.size() >= 6 && points_3d.size() == points_2d.size()){
        return true;
    }
    else{
        return false;
    }
}

bool Calibration::calibration(
        const std::vector<Vector3D>& points_3d, /// input: An array of 3D points.
        const std::vector<Vector2D>& points_2d, /// input: An array of 2D image points.
        double& fx, double& fy,    /// output: the focal length (in our slides, we use 'alpha' and 'beta'),
        double& cx, double& cy,    /// output: the principal point (in our slides, we use 'u0' and 'v0'),
        double& skew,              /// output: the skew factor ('-alpha * cot_theta')
        Matrix33& R,               /// output: the 3x3 rotation matrix encoding camera orientation.
        Vector3D& t)               /// output：a 3D vector encoding camera translation.
{
    // TODO: the above code just demonstrates some useful data structures and APIs. Please remove all above code in your
    //       final submission.

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

//
//    for (const auto&point: points_3d){
//        std::cout<<point<<std::endl;
//    }

//    std::cout << "\n[Liangliang]:\n"
//                 "\tThe input parameters of this function are:\n"
//                 "\t\t- points_3d: An array of 3D points (input to this function)\n"
//                 "\t\t- points_2d: An array of 2D image points (input to this function)\n"
//                 "\tThis function must return either 'true' on success or 'false' otherwise. On success, the camera\n"
//                 "\tparameters are returned by the following variables:\n"
//                 "\t\t- fx and fy: the focal lengths (in our slides, we use 'alpha' and 'beta')\n"
//                 "\t\t- cx and cy: the principal point (in our slides, we use 'u0' and 'v0')\n"
//                 "\t\t- skew:      the skew factor ('-alpha * cot_theta')\n"
//                 "\t\t- R:         the 3x3 rotation matrix encoding camera orientation\n"
//                 "\t\t- t:         a 3D vector encoding camera location.\n"
//                 "\tIMPORTANT: don't forget to write your recovered parameters to the above variables." << std::endl;



    // TODO: check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)
    if (!isvalid(points_3d, points_2d)){
        throw std::invalid_argument( "invalid input" );
    }
    std::vector<std::vector<double>> points_3d_and_2d;
    std::vector<double> all_point;

    int j = 0;
    for (const auto & point: points_3d ){
        all_point.push_back(point.x());
        all_point.push_back(point.y());
        all_point.push_back(point.z());
        all_point.push_back(points_2d[j].x());
        all_point.push_back(points_2d[j].y());
        j = j +1;
    }

    //to print the matrix
    for (const auto & point: points_3d_and_2d ){
        std::cout<<point<<std::endl;
    }


    // TODO: construct the P matrix (so P * m = 0).


    Matrix mat = Matrix(points_3d.size(), 5, all_point);

    std::cout << mat << std::endl;

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


    Vector M = Vector(n, 0.0);
    Matrix U = Matrix(m,m,0.0);
    Matrix S = Matrix(m,n, 0.0);
    Matrix V = Matrix(n,n,0.0);

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
    std::cout<<"we here"<<std::endl;
    std::cout << N << "\n" << std::endl;

    std::cout<<std::endl<<"M matric"<< std::endl;
    std::cout<<M<<std::endl;

    for (int i = 0; i < M.size(); i++){
        std::cout<<"test"<<std::endl;
        std::cout<<M[i]<<std::endl;
        if (M[i]> 1){
            throw std::invalid_argument( "M too big input" );
        }
    }

    //add check to see if close to zero



    // TODO: solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
    //   Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point
    //             should be very close to your input images points.



    // TODO: extract intrinsic parameters from M.

    Vector a1 = Vector(M[1][1], M[1][2]);
    Vector a2 = Vector(M[2][1], M[2][2]);
    Vector a3 = Vector(M[3][1], M[3][2]);
    Vector b = Vector(M[1][3], M[2][3], M[3][3]);

    double ro = 0.0;
    double u0 = 0.0;
    double v0 = 0.0;
    double alpha = 0.0;
    double beta = 0.0;
    double costheta = 0.0;
    double sintheta = 0.0; // computation of sintheta
    ro = 1.0/(a3.norm());
    u0 = pow(ro,2)*dot(a1,a3);
    v0 = pow(ro,2)*dot(a2,a3);
    costheta = -(dot(cross(a1,a3),(cross(a2,a3)))/(dot(norm(cross(a1,a3))),norm(cross(a2,a3))));
    double theta = acos(costheta);
    double sintheta = sin(theta);
    alpha = pow(ro,2)*norm(cross(a1,a3))*sintheta;
    beta = pow(ro,2)*norm(cross(a1,a3))*sintheta;

    Matrix K = (3,3,0);
    K[1][1]= alpha;
    K[1][2]= -alpha(costheta*sintheta);
    K[1][3]= u0;
    K[2][1]= 0;
    K[2][2]= beta/sintheta;
    K[2][3]= v0;
    K[3][1]= 0;
    K[3][2]= 0;
    K[3][3]= 1;

    // TODO: extract extrinsic parameters from M.

    Vector2D r1 = (cross(a2,a3))/(norm(cross(a2,a1)));
    Vector2D r3 = ro*(a3);
    Vector2D r2 = cross(r3,r1);

    Vector3D transl = ro* mult(inverse(K),b);




    std::cout << "\n\tTODO: After you implement this function, please return 'true' - this will trigger the viewer to\n"
                 "\t\tupdate the rendering using your recovered camera parameters. This can help you to visually check\n"
                 "\t\tif your calibration is successful or not.\n\n" << std::flush;
    return false;
}










