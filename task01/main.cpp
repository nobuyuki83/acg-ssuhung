#include <cstdio>
#include <cstdlib>
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>
#include <filesystem>
#include <Eigen/Dense>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "delfem2/glfw/viewer2.h"

Eigen::Matrix<double,4,4,Eigen::RowMajor> GetHomographicTransformation(
    const double c1[4][2])
{
  // const double c0[4][2] = {
  //     {0.0,0.0},
  //     {0.0,1.0},
  //     {1.0, 1.0},
  //     {1.0, 0.0} };
  const double c0[4][2] = {
      {-0.5,-0.5},
      {+0.5,-0.5},
      {+0.5,+0.5},
      {-0.5,+0.5} };
  

  Eigen::Matrix<double,4,4,Eigen::RowMajor> m;
  // write some code to compute the 4x4 Homographic transformation matrix `m`;
  // `m` should transfer :
  // (c0[0][0],c0[][1],z) -> (c1[0][0],c1[0][1],z)
  // (c0[1][0],c0[][1],z) -> (c1[1][0],c1[1][1],z)
  // (c0[2][0],c0[][1],z) -> (c1[2][0],c1[2][1],z)
  // (c0[3][0],c0[][1],z) -> (c1[3][0],c1[3][1],z)

  //* Pre-solved linear equation method
  // double a2 = 0.5 * (c1[0][0] + c1[2][0]);
  // double a0 = c1[2][0] + c1[1][0] - 2 * a2;
  // double a1 = 2 * c1[2][0] - 2 * a2 - a0;
  // double a5 = 0.5 * (c1[0][1] + c1[2][1]);
  // double a3 = c1[2][1] + c1[1][1] - 2 * a5;
  // double a4 = 2 * c1[2][1] - 2 * a5 - a3;
  // double a6 = 0, a7 = 0;

  // m <<
  //   a0, a1, a2, 0,
  //   a3, a4, a5, 0,
  //   a6, a7, 1, 0,
  //   0, 0, 0, 1;

  // Eigen::Matrix<double,3,3,Eigen::RowMajor> m_a;
  // m_a <<
  //   a0, a1, a2,
  //   a3, a4, a5,
  //   a6, a7, 1;
  // Eigen::Matrix<double,3,4,Eigen::RowMajor> s;
  // s <<
  //   c0[0][0], c0[1][0], c0[2][0], c0[3][0],
  //   c0[0][1], c0[1][1], c0[2][1], c0[3][1],
  //   1, 1, 1, 1;
  // Eigen::Matrix<double,3,4,Eigen::RowMajor> d;
  // d <<
  //   c1[0][0], c1[1][0], c1[2][0], c1[3][0],
  //   c1[0][1], c1[1][1], c1[2][1], c1[3][1],
  //   1, 1, 1, 1;

  // std::cout << "m_a * s = \n" << m_a * s << std::endl << "d = \n" << d << std::endl;

  //* Solve Linear Equation Method
  Eigen::Matrix<double, 8, 8> A;
  Eigen::Matrix<double, 8, 1> b;

  for(int i = 0; i < 4; i++){
      double s_x = c0[i][0];
      double s_y = c0[i][1];
      double d_x = c1[i][0];
      double d_y = c1[i][1];

      A.row(2 * i) << s_x, s_y, 1.0, 0.0, 0.0, 0.0, (-d_x)*(s_x), (-d_x)*(s_y);
      A.row(2 * i + 1) << 0.0, 0.0, 0.0, s_x, s_y, 1.0, (-d_y)*(s_x), (-d_y)*(s_y);
      b.row(2 * i) << d_x;
      b.row(2 * i + 1) << d_y;
  }
  // std::cout << "A = \n" << A << std::endl << "b = \n" << b << std::endl;
  Eigen::Matrix<double, 8, 1> h = A.colPivHouseholderQr().solve(b);

  // std::cout << "A * h = " << A * h << std::endl << "b = " << b << std::endl;
  // std::cout << "h = " << std::endl << h << std::endl;
  m << 
    h[0], h[1], h[2], 0,
    h[3], h[4], h[5], 0,
    h[6], h[7], 1, 0,
    0, 0, 0, 1;
  
  // std::cout << "m = " << std::endl << m << std::endl;

  // Eigen::Matrix<double, 4, 1> test;
  // test << -0.5, -0.5, 1, 1;
  // std::cout << "first point is " << m * test << std::endl;
  // test << -0.5, 0.5, 1, 1;
  // test = m * test;
  // std::cout << "second point is " << test << std::endl;
  // test[0] /= test[2];
  // test[1] /= test[2];
  // test[2] = 1;
  // std::cout << "but it should be " << test << std::endl;
  
  //* Least Square Root Method
  // Eigen::Matrix<double, 8, 9> A;
  // Eigen::Matrix<double, 8, 1> b;

  // for(int i = 0; i < 4; i++){
  //     double s_x = c0[i][0];
  //     double s_y = c0[i][1];
  //     double d_x = c1[i][0];
  //     double d_y = c1[i][1];

  //     A.row(2 * i) << s_x, s_y, 1.0, 0.0, 0.0, 0.0, (-d_x)*(s_x), (-d_x)*(s_y), -d_x;
  //     A.row(2 * i + 1) << 0.0, 0.0, 0.0, s_x, s_y, 1.0, (-d_y)*(s_x), (-d_y)*(s_y), -d_y;
  //     b.row(2 * i) << 0;
  //     b.row(2 * i + 1) << 0;
  // }
  // // Eigen::Matrix<double, 9, 1> h = A.colPivHouseholderQr().solve(b);
  // Eigen::Matrix<double, 9, 1> h = A.fullPivLu().kernel();

  // std::cout << "A * h = " << A * h << std::endl << "b = " << b << std::endl;
  // // std::cout << "A = " << std::endl << A << std::endl << "b = " << std::endl << b << std::endl;
  // std::cout << "h = " << h << std::endl;

  // m << 
  //   h[0], h[1], h[2], 0,
  //   h[3], h[4], h[5], 0,
  //   h[6], h[7], h[8], 0,
  //   0, 0, 0, 1;

  // Eigen::Matrix<double, 4, 1> test;
  // test << -0.5, -0.5, 1, 1;
  // std::cout << "first point is " << m * test << std::endl;
  // test << -0.5, 0.5, 1, 1;
  // test = m * test;
  // std::cout << "second point is " << test << std::endl;
  // test[0] /= test[2];
  // test[1] /= test[2];
  // test[2] = 1;
  // std::cout << "but it should be " << test << std::endl;

  return m;
}

int main() {

  std::string path = std::string(SOURCE_DIR) + "/../assets/ada.png";
  std::vector<char> img_data;
  int img_width, img_height, img_channels;
  { // load image data using stb library
    stbi_set_flip_vertically_on_load(true);
    assert( std::filesystem::exists(path.c_str()) );
    unsigned char *img = stbi_load(
        path.c_str(),
        &img_width, &img_height, &img_channels, 0);
    assert(img_width > 0 && img_height > 0);
    std::cout << "image size: " << img_width << " " << img_height << " " << img_channels << std::endl;
    img_data.assign(img,img+img_width*img_height*img_channels);
    stbi_image_free(img);
  }

  if (!glfwInit()) { exit(EXIT_FAILURE); }
  // set OpenGL's version (note: ver. 2.1 is very old, but I chose because it's simple)
  ::glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
  ::glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
  GLFWwindow *window = ::glfwCreateWindow(500, 500, "task01", nullptr, nullptr);
  if (!window) { // exit if failed to create window
    ::glfwTerminate();
    exit(EXIT_FAILURE);
  }
  ::glfwMakeContextCurrent(window); // working on this window below

  // set image to texture calling OpenGL's function
  // if you are interested in OpenGL's texture, look: https://learnopengl.com/Getting-started/Textures
  unsigned int texture;
  glGenTextures(1, &texture);
  glBindTexture(GL_TEXTURE_2D, texture);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
               static_cast<int>(img_width),
               static_cast<int>(img_height),
               0, GL_RGBA, GL_UNSIGNED_BYTE,
               img_data.data());


  ::glClearColor(1, 1, 1, 1);
  ::glEnable(GL_DEPTH_TEST);
  ::glEnable(GL_POLYGON_OFFSET_FILL);
  ::glPolygonOffset(1.1f, 4.0f);
  while (!::glfwWindowShouldClose(window)) {
    double time = glfwGetTime();
    const double corners[4][2] = {
        {-0.5, -0.5},
        {+0.5, -0.5},
        {+0.5 - 0.4*cos(1*time), +0.5 - 0.4*sin(3*time)},
        {-0.5 + 0.4*sin(2*time), +0.5 + 0.4*cos(5*time)} };

    Eigen::Matrix<double,4,4,Eigen::ColMajor> modelview_matrix = GetHomographicTransformation(corners);

    ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // set projection matrix
    ::glMatrixMode(GL_PROJECTION);
    ::glLoadIdentity();

    // draw red points
    ::glDisable(GL_TEXTURE_2D);
    ::glDisable(GL_LIGHTING);
    ::glMatrixMode(GL_MODELVIEW);
    ::glLoadIdentity();
    ::glColor3d(1,0,0);
    ::glPointSize(10);
    ::glBegin(GL_POINTS);
    ::glVertex2dv(corners[0]);
    ::glVertex2dv(corners[1]);
    ::glVertex2dv(corners[2]);
    ::glVertex2dv(corners[3]);
    ::glEnd();

    // set model view matrix
    ::glMatrixMode(GL_MODELVIEW);
    ::glLoadIdentity();
    ::glMultMatrixd(modelview_matrix.data());
    ::glEnable(GL_TEXTURE_2D);
    ::glColor3d(1,1,1);
    ::glBegin(GL_QUADS);
    ::glTexCoord2d(0,0);
    ::glVertex2d(-0.5,-0.5);
    ::glTexCoord2d(1,0);
    ::glVertex2d(+0.5,-0.5);
    ::glTexCoord2d(1,1);
    ::glVertex2d(+0.5,+0.5);
    ::glTexCoord2d(0,1);
    ::glVertex2d(-0.5,+0.5);
    ::glEnd();
    //
    ::glfwSwapBuffers(window);
    ::glfwPollEvents();
  }
  ::glfwDestroyWindow(window);
  ::glfwTerminate();
  exit(EXIT_SUCCESS);
}
