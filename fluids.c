// Usage: Drag with the mouse to add smoke to the fluid. This will also move a "rotor" that disturbs
//        the velocity field at the mouse location. Press the indicated keys to change options
//--------------------------------------------------------------------------------------------------

#include <rfftw.h>              //the numerical simulation FFTW library
#include <stdio.h>              //for printing the help text
#include <math.h>               //for various math functions
#include <GLUT/glut.h>            //the GLUT graphics library
#include <string.h>
#include <glui.h>
#include <vector>

#define PI 3.14159265
#define WINDOW_TITLE_PREFIX "Real Time Fluid Flow Simulation Step3"
/** These are the live variables passed into GLUI ***/
int main_window;
GLUI_RadioGroup *radio;
GLUI_RadioGroup *radio2;
GLUI_RadioGroup *radio3;
GLUI_Spinner *min_spinner, *max_spinner;
int minimal = 1, maximal = 256;


//--- SIMULATION PARAMETERS ------------------------------------------------------------------------
const int DIM = 60;            //size of simulation grid
double dt = 0.4;                //simulation time step
float visc = 0.001;                //fluid viscosity
fftw_real *vx, *vy;             //(vx,vy)   = velocity field at the current moment
fftw_real *vx0, *vy0;           //(vx0,vy0) = velocity field at the previous moment
fftw_real *fx, *fy;                //(fx,fy)   = user-controlled simulation forces, steered with the mouse
fftw_real *rho, *rho0;            //smoke density at the current (rho) and previous (rho0) moment
rfftwnd_plan plan_rc, plan_cr;  //simulation domain discretization


//--- VISUALIZATION PARAMETERS ---------------------------------------------------------------------
int winWidth, winHeight;      //size of the graphics window, in pixels
int color_dir = 1;            //use direction color-coding or not
float vec_scale = 1000;            //scaling of hedgehogs
int draw_smoke = 0;           //draw the smoke or not
float clamp_range = 1;
int draw_vecs = 5;//draw the vector field or not
int color_map_dataset = 2;
const int color_rho = 2;//use density datasets
const int color_v = 3;//use fluid velocity magnitude datasets
const int color_f = 1;//use force field magnitude datasets

// streamlines properties
const int streamlines = 5;
float seed_point_x = 300; // set streamline starting points
float seed_point_y = 300;

const int COLOR_BLACKWHITE = 0;   //different types of color mapping: black-and-white, rainbow, banded
const int COLOR_RAINBOW = 1;
const int RED_RAINBOW = 2;
const int COLOR_BANDS = 3;
int scalar_col = 1;           //method for scalar coloring
int frozen = 0;               //toggles on/off the animation



float get_velocity_x(float current_loca_x, float wn, float hn);
float get_velocity_y(float current_loca_y, float wn, float hn);

//------ SIMULATION CODE STARTS HERE -----------------------------------------------------------------

//init_simulation: Initialize simulation data structures as a function of the grid size 'n'.
//                 Although the simulation takes place on a 2D grid, we allocate all data structures as 1D arrays,
//                 for compatibility with the FFTW numerical library.
void init_simulation(int n) {
    int i;
    size_t dim;

    dim = n * 2 * (n / 2 + 1) * sizeof(fftw_real);        //Allocate data structures
    vx = (fftw_real *) malloc(dim);
    vy = (fftw_real *) malloc(dim);
    vx0 = (fftw_real *) malloc(dim);
    vy0 = (fftw_real *) malloc(dim);
    dim = n * n * sizeof(fftw_real);
    fx = (fftw_real *) malloc(dim);
    fy = (fftw_real *) malloc(dim);
    rho = (fftw_real *) malloc(dim);
    rho0 = (fftw_real *) malloc(dim);
    plan_rc = rfftw2d_create_plan(n, n, FFTW_REAL_TO_COMPLEX, FFTW_IN_PLACE);
    plan_cr = rfftw2d_create_plan(n, n, FFTW_COMPLEX_TO_REAL, FFTW_IN_PLACE);

    for (i = 0; i < n * n; i++)                      //Initialize data structures to 0
    { vx[i] = vy[i] = vx0[i] = vy0[i] = fx[i] = fy[i] = rho[i] = rho0[i] = 0.0f; }
}


//FFT: Execute the Fast Fourier Transform on the dataset 'vx'.
//     'dirfection' indicates if we do the direct (1) or inverse (-1) Fourier Transform
void FFT(int direction, void *vx) {
    if (direction == 1) rfftwnd_one_real_to_complex(plan_rc, (fftw_real *) vx, (fftw_complex *) vx);
    else rfftwnd_one_complex_to_real(plan_cr, (fftw_complex *) vx, (fftw_real *) vx);
}

int clamp(float x) { return ((x) >= 0.0 ? ((int) (x)) : (-((int) (1 - (x))))); }

float max(float x, float y) { return x > y ? x : y; }

//solve: Solve (compute) one step of the fluid flow simulation
void solve(int n, fftw_real *vx, fftw_real *vy, fftw_real *vx0, fftw_real *vy0, fftw_real visc, fftw_real dt) {
    fftw_real x, y, x0, y0, f, r, U[2], V[2], s, t;
    int i, j, i0, j0, i1, j1;

    for (i = 0; i < n * n; i++) {
        vx[i] += dt * vx0[i];
        vx0[i] = vx[i];
        vy[i] += dt * vy0[i];
        vy0[i] = vy[i];
    }

    for (x = 0.5f / n, i = 0; i < n; i++, x += 1.0f / n)
        for (y = 0.5f / n, j = 0; j < n; j++, y += 1.0f / n) {
            x0 = n * (x - dt * vx0[i + n * j]) - 0.5f;
            y0 = n * (y - dt * vy0[i + n * j]) - 0.5f;
            i0 = clamp(x0);
            s = x0 - i0;
            i0 = (n + (i0 % n)) % n;
            i1 = (i0 + 1) % n;
            j0 = clamp(y0);
            t = y0 - j0;
            j0 = (n + (j0 % n)) % n;
            j1 = (j0 + 1) % n;
            vx[i + n * j] = (1 - s) * ((1 - t) * vx0[i0 + n * j0] + t * vx0[i0 + n * j1]) +
                            s * ((1 - t) * vx0[i1 + n * j0] + t * vx0[i1 + n * j1]);
            vy[i + n * j] = (1 - s) * ((1 - t) * vy0[i0 + n * j0] + t * vy0[i0 + n * j1]) +
                            s * ((1 - t) * vy0[i1 + n * j0] + t * vy0[i1 + n * j1]);
        }

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            vx0[i + (n + 2) * j] = vx[i + n * j];
            vy0[i + (n + 2) * j] = vy[i + n * j];
        }

    FFT(1, vx0);
    FFT(1, vy0);

    for (i = 0; i <= n; i += 2) {
        x = 0.5f * i;
        for (j = 0; j < n; j++) {
            y = j <= n / 2 ? (fftw_real) j : (fftw_real) j - n;
            r = x * x + y * y;
            if (r == 0.0f) continue;
            f = (fftw_real) exp(-r * dt * visc);
            U[0] = vx0[i + (n + 2) * j];
            V[0] = vy0[i + (n + 2) * j];
            U[1] = vx0[i + 1 + (n + 2) * j];
            V[1] = vy0[i + 1 + (n + 2) * j];

            vx0[i + (n + 2) * j] = f * ((1 - x * x / r) * U[0] - x * y / r * V[0]);
            vx0[i + 1 + (n + 2) * j] = f * ((1 - x * x / r) * U[1] - x * y / r * V[1]);
            vy0[i + (n + 2) * j] = f * (-y * x / r * U[0] + (1 - y * y / r) * V[0]);
            vy0[i + 1 + (n + 2) * j] = f * (-y * x / r * U[1] + (1 - y * y / r) * V[1]);
        }
    }

    FFT(-1, vx0);
    FFT(-1, vy0);

    f = 1.0 / (n * n);
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            vx[i + n * j] = f * vx0[i + (n + 2) * j];
            vy[i + n * j] = f * vy0[i + (n + 2) * j];
        }
}


// diffuse_matter: This function diffuses matter that has been placed in the velocity field. It's almost identical to the
// velocity diffusion step in the function above. The input matter densities are in rho0 and the result is written into rho.
void diffuse_matter(int n, fftw_real *vx, fftw_real *vy, fftw_real *rho, fftw_real *rho0, fftw_real dt) {
    fftw_real x, y, x0, y0, s, t;
    int i, j, i0, j0, i1, j1;

    for (x = 0.5f / n, i = 0; i < n; i++, x += 1.0f / n)
        for (y = 0.5f / n, j = 0; j < n; j++, y += 1.0f / n) {
            x0 = n * (x - dt * vx[i + n * j]) - 0.5f;
            y0 = n * (y - dt * vy[i + n * j]) - 0.5f;
            i0 = clamp(x0);
            s = x0 - i0;
            i0 = (n + (i0 % n)) % n;
            i1 = (i0 + 1) % n;
            j0 = clamp(y0);
            t = y0 - j0;
            j0 = (n + (j0 % n)) % n;
            j1 = (j0 + 1) % n;
            rho[i + n * j] = (1 - s) * ((1 - t) * rho0[i0 + n * j0] + t * rho0[i0 + n * j1]) +
                             s * ((1 - t) * rho0[i1 + n * j0] + t * rho0[i1 + n * j1]);
        }
}

//set_forces: copy user-controlled forces to the force vectors that are sent to the solver.
//            Also dampen forces and matter density to get a stable simulation.
void set_forces(void) {
    int i;
    for (i = 0; i < DIM * DIM; i++) {
        rho0[i] = 0.995 * rho[i];
        fx[i] *= 0.85;
        fy[i] *= 0.85;
        vx0[i] = fx[i];
        vy0[i] = fy[i];
    }
}


//do_one_simulation_step: Do one complete cycle of the simulation:
//      - set_forces:
//      - solve:            read forces from the user
//      - diffuse_matter:   compute a new set of velocities
//      - gluPostRedisplay: draw a new visualization frame
void do_one_simulation_step(void) {
    if (!frozen) {
        set_forces();
        solve(DIM, vx, vy, vx0, vy0, visc, dt);
        diffuse_matter(DIM, vx, vy, rho, rho0, dt);
        glutPostRedisplay();
    }
}


//------ VISUALIZATION CODE STARTS HERE -----------------------------------------------------------------

//rainbow: Implements a color palette, mapping the scalar 'value' to a rainbow color RGB
void rainbow(float value, float *R, float *G, float *B) {
    const float dx = 0.8f;
//    printf("%f\n", dx);

    if (value < 0) value = 0;
    if (value > clamp_range) value = clamp_range;
    value = (6 - 2 * dx) * value + dx;

    *R = max(0.0f, (3 - (float) fabs(value - 4) - (float) fabs(value - 5)) / 2);
    *G = max(0.0f, (4 - (float) fabs(value - 2) - (float) fabs(value - 4)) / 2);
    *B = max(0.0f, (3 - (float) fabs(value - 1) - (float) fabs(value - 2)) / 2);
}

void redwhite(float value, float *R, float *G, float *B) {
    const float dx = 0.2f;
//    printf("%f\n", dx);

    if (value < 0) value = 0;
    if (value > clamp_range) value = clamp_range;
    value = (6 - 2 * dx) * value + dx;

    *R = max(0.0f, (2 - (float) fabs(value - 1) + (float) fabs(value - 2)) / 2);
    *G = max(0.0f, (2 - (float) fabs(value - 1) + (float) fabs(value - 2)) / 2);
    *B = 1;

}

float scaleVelocity(float v, float min, float max) {
    return (v - min) / (max - min);
}

//set_colormap: Sets three different types of colormaps
void set_colormap(float vy) {
    float R, G, B;
//    printf("%d\n", scalar_col);

    if (scalar_col == COLOR_BANDS) {
        float color_clamp_min = (float) minimal / 256;
        float color_clamp_max = (float) maximal / 256;

        if (vy < color_clamp_min) {
            vy = color_clamp_min;
        }
        if (vy > color_clamp_max) {
            vy = color_clamp_max;
        }
    }

    if (scalar_col == COLOR_BLACKWHITE) {
        R = G = B = vy;
    } else if (scalar_col == COLOR_RAINBOW) {
        rainbow(vy, &R, &G, &B);
    } else if (scalar_col == RED_RAINBOW) {
        redwhite(vy, &R, &G, &B);
    } else if (scalar_col == COLOR_BANDS) {
//        const int NLEVELS = 7;
//        printf("before %f\n", vy);
//        vy *= NLEVELS;
//        vy = (int) (vy);
//        vy /= NLEVELS;
//        printf("after %f\n", vy);
        rainbow(vy, &R, &G, &B);
    }
    glColor3f(R, G, B);
}


//direction_to_color: Set the current color by mapping a direction vector (x,y), using
//                    the color mapping method 'method'. If method==1, map the vector direction
//                    using a rainbow colormap. If method==0, simply use the white color
void direction_to_color(float x, float y, int method) {
    float r, g, b, f;
    if (method) {
        f = atan2(y, x) / 3.1415927 + 1;
        r = f;
        if (r > 1) r = 2 - r;
        g = f + .66667;
        if (g > 2) g -= 2;
        if (g > 1) g = 2 - g;
        b = f + 2 * .66667;
        if (b > 2) b -= 2;
        if (b > 1) b = 2 - b;
    } else { r = g = b = 1; }
    glColor3f(r, g, b);
}

void drawLegends() {
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glBegin(GL_QUADS);
    float pTop = 165;
    float pLeft = 715;
    float pRight = 750;

    for (int i = 0; i < 1001; i = i + 1) {
        float vy = 0.001 * i;

        if (scalar_col == COLOR_BANDS) {
            float color_clamp_min = (float) minimal / 256;
            float color_clamp_max = (float) maximal / 256;

            if (vy < color_clamp_min) {
                vy = color_clamp_min;
            }
            if (vy > color_clamp_max) {
                vy = color_clamp_max;
            }
        }

        set_colormap(vy);
        glVertex2f(pLeft, (0.5 * i) + pTop); //(x,y top left)
        glVertex2f(pRight, (0.5 * i) + pTop); //(x,y bottom left)
        glVertex2f(pLeft, (0.5 * (i + 1)) + pTop); //(x,y bottom right)
        glVertex2f(pRight, (0.5 * (i + 1)) + pTop); //(x,y top right)
    }

    glEnd();
    glPopMatrix();

    glBegin(GL_LINES);
    glColor3f(255, 255, 255);
    glVertex2f(pLeft, pTop); //(x,y top left)
    glVertex2f(pRight, pTop); //(x,y bottom left)
    glVertex2f(pLeft, 500 + pTop); //(x,y bottom right)
    glVertex2f(pRight, 500 + pTop); //(x,y top right)
    glVertex2f(pLeft, pTop); //(x,y bottom right)
    glVertex2f(pLeft, 500 + pTop); //(x,y bottom right)
    glVertex2f(pRight, pTop); //(x,y bottom left)
    glVertex2f(pRight, 500 + pTop); //(x,y bottom left)
    glEnd();
}

//visualize: This is the main visualization function
void visualize(void) {
    int i, j, idx, idx0, idx1, idx2, idx3;
    double px0, py0, px1, py1, px2, py2, px3, py3;
    fftw_real wn = (fftw_real) winWidth / (fftw_real) (DIM + 1);   // Grid cell width
    fftw_real hn = (fftw_real) winHeight / (fftw_real) (DIM + 1);  // Grid cell heigh


    if (draw_smoke == 1) {
        if (color_map_dataset == 1) {
            //when dataset is velocity magnitude
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glBegin(GL_TRIANGLES);


            for (j = 0; j < DIM - 1; j++)            //draw smoke
            {
                for (i = 0; i < DIM - 1; i++) {
                    px0 = wn + (fftw_real) i * wn;
                    py0 = hn + (fftw_real) j * hn;

                    idx0 = (j * DIM) + i;

                    px1 = wn + (fftw_real) i * wn;
                    py1 = hn + (fftw_real) (j + 1) * hn;
                    idx1 = ((j + 1) * DIM) + i;

                    px2 = wn + (fftw_real) (i + 1) * wn;
                    py2 = hn + (fftw_real) (j + 1) * hn;
                    idx2 = ((j + 1) * DIM) + (i + 1);

                    px3 = wn + (fftw_real) (i + 1) * wn;
                    py3 = hn + (fftw_real) j * hn;
                    idx3 = (j * DIM) + (i + 1);


                    float f_mag0 = (sqrt(pow(fx[idx0], 2) + pow(fy[idx0], 2))) * 50;
                    float f_mag1 = (sqrt(pow(fx[idx1], 2) + pow(fy[idx1], 2))) * 50;
                    float f_mag2 = (sqrt(pow(fx[idx2], 2) + pow(fy[idx2], 2))) * 50;
                    float f_mag3 = (sqrt(pow(fx[idx3], 2) + pow(fy[idx3], 2))) * 50;


                    set_colormap(f_mag0);
                    glVertex2f(px0, py0);
                    set_colormap(f_mag1);
                    glVertex2f(px1, py1);
                    set_colormap(f_mag2);
                    glVertex2f(px2, py2);

                    set_colormap(f_mag0);
                    glVertex2f(px0, py0);
                    set_colormap(f_mag2);
                    glVertex2f(px2, py2);
                    set_colormap(f_mag3);
                    glVertex2f(px3, py3);
                }
            }
            glEnd();
        } else if (color_map_dataset == 2) {
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glBegin(GL_TRIANGLES);

            for (j = 0; j < DIM - 1; j++)            //draw smoke
            {
                for (i = 0; i < DIM - 1; i++) {
                    px0 = wn + (fftw_real) i * wn;
                    py0 = hn + (fftw_real) j * hn;

                    idx0 = (j * DIM) + i;

                    px1 = wn + (fftw_real) i * wn;
                    py1 = hn + (fftw_real) (j + 1) * hn;
                    idx1 = ((j + 1) * DIM) + i;

                    px2 = wn + (fftw_real) (i + 1) * wn;
                    py2 = hn + (fftw_real) (j + 1) * hn;
                    idx2 = ((j + 1) * DIM) + (i + 1);

                    px3 = wn + (fftw_real) (i + 1) * wn;
                    py3 = hn + (fftw_real) j * hn;
                    idx3 = (j * DIM) + (i + 1);

                    set_colormap(rho[idx0]);
                    glVertex2f(px0, py0);
                    set_colormap(rho[idx1]);
                    glVertex2f(px1, py1);
                    set_colormap(rho[idx2]);
                    glVertex2f(px2, py2);

                    set_colormap(rho[idx0]);
                    glVertex2f(px0, py0);
                    set_colormap(rho[idx2]);
                    glVertex2f(px2, py2);
                    set_colormap(rho[idx3]);
                    glVertex2f(px3, py3);
                }
            }
            glEnd();
        } else if (color_map_dataset == 3) {
            //when dataset is velocity magnitude
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glBegin(GL_TRIANGLES);

            for (j = 0; j < DIM - 1; j++)            //draw smoke
            {
                for (i = 0; i < DIM - 1; i++) {
                    px0 = wn + (fftw_real) i * wn;
                    py0 = hn + (fftw_real) j * hn;

                    idx0 = (j * DIM) + i;

                    px1 = wn + (fftw_real) i * wn;
                    py1 = hn + (fftw_real) (j + 1) * hn;
                    idx1 = ((j + 1) * DIM) + i;

                    px2 = wn + (fftw_real) (i + 1) * wn;
                    py2 = hn + (fftw_real) (j + 1) * hn;
                    idx2 = ((j + 1) * DIM) + (i + 1);

                    px3 = wn + (fftw_real) (i + 1) * wn;
                    py3 = hn + (fftw_real) j * hn;
                    idx3 = (j * DIM) + (i + 1);

                    fftw_real v_mag0 = (sqrt(pow(vx[idx0], 2) + pow(vy[idx0], 2))) * 25;
                    fftw_real v_mag1 = (sqrt(pow(vx[idx1], 2) + pow(vy[idx1], 2))) * 25;
                    fftw_real v_mag2 = (sqrt(pow(vx[idx2], 2) + pow(vy[idx2], 2))) * 25;
                    fftw_real v_mag3 = (sqrt(pow(vx[idx3], 2) + pow(vy[idx3], 2))) * 25;

                    set_colormap(v_mag0);
                    glVertex2f(px0, py0);
                    set_colormap(v_mag1);
                    glVertex2f(px1, py1);
                    set_colormap(v_mag2);
                    glVertex2f(px2, py2);

                    set_colormap(v_mag0);
                    glVertex2f(px0, py0);
                    set_colormap(v_mag2);
                    glVertex2f(px2, py2);
                    set_colormap(v_mag3);
                    glVertex2f(px3, py3);
                }
            }
            glEnd();
        }

    }


    if (draw_vecs == 1) {

////---- define variables for generating arrow head ----//
        double theta = PI / 4;
//        //create a matrix to rotate the velocity vector to get a arrow vector which is pi/4 angular with the velocity
        double matrixA[4] = {cos(theta), -sin(theta), sin(theta), cos(theta)};//counter clockwise rotation matrix
        double matrixB[4] = {cos(theta), sin(theta), -sin(theta), cos(theta)};//clockwise rotation matrix


//        //get start points x_min, y_min and end points x_max,y_max
//        float x_min = 200;
//        float y_min = 200;
//        float x_max = 500;
//        float y_max = 600;
//        int number_points_xdir = 10;// set the number of points in x directions
//        int number_points_ydir = 10; // set the number of points in y directions
////        int number_whole direction =
//        float sample_locations_x [1000];
//        float sample_locations_y [1000];
//
//        float units_distance_x = (x_max - x_min)/number_points_xdir; // get the distance of cells based on subsample range and number of points
//        float units_distance_y = (y_max - y_min)/number_points_ydir;
//        // get the cells between the start points and ending point loat
//        float sample_current_loca_x = x_min;
//        float sample_current_loca_y = y_min;
//        //
//        float sample_vx[1000];
//        float sample_vy[1000];
//
//
//        sample_locations_x[0] = sample_current_loca_x;
//        sample_locations_y[0] = sample_current_loca_y;
//
//
//        for(i=0; i < (number_points_ydir - 1); i++){
//            sample_current_loca_x = x_min;
//            for(j = 0; j <  (number_points_xdir -1); j ++){
//                sample_current_loca_x = sample_current_loca_x + units_distance_x;
//                sample_locations_x[i+j+1] = sample_current_loca_x;
//            }
//        }
//
//
//        for(i=0; i < (number_points_xdir - 1); i++){
//            sample_current_loca_y = y_min;
//            for(j = 0; j <  (number_points_ydir -1); j ++){
//                sample_current_loca_y = sample_current_loca_y + units_distance_y;
//                sample_locations_y[i+j+1] = sample_current_loca_y;
//            }
//        }
//
//        printf("unit %f\n x2 %f \n",units_distance_x,sample_locations_x[1]);
//        printf("x3 %f\n x4 %f \n",sample_locations_x[3],sample_locations_x[4]);
//        printf("x5 %f\n x6 %f \n",sample_locations_x[5],sample_locations_x[6]);
//        printf("x9 %f\n x10 %f \n",sample_locations_x[8],sample_locations_x[9]);
//        printf("x11 %f\n x12 %f \n",sample_locations_x[10],sample_locations_x[11]);
//        printf("y9 %f\n y10 %f \n",sample_locations_y[8],sample_locations_y[9]);
//        printf("y11 %f\n y12 %f \n",sample_locations_y[10],sample_locations_y[11]);

//// ---------- get the velocity of each sample points with bilinear interpolation---- ////

 // for loop to save the velocitu value into array


            //compute the cell corners for the current point
//            int x_corn1 = (int) floor(sample_locations_x[0] / wn - 1);
//            int x_corn2 = (int) ceil(sample_locations_x[0] / wn - 1);
//            int x_corn3 = (int) floor(sample_locations_x[0] / wn - 1);
//            int x_corn4 = (int) ceil(sample_locations_x[0] / wn - 1);
//            int y_corn1 = (int) floor(sample_locations_y[0] / hn - 1);
//            int y_corn2 = (int) ceil(sample_locations_y[0] / hn - 1);
//            int y_corn3 = (int) floor(sample_locations_y[0] / hn - 1);
//            int y_corn4 = (int) ceil(sample_locations_y[0] / hn - 1);
//
//            //compute the index of the corners of a sample point
//            int sample_index_corn1 = y_corn1 * DIM + x_corn1;
//            int sample_index_corn2 = y_corn2 * DIM + x_corn2;
//            int sample_index_corn3 = y_corn3 * DIM + x_corn3;
//            int sample_index_corn4 = y_corn4 * DIM + x_corn4;
//
//            //compute the sample point velocity of x,y directions
//            float samplev_x_corn1 = vx[sample_index_corn1];
//            float samplev_x_corn2 = vx[sample_index_corn2];
//            float samplev_x_corn3 = vx[sample_index_corn3];
//            float samplev_x_corn4 = vx[sample_index_corn4];
//            float samplev_y_corn1 = vy[sample_index_corn1];
//            float samplev_y_corn2 = vy[sample_index_corn2];
//            float samplev_y_corn3 = vy[sample_index_corn3];
//            float samplev_y_corn4 = vy[sample_index_corn4];
//            //compute the distance of sample points to the cell corns
//            float sample_d2corn1x = sample_locations_x[0] - x_corn1 * wn - wn;
//            float sample_d2corn2x = wn - sample_d2corn1x;
//
//            float sample_d2corn1y = sample_locations_y[0] - y_corn1 * hn - hn;
//            float sample_d2corn3y = hn - sample_d2corn1y;
//
//            //compute the interpolartion x,y velocity of sample with space weighted
//            sample_vx[0] = (samplev_x_corn1 * sample_d2corn1x * sample_d2corn1y +
//                               samplev_x_corn2 * sample_d2corn2x * sample_d2corn1y +
//                               samplev_x_corn3 * sample_d2corn1x * sample_d2corn3y +
//                               samplev_x_corn4 * sample_d2corn2x * sample_d2corn3y) / (hn * wn);
//
//            sample_vy[0] = (samplev_y_corn1 * sample_d2corn1x * sample_d2corn1y +
//                               samplev_y_corn2 * sample_d2corn2x * sample_d2corn1y +
//                               samplev_y_corn3 * sample_d2corn1x * sample_d2corn3y +
//                               samplev_y_corn4 * sample_d2corn2x * sample_d2corn3y) / (hn * wn);
//
//
//        glBegin(GL_LINES);//draw velocities
//        glColor3f(1,1,1);
//        glVertex2f(sample_locations_x[0],sample_locations_y[0]);
//        glVertex2f(500,500);
//
//
//        glEnd();
//
////

//
//
//        glBegin(GL_LINES);//draw velocities
//
//        for (i = 0; i < DIM; i += 2)
//            for (j = 0; j < DIM; j += 2) {
//                idx = (j * DIM) + i;
//                direction_to_color(vx[idx], vy[idx], color_dir);
//                float x1 = wn + (fftw_real) i * wn;
//                float y1 = hn + (fftw_real) j * hn;
//                float x2 = (wn + (fftw_real) i * wn) + vec_scale * vx[idx];
//                float y2 = (hn + (fftw_real) j * hn) + vec_scale * vy[idx];
//
//                float arrowhead[2] = {(wn + (fftw_real) i * wn + vec_scale * vx[idx]),
//                                      (hn + (fftw_real) j * hn + vec_scale * vy[idx])};
//                float arrow[2] = {x1 - x2, y1 - y2};
//
//                float rotate_head_l1_x = arrow[0] * matrixA[0] + arrow[1] * matrixA[1];
//                float rotate_head_l1_y = arrow[0] * matrixA[2] + arrow[1] * matrixA[3];
//                float rotate_head_l1[2] = {(rotate_head_l1_x) / 3, (rotate_head_l1_y) / 3};
//
//                float rotate_head_l2_x = arrow[0] * matrixB[0] + arrow[1] * matrixB[1];
//                float rotate_head_l2_y = arrow[0] * matrixB[2] + arrow[1] * matrixB[3];
//                float rotate_head_l2[2] = {(rotate_head_l2_x) / 3, (rotate_head_l2_y) / 3};
//
//                float headvertex1x = (x2) + rotate_head_l1[0];
//                float headvertex1y = (y2) + rotate_head_l1[1];
//                float headvertex2x = (x2) + rotate_head_l2[0];
//                float headvertex2y = (y2) + rotate_head_l2[1];
//
//
//                glVertex2f(wn + (fftw_real) i * wn, hn + (fftw_real) j * hn);
//                glVertex2f((wn + (fftw_real) i * wn) + vec_scale * vx[idx],
//                           (hn + (fftw_real) j * hn) + vec_scale * vy[idx]);
//
//                glVertex2f((wn + (fftw_real) i * wn) + vec_scale * vx[idx],
//                           (hn + (fftw_real) j * hn) + vec_scale * vy[idx]);
//                glVertex2f(headvertex1x, headvertex1y);
////
//                glVertex2f((wn + (fftw_real) i * wn) + vec_scale * vx[idx],
//                           (hn + (fftw_real) j * hn) + vec_scale * vy[idx]);
//                glVertex2f(headvertex2x, headvertex2y);
//
//            }
//        glEnd();

    } else if (draw_vecs == 2) {
        /** Gradient **/

        double theta = PI / 4;
        double matrixA[4] = {cos(theta), -sin(theta), sin(theta), cos(theta)};//counter clockwise rotation matrix
        double matrixB[4] = {cos(theta), sin(theta), -sin(theta), cos(theta)};//clockwise rotation matrix

        glBegin(GL_LINES);//draw velocities

        for (j = 0; j < DIM - 1; j += 2) {
            for (i = 0; i < DIM - 1; i += 2) {
                px0 = wn + (fftw_real) i * wn;
                py0 = hn + (fftw_real) j * hn;

                idx0 = (j * DIM) + i;

                px1 = wn + (fftw_real) i * wn;
                py1 = hn + (fftw_real) (j + 1) * hn;

                idx1 = ((j + 1) * DIM) + i;

                px2 = wn + (fftw_real) (i + 1) * wn;
                py2 = hn + (fftw_real) (j + 1) * hn;

                idx2 = ((j + 1) * DIM) + (i + 1);

                px3 = wn + (fftw_real) (i + 1) * wn;
                py3 = hn + (fftw_real) j * hn;

                idx3 = (j * DIM) + (i + 1);

                fftw_real d_x = 100 * (-rho[idx0] + rho[idx3] - rho[idx1] + rho[idx2]);
                fftw_real d_y = 100 * (rho[idx1] - rho[idx0] + rho[idx2] - rho[idx3]);
                fftw_real threshold = 10;

                if (fabs(d_x) >= fabs(d_y)) {
                    if (d_x >= threshold) {
                        d_y = d_y * threshold / d_x;
                        d_x = threshold;
                    } else if (d_x <= -threshold) {
                        d_y = -d_y * threshold / d_x;
                        d_x = -threshold;
                    }
                } else {
                    if (d_y >= threshold) {
                        d_x = d_x * threshold / d_y;
                        d_y = threshold;
                    } else if (d_y <= -threshold) {
                        d_x = -d_x * threshold / d_y;
                        d_y = -threshold;
                    }
                }

                fftw_real pxm = (px0 + px1 + px2 + px3) / 4;
                fftw_real pym = (py0 + py1 + py2 + py3) / 4;

                set_colormap(rho[idx0]);
                glVertex2f(pxm, pym);

                set_colormap(rho[idx1]);
                glVertex2f(pxm + d_x, pym + d_y);

                fftw_real pxn = pxm + d_x;
                fftw_real pyn = pym + d_y;

                fftw_real arrow[2] = {pxm - pxn, pym - pyn};

                fftw_real rotate_head_l1_x = arrow[0] * matrixA[0] + arrow[1] * matrixA[1];
                fftw_real rotate_head_l1_y = arrow[0] * matrixA[2] + arrow[1] * matrixA[3];
                fftw_real rotate_head_l1[2] = {(rotate_head_l1_x) / 3, (rotate_head_l1_y) / 3};

                fftw_real rotate_head_l2_x = arrow[0] * matrixB[0] + arrow[1] * matrixB[1];
                fftw_real rotate_head_l2_y = arrow[0] * matrixB[2] + arrow[1] * matrixB[3];
                fftw_real rotate_head_l2[2] = {(rotate_head_l2_x) / 3, (rotate_head_l2_y) / 3};

                fftw_real headvertex1x = (pxn) + rotate_head_l1[0];
                fftw_real headvertex1y = (pyn) + rotate_head_l1[1];
                fftw_real headvertex2x = (pxn) + rotate_head_l2[0];
                fftw_real headvertex2y = (pyn) + rotate_head_l2[1];

                glVertex2f(pxm + d_x, pym + d_y);
                glVertex2f(headvertex1x, headvertex1y);

                glVertex2f(pxm + d_x, pym + d_y);
                glVertex2f(headvertex2x, headvertex2y);
            }
        }
        glEnd();

    } else if (draw_vecs == 3) {
        /** Gradient of velocity magnitude **/

        double theta = PI / 4;
        double matrixA[4] = {cos(theta), -sin(theta), sin(theta), cos(theta)};//counter clockwise rotation matrix
        double matrixB[4] = {cos(theta), sin(theta), -sin(theta), cos(theta)};//clockwise rotation matrix

        glBegin(GL_LINES);//draw velocities

        for (j = 0; j < DIM - 1; j += 2) {
            for (i = 0; i < DIM - 1; i += 2) {
                px0 = wn + (fftw_real) i * wn;
                py0 = hn + (fftw_real) j * hn;

                idx0 = (j * DIM) + i;

                px1 = wn + (fftw_real) i * wn;
                py1 = hn + (fftw_real) (j + 1) * hn;

                idx1 = ((j + 1) * DIM) + i;

                px2 = wn + (fftw_real) (i + 1) * wn;
                py2 = hn + (fftw_real) (j + 1) * hn;

                idx2 = ((j + 1) * DIM) + (i + 1);

                px3 = wn + (fftw_real) (i + 1) * wn;
                py3 = hn + (fftw_real) j * hn;

                idx3 = (j * DIM) + (i + 1);

//                fftw_real v_mag0 = (sqrt(pow(vx[idx0], 2) + pow(vy[idx0], 2))) * 25;
//                fftw_real v_mag1 = (sqrt(pow(vx[idx1], 2) + pow(vy[idx1], 2))) * 25;
//                fftw_real v_mag2 = (sqrt(pow(vx[idx2], 2) + pow(vy[idx2], 2))) * 25;
//                fftw_real v_mag3 = (sqrt(pow(vx[idx3], 2) + pow(vy[idx3], 2))) * 25;
//
//                fftw_real d_x = 100 * (-v_mag0 + v_mag1 - v_mag2 + v_mag3);
//                fftw_real d_y = 100 * (v_mag0 - v_mag1 + v_mag2 - v_mag3);

                fftw_real d_x = 2500 * (-vx[idx0] + vx[idx3] - vx[idx1] + vx[idx2]);
                fftw_real d_y = 2500 * (vy[idx1] - vy[idx0] + vy[idx2] - vy[idx3]);
                fftw_real threshold = 10;

                if (fabs(d_x) >= fabs(d_y)) {
                    if (d_x >= threshold) {
                        d_y = d_y * threshold / d_x;
                        d_x = threshold;
                    } else if (d_x <= -threshold) {
                        d_y = -d_y * threshold / d_x;
                        d_x = -threshold;
                    }
                } else {
                    if (d_y >= threshold) {
                        d_x = d_x * threshold / d_y;
                        d_y = threshold;
                    } else if (d_y <= -threshold) {
                        d_x = -d_x * threshold / d_y;
                        d_y = -threshold;
                    }
                }

                fftw_real pxm = (px0 + px1 + px2 + px3) / 4;
                fftw_real pym = (py0 + py1 + py2 + py3) / 4;

                set_colormap(rho[idx0]);
                glVertex2f(pxm, pym);

                set_colormap(rho[idx1]);
                glVertex2f(pxm + d_x, pym + d_y);

                fftw_real pxn = pxm + d_x;
                fftw_real pyn = pym + d_y;

                fftw_real arrow[2] = {pxm - pxn, pym - pyn};

                fftw_real rotate_head_l1_x = arrow[0] * matrixA[0] + arrow[1] * matrixA[1];
                fftw_real rotate_head_l1_y = arrow[0] * matrixA[2] + arrow[1] * matrixA[3];
                fftw_real rotate_head_l1[2] = {(rotate_head_l1_x) / 3, (rotate_head_l1_y) / 3};

                fftw_real rotate_head_l2_x = arrow[0] * matrixB[0] + arrow[1] * matrixB[1];
                fftw_real rotate_head_l2_y = arrow[0] * matrixB[2] + arrow[1] * matrixB[3];
                fftw_real rotate_head_l2[2] = {(rotate_head_l2_x) / 3, (rotate_head_l2_y) / 3};

                fftw_real headvertex1x = (pxn) + rotate_head_l1[0];
                fftw_real headvertex1y = (pyn) + rotate_head_l1[1];
                fftw_real headvertex2x = (pxn) + rotate_head_l2[0];
                fftw_real headvertex2y = (pyn) + rotate_head_l2[1];

                glVertex2f(pxm + d_x, pym + d_y);
                glVertex2f(headvertex1x, headvertex1y);

                glVertex2f(pxm + d_x, pym + d_y);
                glVertex2f(headvertex2x, headvertex2y);
            }
        }
        glEnd();
    } else if (draw_vecs == 5){
        int T = 100;

        float current_vx_list[1000]; // list of the velocity of all streamlines points
        float current_vy_list[1000];


        float current_vx; // current point's velocity
        float current_vy;
        float current_v[2];

        float points_x_list[1000]; // lists of points location
        float points_y_list[1000];

        float current_v_mag[1000];//

        float current_loca_x = seed_point_x;
        float current_loca_y = seed_point_y;
        float cell_diag_distance = sqrt(pow(wn,2) + pow(hn,2));
        float dt = 0.25 * cell_diag_distance;

        float v_x_corn1, v_x_corn2, v_x_corn3, v_x_corn4, v_y_corn1, v_y_corn2,v_y_corn3, v_y_corn4;
//        printf("DIM %i", DIM);

//        wn + (fftw_real) DIM * wn;
//        hn + (fftw_real) DIM * hn;
////
//        float interval1 = (fftw_real) 21 * wn - (fftw_real) 20 * wn;
//        float interval2 = (fftw_real)21;
//        int position = (int) floor( current_loca_x / wn - 1);
//        printf("interval1 %f \n interval 2 %f \n position %f\n", interval1,interval2, position);

     //COMPUTING THE xy index for 4 corners


        points_x_list[0] = current_loca_x;
        points_y_list[0] = current_loca_y;

        for (i = 0; i < T; i++){
            //COMPUTING THE xy for 4 corners
            int x_index_corn1 = (int) floor(current_loca_x / wn -1);
            int x_index_corn2 = (int) ceil(current_loca_x / wn -1);
            int x_index_corn3 = (int) floor(current_loca_x / wn -1);
            int x_index_corn4 = (int) ceil(current_loca_x / wn -1);
            int y_index_corn1 = (int) floor(current_loca_y / hn -1);
            int y_index_corn2 = (int) ceil(current_loca_y / hn -1);
            int y_index_corn3 = (int) floor(current_loca_y / hn -1);
            int y_index_corn4 = (int) ceil(current_loca_y / hn -1);
            //COMPUTING THE index for 4 corners
            int index_corn1 = y_index_corn1 * DIM + x_index_corn1;
            int index_corn2 = y_index_corn2 * DIM + x_index_corn2;
            int index_corn3 = y_index_corn3 * DIM + x_index_corn3;
            int index_corn4 = y_index_corn4 * DIM + x_index_corn4;
            //compute the velocity of x,y directions
            v_x_corn1 = vx[index_corn1];
            v_x_corn2 = vx[index_corn2];
            v_x_corn3 = vx[index_corn3];
            v_x_corn4 = vx[index_corn4];
            v_y_corn1 = vy[index_corn1];
            v_y_corn2 = vy[index_corn2];
            v_y_corn3 = vy[index_corn3];
            v_y_corn4 = vy[index_corn4];

            float d2corn1x = current_loca_x - x_index_corn1 * wn - wn;
            float d2corn2x = wn - d2corn1x;

            float d2corn1y = current_loca_y - y_index_corn1 * hn - hn;
            float d2corn3y = hn - d2corn1y;

            //interpolation with space weighted
            current_vx = (v_x_corn1 * d2corn1x * d2corn1y + v_x_corn2 * d2corn2x * d2corn1y +
                    v_x_corn3 * d2corn1x * d2corn3y + v_x_corn4 * d2corn2x * d2corn3y) / (hn * wn);

            current_vy = (v_y_corn1 * d2corn1x * d2corn1y + v_y_corn2 * d2corn2x * d2corn1y +
                          v_y_corn3 * d2corn1x * d2corn3y + v_y_corn4 * d2corn2x * d2corn3y) / (hn * wn);

            current_vx_list[i] = current_vx;
            current_vy_list[i] = current_vy;
//another interpolation
//            float d2corn1x = current_loca_x - x_index_corn1 * wn - wn;
//            float d2corn2x = wn - d2corn1x;
//            current_vx = (v_x_corn1 * d2corn1x + v_x_corn2 * d2corn2x + v_x_corn3 * d2corn1x + v_x_corn4 * d2corn2x ) / (2 * wn);

//            float d2corn1y = current_loca_y - y_index_corn1 * hn - hn;
//            float d2corn3y = hn - d2corn1y;
//            current_vy = (v_y_corn1 * d2corn1y + v_y_corn3 * d2corn3y + v_y_corn2 * d2corn1y + v_y_corn4 * d2corn3y) / (2 * hn);
            current_v_mag[i] = sqrt(pow(current_vx_list[i],2) + pow(current_vy_list[i],2));

            float norm_current_vx = current_vx / current_v_mag[i];
            float norm_current_vy = current_vy / current_v_mag[i];


            current_loca_x = current_loca_x * (1 + current_vx * dt );
            current_loca_y = current_loca_y * (1 + current_vy * dt );

            points_x_list[i+1] = current_loca_x;
            points_y_list[i+1] = current_loca_y;

        }

        glBegin(GL_LINES);

            for(i = 0; i < T; i++) {

                set_colormap(current_v_mag[i]*50);
                glVertex2f(points_x_list[i], points_y_list[i]);
                glVertex2f(points_x_list[i+1], points_y_list[i+1]);
            }

         glEnd();


    }


    drawLegends();
}
//
float get_velocity_x(float current_loca_x, float wn, float hn){
    float v1x;
    int vx_index_corn1 = (int) floor(current_loca_x / wn -1);
    int vx_index_corn2 = (int) ceil(current_loca_x / wn -1);
    int vx_index_corn3 = (int) floor(current_loca_x / wn -1);
    int vx_index_corn4 = (int) ceil(current_loca_x / wn -1);



    v1x = vx[vx_index_corn1];

    return v1x;
//
//
}

float get_velocity_y(float current_loca_y, float wn, float hn){
    float v1y;
    int vy_index_corn1 = (int) floor(current_loca_y / hn -1);
    int vy_index_corn2 = (int) ceil(current_loca_y / hn -1);
    int vy_index_corn3 = (int) floor(current_loca_y / hn -1);
    int vy_index_corn4 = (int) ceil(current_loca_y / hn -1);

    v1y = vy[vy_index_corn1];

    return v1y;
}


//------ INTERACTION CODE STARTS HERE -----------------------------------------------------------------

//display: Handle window redrawing events. Simply delegates to visualize().
void display(void) {
//    glEnable(GL_DEPTH_TEST);
//    glEnable(GL_COLOR_TABLE);
//    glEnable(GL_SMOOTH);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

//    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//    glEnable(GL_BLEND);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    visualize();
    glFlush();
    glutSwapBuffers();
}

//reshape: Handle window resizing (reshaping) events
void reshape(int w, int h) {
    float panel_length = 215;
    glViewport(0.0f, 0.0f, (GLfloat) w - panel_length, (GLfloat) h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, (GLdouble) w - panel_length, 0.0, (GLdouble) h);
    winWidth = w - panel_length;
    winHeight = h;
}

//keyboard: Handle key presses
void keyboard(unsigned char key, int x, int y) {
    switch (key) {
        case 't':
            dt -= 0.001;
            break;
        case 'T':
            dt += 0.001;
            break;
        case 'c':
            color_dir = 1 - color_dir;
            break;
        case 'S':
            vec_scale *= 1.2;
            break;
        case 's':
            vec_scale *= 0.8;
            break;
        case 'V':
            visc *= 5;
            break;
        case 'v':
            visc *= 0.2;
            break;
        case 'x':
            draw_smoke = 1 - draw_smoke;
            if (draw_smoke == 0) draw_vecs = 1;
            break;
        case 'y':
            draw_vecs = 1 - draw_vecs;
            if (draw_vecs == 0) draw_smoke = 1;
            break;
        case 'm':
            scalar_col++;
            if (scalar_col > COLOR_BANDS) scalar_col = COLOR_BLACKWHITE;
            break;
        case 'a':
            frozen = 1 - frozen;
            break;
        case 'q':
            exit(0);
    }
}


// drag: When the user drags with the mouse, add a force that corresponds to the direction of the mouse
//       cursor movement. Also inject some new matter into the field at the mouse location.
void drag(int mx, int my) {
    int xi, yi, X, Y;
    double dx, dy, len;
    static int lmx = 0, lmy = 0;                //remembers last mouse location

    // Compute the array index that corresponds to the cursor location
    xi = (int) clamp((double) (DIM + 1) * ((double) mx / (double) winWidth));
    yi = (int) clamp((double) (DIM + 1) * ((double) (winHeight - my) / (double) winHeight));

    X = xi;
    Y = yi;
//
    seed_point_x = mx;
    seed_point_y = my;


    if (X > (DIM - 1)) X = DIM - 1;
    if (Y > (DIM - 1)) Y = DIM - 1;
    if (X < 0) X = 0;
    if (Y < 0) Y = 0;

    // Add force at the cursor location
    my = winHeight - my;
    dx = mx - lmx;
    dy = my - lmy;
    len = sqrt(dx * dx + dy * dy);
    if (len != 0.0) {
        dx *= 0.1 / len;
        dy *= 0.1 / len;
    }
    fx[Y * DIM + X] += dx;
    fy[Y * DIM + X] += dy;
    rho[Y * DIM + X] = 10.0f;
    lmx = mx;
    lmy = my;
}


void control_cb(int control) {
    printf("callback: %d\n", control);
//    printf( "             checkbox: %d\n", checkbox->get_int_val() );
    printf("              spinner: %d\n", min_spinner->get_int_val());
    printf("              spinner: %d\n", max_spinner->get_int_val());
    printf("          radio group: %d\n", radio->get_int_val());
//    printf( "                 text: %s\n", edittext->get_text().c_str() );
}

void control_radio(int control) {
    int val = radio->get_int_val();
    scalar_col = val;

    if (val == 0) {
        color_dir = 0;
        min_spinner->disable();
        max_spinner->disable();
    } else if (val == 3) {
        color_dir = 1;
        min_spinner->enable();
        max_spinner->enable();
    } else {
        color_dir = 1;
        min_spinner->disable();
        max_spinner->disable();
    }
}

//main: The main program
int main(int argc, char **argv) {
    printf("Fluid Flow Simulation and Visualization\n");
    printf("=======================================\n");
    printf("Click and drag the mouse to steer the flow!\n");
    printf("T/t:   increase/decrease simulation timestep\n");
    printf("S/s:   increase/decrease hedgehog scaling\n");
    printf("c:     toggle direction coloring on/off\n");
    printf("V/v:   increase decrease fluid viscosity\n");
    printf("x:     toggle drawing matter on/off\n");
    printf("y:     toggle drawing hedgehogs on/off\n");
    printf("m:     toggle thru scalar coloring\n");
    printf("a:     toggle the animation on/off\n");
    printf("q:     quit\n\n");

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(1000, 800);

    main_window = glutCreateWindow(WINDOW_TITLE_PREFIX);

    /*** Create the side subwindow ***/
    GLUI *glui = GLUI_Master.create_glui_subwindow(main_window, GLUI_SUBWINDOW_RIGHT);
    GLUI_Panel *obj_panel = new
            GLUI_Rollout(glui, "Step 2", true);

    /***** Control for colormap *****/
    GLUI_Panel *type_panel = new
            GLUI_Panel(obj_panel, "Colormap");
    radio = new
            GLUI_RadioGroup(type_panel, &scalar_col, 1, control_radio);
    new
            GLUI_RadioButton(radio, "Black and white");
    new
            GLUI_RadioButton(radio, "Rainbow");
    new
            GLUI_RadioButton(radio, "Fantasy");
    new
            GLUI_RadioButton(radio, "Color band");

    min_spinner = new
            GLUI_Spinner(type_panel, "Min:", &minimal, 2, control_cb);
    min_spinner->set_int_limits(2, 256);
//    min_spinner->set_alignment(GLUI_ALIGN_LEFT);
    min_spinner->disable();

    max_spinner = new
            GLUI_Spinner(type_panel, "Max:", &maximal, 3, control_cb);
    max_spinner->set_int_limits(2, 256);
//    max_spinner->set_alignment(GLUI_ALIGN_LEFT);
    max_spinner->disable();

    /***** Control for dataset *****/
    GLUI_Panel *set_panel = new
            GLUI_Panel(obj_panel, "Colormap dataset");
    radio2 = new
            GLUI_RadioGroup(set_panel, &color_map_dataset, 1, control_radio);
    new
            GLUI_RadioButton(radio2, "Black and white");
    new
            GLUI_RadioButton(radio2, "Force field");
    new
            GLUI_RadioButton(radio2, "Fluid density");
    new
            GLUI_RadioButton(radio2, "Fluid velocity magnitude");

    GLUI_Panel *vecs_panel = new
            GLUI_Panel(obj_panel, "Draw vectors");
    radio3 = new
            GLUI_RadioGroup(vecs_panel, &draw_vecs, 1, control_radio);
    new
            GLUI_RadioButton(radio3, "None");
    new
            GLUI_RadioButton(radio3, "Normal");
    new
            GLUI_RadioButton(radio3, "Gradient of fluid density");
    new
            GLUI_RadioButton(radio3, "Gradient of fluid velocity");

    GLUI_Panel *obj_panel2 = new
            GLUI_Rollout(glui, "Step 4", true);

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);

    glutKeyboardFunc(keyboard);
    glutMotionFunc(drag);
    init_simulation(DIM);

//    new GLUI_Checkbox(obj_panel, "color_dir", &scalar_col);
//    (new GLUI_Spinner(obj_panel, "Segments", &segments))->set_int_limits(3, 60);
//    printf("%i\n", wireframe);
//    printf("%i\n", segments);

    glui->set_main_gfx_window(main_window);
    //calls do_one_simulation_step, keyboard, display, drag, reshape
    GLUI_Master.set_glutIdleFunc(do_one_simulation_step);
    glutMainLoop();
    return 0;
}