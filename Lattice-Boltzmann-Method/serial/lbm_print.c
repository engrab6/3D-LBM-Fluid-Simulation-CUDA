#include <stdio.h>
static const int lbm_velocities[27][3] = {{-1, -1, -1}, { 0, -1, -1}, {+1, -1, -1},
                                          {-1,  0, -1}, { 0,  0, -1}, {+1,  0, -1},
                                          {-1,  1, -1}, { 0,  1, -1}, {+1,  1, -1},

                                          {-1, -1,  0}, { 0, -1,  0}, {+1, -1,  0},
                                          {-1,  0,  0}, { 0,  0,  0}, {+1,  0,  0},
                                          {-1,  1,  0}, { 0,  1,  0}, {+1,  1,  0},

                                          {-1, -1, +1}, { 0, -1, +1}, {+1, -1, +1},
                                          {-1,  0, +1}, { 0,  0, +1}, {+1,  0, +1},
                                          {-1,  1, +1}, { 0,  1, +1}, {+1,  1, +1}};

#define lbm_wc (8./27.)
#define lbm_wf (2./27.)
#define lbm_we (1./54.)
#define lbm_wv (1./216.)

static double lbm_weights[27] = {lbm_wv, lbm_we, lbm_wv,
                                 lbm_we, lbm_wf, lbm_we,
                                 lbm_wv, lbm_we, lbm_wv,

                                 lbm_we, lbm_wf, lbm_we,
                                 lbm_wf, lbm_wc, lbm_wf,
                                 lbm_we, lbm_wf, lbm_we,

                                 lbm_wv, lbm_we, lbm_wv,
                                 lbm_we, lbm_wf, lbm_we,
                                 lbm_wv, lbm_we, lbm_wv};

int main(int argc, char **argv) {
    for (int i = 0; i < 27; i++) {
        if (lbm_velocities[i][0] != 0) {
            printf("x[%d] = %d\n",i, lbm_velocities[i][0]); 
        }
    }
    for (int i = 0; i < 27; i++) {
        if (lbm_velocities[i][1] != 0) {
            printf("y[%d] = %d\n",i, lbm_velocities[i][1]); 
        }
    }
    for (int i = 0; i < 27; i++) {
        if (lbm_velocities[i][2] != 0) {
            printf("z[%d] = %d\n",i, lbm_velocities[i][2]); 
        }
    }
    for (int i = 0; i < 27; i++) {
        // if (lbm_velocities[i][2] != 0) {
            printf("w[%d] = %f\n",i, lbm_weights[i]); 
    }     

    for (int i = 0; i < 27; i++) {
        printf("double u_alpha_%d = ", i);
        // if (lbm_velocities[i][0] != 0) {
            // printf("%d * u + ", lbm_velocities[i][0]);
        // }
        if (lbm_velocities[i][0] == 1) {
            printf("u + "); 
        } else if (lbm_velocities[i][0] == -1) {
            printf("-u + ");
        }
        // if (lbm_velocities[i][1] != 0) {
        //     printf("%d * v + ", lbm_velocities[i][1]);
        // }
        if (lbm_velocities[i][1] == 1) {
            printf("v + " );
        } else if (lbm_velocities[i][1] == -1) {
            printf("(-v) + ");
        }
        if (lbm_velocities[i][2] == 1) {
            printf("w;\n");
        } else if (lbm_velocities[i][2] == -1) {
            printf("(-w);\n");
        } else {
            printf(";\n");
        }
        // printf("double u_alpha = lbm_velocities[dir][0] * u + lbm_velocities[dir][1] * 
        // v + lbm_velocities[dir][2] * w;\n");
    }
    double u_alpha_0 = -(u + v + w);
    double u_alpha_1 = -(v + w);
    double u_alpha_2 = u -(v + w);
    double u_alpha_3 = -(u + w);
    double u_alpha_4 = -w;
    double u_alpha_5 = u - w;
    double u_alpha_6 = -u + v - w;
    double u_alpha_7 = v - w;
    double u_alpha_8 = u + v - w;
    double u_alpha_9 = -(u + v);
    double u_alpha_10 = -v;
    double u_alpha_11 = u - v;
    double u_alpha_12 = -u;
    double u_alpha_13 = 0;
    double u_alpha_14 = u;
    double u_alpha_15 = -u + v;
    double u_alpha_16 = v;
    double u_alpha_17 = u + v;
    double u_alpha_18 = -(u + v) + w;
    double u_alpha_19 = -v + w;
    double u_alpha_20 = u - v + w;
    double u_alpha_21 = -u + w;
    double u_alpha_22 = w;
    double u_alpha_23 = u + w;
    double u_alpha_24 = -u + v + w;
    double u_alpha_25 = v + w;
    double u_alpha_26 = u + v + w;

    return 0;   
}